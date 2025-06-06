#![no_std]

//! This module contains code for determing attitude and heading, as well a sensor
//! fusion with locaiton.
//!
//! To determine attitude, it fuses acceleeration, gyro, and mag data, 3D each.
//!
//!Conventions We define +X to be right, +Y to be forward and
//! +Z to be u. We use the right hand rule for rotations along these axes.
//!
//! Magnetic inclination chart: https://upload.wikimedia.org/wikipedia/commons/d/de/World_Magnetic_Inclination_2015.pdf
//! https://www.magnetic-declination.com/

mod acc;
pub mod attitude;
mod linear_acc;
mod mag;
mod mag_decl_table;
mod mag_ellipsoid_fitting;
pub mod params;
pub mod ppks;
mod util;

use core::sync::atomic::{AtomicU16, Ordering};

use chrono::NaiveDateTime;
use lin_alg::f32::{Quaternion, Vec3};
use num_enum::TryFromPrimitive;
use num_traits::Float;

#[cfg(feature = "defmt")]
use defmt::{Format, println};

pub use crate::{attitude::Ahrs, params::Params};

// C file with impl of EKF for quaternion rotation:
// https://github.com/pms67/EKF-Quaternion-Attitude-Estimation/blob/master/EKF.h
// https://github.com/pms67/EKF-Quaternion-Attitude-Estimation/blob/master/updateEKFQuatAtt.m

// Rotors: Quaternion generalization?
// https://marctenbosch.com/quaternions/

// Re quaternions, euler angles, and error state;
// Euler angles are mostly used to synthesize controllers such as PIDs.
// When it comes to state/attitude estimation, quaternions are more common.
// However you can't use a quaternion directly as a state in your EKF (technically you can with
// some care) since the the unit constraint will be violated. You can instead use an
// error-state EKF
// http://www.iri.upc.edu/people/jsola/JoanSola/objectes/notes/kinematics.pdf

pub const G: f32 = 9.80665; // Gravity, in m/s^2

pub const DEG_SCALE_1E8: f32 = 100_000_000.;

pub const UP: Vec3 = Vec3 {
    x: 0.,
    y: 0.,
    z: 1.,
};
pub const FORWARD: Vec3 = Vec3 {
    x: 0.,
    y: 1.,
    z: 0.,
};
pub const RIGHT: Vec3 = Vec3 {
    x: 1.,
    y: 0.,
    z: 0.,
};

#[derive(Clone, Copy)]
pub enum DeviceOrientation {
    YFwdXRight,
    YLeftXFwd,
}

impl Default for DeviceOrientation {
    fn default() -> Self {
        Self::YFwdXRight
    }
}

#[cfg_attr(feature = "defmt", derive(Format))]
#[cfg_attr(not(feature = "defmt"), derive(Debug))]
#[derive(Clone, Copy, Eq, PartialEq, TryFromPrimitive)]
#[repr(u8)]
pub enum FixType {
    // These values are from the UBLOX protocol.
    NoFix = 0,
    DeadReckoning = 1,
    Fix2d = 2,
    Fix3d = 3,
    /// GNSS + dead reckoning combined
    Combined = 4,
    TimeOnly = 5,
}

impl Default for FixType {
    fn default() -> Self {
        Self::NoFix
    }
}

#[derive(Clone, Default)]
// todo: Format unavail on NaiveDateTime
// #[cfg_attr(feature = "defmt", derive(Format))]
#[cfg_attr(not(feature = "defmt"), derive(Debug))]
/// In a format conducive to being parsed from the UBX PVT. (`UBX-NAV-PVT`)
/// Note: For position and elevation, we use the same units as Ublox reports; we
/// can convert to floats as required downstream.
pub struct Fix {
    /// This timestamp is local, eg synced from CAN bus.
    pub timestamp_s: f32,
    pub datetime: NaiveDateTime,
    pub type_: FixType,
    /// Degrees x 1e7
    pub lat_e7: i32,
    /// Degrees x 1e7
    pub lon_e7: i32,
    /// mm
    pub elevation_hae: i32,
    /// mm
    pub elevation_msl: i32,
    /// mm/s
    pub ground_speed: i32,
    /// North, east, down velocity; in that order. In mm/s.
    pub ned_velocity: [i32; 3],
    /// Degrees
    pub heading: Option<f32>, // only when valid.
    pub sats_used: u8,
    /// Position dilution of precision. Divided by 100.
    pub pdop: u16,
}

/// Represents sensor readings from a 6-axis accelerometer + gyro.
/// Accelerometer readings are in m/2^2. Gyroscope readings are in radians/s.
#[derive(Default, Clone)]
#[cfg_attr(feature = "defmt", derive(Format))]
#[cfg_attr(not(feature = "defmt"), derive(Debug))]
pub struct ImuReadings {
    /// Positive X: Accel towards right wing
    pub a_x: f32,
    /// Positive Y: Accel forwards
    pub a_y: f32,
    /// Positive X: Accel up
    pub a_z: f32,
    /// Positive pitch: nose up
    pub v_pitch: f32,
    /// Positive roll: left wing up
    pub v_roll: f32,
    /// Positive yaw: CW rotation.
    pub v_yaw: f32,
}

impl ImuReadings {
    /// We use this to assemble readings from the DMA buffer.
    pub fn from_buffer(buf: &[u8], accel_fullscale: f32, gyro_fullscale: f32) -> Self {
        // todo: Note: this mapping may be different for diff IMUs, eg if they use a different reading register ordering.
        // todo: Currently hard-set for ICM426xx.

        // Ignore byte 0; it's for the first reg passed during the `write` transfer.
        Self {
            a_x: interpret_accel_or_gyro(i16::from_be_bytes([buf[1], buf[2]]), accel_fullscale),
            a_y: interpret_accel_or_gyro(i16::from_be_bytes([buf[3], buf[4]]), accel_fullscale),
            a_z: interpret_accel_or_gyro(i16::from_be_bytes([buf[5], buf[6]]), accel_fullscale),
            v_pitch: interpret_accel_or_gyro(i16::from_be_bytes([buf[7], buf[8]]), gyro_fullscale),
            v_roll: interpret_accel_or_gyro(i16::from_be_bytes([buf[9], buf[10]]), gyro_fullscale),
            v_yaw: interpret_accel_or_gyro(i16::from_be_bytes([buf[11], buf[12]]), gyro_fullscale),
        }
    }
}

/// Convert raw accelerometer or gyro readings to their value in SI units, floating point.
/// Output: m/s^2 for accelerometer; rad/s for gyro.
pub fn interpret_accel_or_gyro(val: i16, fullscale: f32) -> f32 {
    (val as f32 / i16::MAX as f32) * fullscale
}

/// Utility function to print a quaternion by axes.
pub fn print_quat(quat: Quaternion, name: &str) {
    let (x_component, y_component, z_component) = quat.to_axes();

    #[cfg(feature = "defmt")]
    println!(
        "{} -- x{} y{} z{}",
        name, x_component, y_component, z_component
    );
}

/// `amount=1` means val1; amount=0 means val0
/// todo: Make this generic over things that mutiply? Missed opportunity to use this
/// todo in a few places for Vec and Matrix.
pub(crate) fn blend(val0: f32, val1: f32, amount: f32) -> f32 {
    let amount_inv = 1. - amount;
    val0 * amount_inv + val1 * amount
}

pub enum CalResult {
    Incomplete,
    Success((f32, f32, f32)),
    Fail,
    Disabled,
}

static ACC_CAL_VALS_LOGGED: AtomicU16 = AtomicU16::new(0);
static mut ACC_CAL_TOTAL: Vec3 = Vec3::new_zero();

/// Utility function to be called from firmware; designed to be run from a main loop or similar function.
pub fn cal_accel(
    calibrating: &mut bool,
    num_vals: u16,
    vals_to_skip: u16,
    imu_data: &ImuReadings,
) -> CalResult {
    if !*calibrating {
        return CalResult::Disabled;
    }

    // Actually all loops since start, including the skipped ones.
    let vals = ACC_CAL_VALS_LOGGED.fetch_add(1, Ordering::Relaxed);

    if vals > (num_vals + vals_to_skip) {
        *calibrating = false;
        ACC_CAL_VALS_LOGGED.store(0, Ordering::Release);

        let x = unsafe { ACC_CAL_TOTAL.x } / num_vals as f32;
        let y = unsafe { ACC_CAL_TOTAL.y } / num_vals as f32;
        let z = unsafe { ACC_CAL_TOTAL.z } / num_vals as f32 - G;

        const THRESH_A: f32 = 0.45; // m/s^2
        if x.abs() < THRESH_A && y.abs() < THRESH_A && z.abs() < THRESH_A {
            return CalResult::Success((x, y, z));
        } else {
            #[cfg(feature = "defmt")]
            println!(
                "Acc cal failed due to out of bounds value. X: {} Y: {} Z: {} Thresh: {}",
                x, y, z, THRESH_A
            );
            return CalResult::Fail;
        }
    } else if vals > vals_to_skip {
        let acc_data = Vec3::new(imu_data.a_x, imu_data.a_y, imu_data.a_z);

        // If any value is above the thresh, ommit it.
        // todo: You still need to fail the calibration if there are too many failures etc.
        const THRESH: f32 = 0.6; // m/s^2
        if acc_data.x.abs() < THRESH && acc_data.y.abs() < THRESH && (acc_data.z - G).abs() < THRESH
        {
            unsafe {
                ACC_CAL_TOTAL += acc_data;
            }
        } else {
            // Don't count this towards the number of logged values.
            ACC_CAL_VALS_LOGGED.store(vals - 1, Ordering::Release);
        }
    }
    CalResult::Incomplete
}
