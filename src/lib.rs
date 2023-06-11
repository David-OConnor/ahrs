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

pub mod attitude;
pub mod params;
pub mod ppks;

pub use crate::{attitude::Ahrs, params::Params};

use chrono::NaiveDateTime;
use lin_alg2::f32::{Mat3, Quaternion, Vec3};
use num_enum::{TryFromPrimitive, TryFromPrimitiveError};

// todo: Try this : https://github.com/Mayitzin/ahrs/blob/master/ahrs/filters/ekf.py

// use defmt::println;

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

/// Apply this linear map to raw IMU readings to get calibrated ones, that should always
/// return 1G of acceleration if no linear acceleration is applied.
pub struct AccelCal {
    pub slope_x: f32,
    pub intercept_x: f32,
    pub slope_y: f32,
    pub intercept_y: f32,
    pub slope_z: f32,
    pub intercept_z: f32,
}

impl Default for AccelCal {
    fn default() -> Self {
        Self {
            slope_x: 1.,
            intercept_x: 0.,
            slope_y: 1.,
            intercept_y: 0.,
            slope_z: 1.,
            intercept_z: 0.,
        }
    }
}

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

#[derive(Default)]
/// In a format conducive to being parsed from the UBX PVT. (`UBX-NAV-PVT`)
/// Note: For position and elevation, we use the same units as Ublox reports; we
/// can convert to floats as required downstream.
pub struct Fix {
    /// This timestamp is local, eg synced from CAN bus.
    pub timestamp_s: f32,
    pub datetime: NaiveDateTime,
    pub type_: FixType,
    // /// Degrees
    // pub lat: f64,
    // /// Degrees
    // pub lon: f64,
    /// Degrees x 1e7
    pub lat: i32,
    /// Degrees x 1e7
    pub lon: i32,
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
#[derive(Default)]
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

/// Update the attitude from the AHRS.
pub fn get_attitude(
    ahrs: &mut Ahrs,
    imu_readings: &ImuReadings,
    mag_readings: Option<Vec3>,
) -> Quaternion {
    // Gyro measurements - not really a vector.
    // In our IMU interpretation, we use direction references that make sense for our aircraft.
    // See `ImuReadings` field descriptions for this. Here, we undo it: The AHRS
    // fusion algorithm expects raw readings. (See the - signs there and here; they correspond)

    // ICM dirs: X right, Y back, Z up
    // Pitch nose up, roll right wing down, yaw CCW
    // LIS3 mat: X left, Y back, z up

    // Note; The ICM42688 uses the right hand rule for relating accelerometer and gyro readings.

    let accel_data = Vec3 {
        x: imu_readings.a_x,
        y: imu_readings.a_y, // negative due to our IMU's coord system.
        z: imu_readings.a_z,
    };

    let gyro_data = Vec3 {
        x: imu_readings.v_pitch,
        y: imu_readings.v_roll,
        z: imu_readings.v_yaw,
    };
    // let gyro_data = Vec3::new_zero();

    // Apply calibration
    // todo: Come back to this.
    // gyro_data = madgwick::apply_cal_inertial(
    //     gyro_data,
    //     ahrs.calibration.gyro_misalignment.clone(),
    //     ahrs.calibration.gyro_sensitivity,
    //     ahrs.calibration.gyro_offset,
    // );
    // accel_data = madgwick::apply_cal_inertial(
    //     accel_data,
    //     ahrs.calibration.accel_misalignment.clone(),
    //     ahrs.calibration.accel_sensitivity,
    //     ahrs.calibration.accel_offset,
    // );

    // let magnetometer = madgwick::apply_cal_magnetic(magnetometer, softIronMatrix, hardIronOffset);

    // todo: Consider putting this offset bit back later
    // Update gyroscope offset correction algorithm
    // let gyro_data_with_offset = ahrs.offset.update(gyro_data);
    // let gyro_data_with_offset = gyro_data;

    match mag_readings {
        Some(m) => {
            let mag_data = Vec3 {
                x: -m.y, // negative due to our mag's coord system.
                y: -m.x,
                z: m.z,
            };

            ahrs.update(gyro_data, accel_data, Some(mag_data));
        }
        None => {
            ahrs.update(gyro_data, accel_data, None);
        }
    }

    ahrs.attitude
}

/// Output: m/s^2, or Output: rad/s.
pub fn interpret_accel_or_gyro(val: i16, fullscale: f32) -> f32 {
    (val as f32 / i16::MAX as f32) * fullscale
}

// This calibration functionality is from [AHRS](https://github.com/xioTechnologies/Fusion)

pub struct ImuCalibration {
    pub gyro_misalignment: Mat3,
    pub gyro_sensitivity: Vec3,
    pub gyro_offset: Vec3,
    pub accel_misalignment: Mat3,
    pub accel_sensitivity: Vec3,
    pub accel_offset: Vec3,
    pub soft_iron_matrix: Mat3,
    pub hard_iron_offset: Vec3,
}

impl Default for ImuCalibration {
    #[rustfmt::skip]
    fn default() -> Self {
        Self {
            gyro_misalignment: Mat3 {
                data: [
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0
                ],
            },
            gyro_sensitivity: Vec3::new(1.0, 1.0, 1.0),
            gyro_offset: Vec3::new(0.0, 0.0, 0.0),
            accel_misalignment: Mat3 {
                data: [
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0
                ],
            },
            accel_sensitivity: Vec3::new(1.0, 1.0, 1.0),
            accel_offset: Vec3::new(0.0, 0.0, 0.0),
            soft_iron_matrix: Mat3 {
                data: [
                    1.0, 0.0, 0.0,
                    0.0, 1.0, 0.0,
                    0.0, 0.0, 1.0
                ],
            },
            hard_iron_offset: Vec3::new(0.0, 0.0, 0.0),
        }
    }
}

/// Gyroscope and accelerometer calibration model. Returns calibrated measurement.
pub fn apply_cal_inertial(
    uncalibrated: Vec3,
    misalignment: Mat3,
    sensitivity: Vec3,
    offset: Vec3,
) -> Vec3 {
    misalignment * (uncalibrated - offset).hadamard_product(sensitivity)
}

/// Magnetometer calibration model. Returns calibrated measurement.
pub fn apply_cal_magnetic(
    uncalibrated: Vec3,
    soft_iron_matrix: Mat3,
    hard_iron_offset: Vec3,
) -> Vec3 {
    soft_iron_matrix * uncalibrated - hard_iron_offset
}

/// Calibrate the IMU, by taking a series of series while on a level surface.
pub fn calibrate() -> ImuCalibration {
    // todo: average? lowpass?
    Default::default()
}
