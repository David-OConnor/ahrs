//! This module contains code that stores and updates flight parameters: Attitude, angular
//! rates, altitude etc.
//!
//! Sync with Corvus.

// use crate::{ppks::PositionFused};

use crate::{ppks::PositFused, Ahrs, ImuReadings};

use lin_alg2::f32::{Quaternion, Vec3};

/// Aircraft flight parameters, at a given instant. Pitch and roll rates are in the aircraft's
/// frame of reference.
#[derive(Default, Clone)]
pub struct Params {
    pub posit_fused: PositFused,
    // /// Latitude in degrees; fused.
    // pub lat_e8: i64,
    // /// Longitude in degrees; fused.
    // pub lon_e8: i64,
    // /// Altitude fused from GNSS, IMU, and maybe baro. In meters.
    // pub alt_msl_fused: f32,
    // /// MSL altitude in meters QFE (takeoff location is 0), from a barometer.
    pub alt_msl_baro: f32,
    /// AGL altitude in meters, from the Time of flight sensor.
    pub alt_tof: Option<f32>,
    //
    // pub s_pitch: f32,
    // pub s_roll: f32,
    // /// Ie heading
    // pub s_yaw_heading: f32,
    /// Quaternion of the attitude.
    pub attitude: Quaternion,
    // todo: AHRS quaternion field, or leave that as part of the `AHRS` struct?
    pub accel_linear: Vec3,

    // Velocity
    pub v_x: f32,
    pub v_y: f32,
    pub v_z: f32,

    pub v_pitch: f32,
    pub v_roll: f32,
    pub v_yaw: f32,

    // Acceleration
    pub a_x: f32,
    pub a_y: f32,
    pub a_z: f32,

    pub a_pitch: f32,
    pub a_roll: f32,
    pub a_yaw: f32,

    pub mag_x: f32,
    pub mag_y: f32,
    pub mag_z: f32,
}

impl Params {
    /// Update params with IMU readings, and attitude. If filtering IMU readings, to so before
    /// running this.
    pub fn update_from_imu_readings(
        &mut self,
        imu_data: &ImuReadings,
        mag_data: Option<Vec3>,
        // attitude: Quaternion,
        // accel_linear: Vec3,
        // acc_cal: &crate::ImuCalibration,
        ahrs: &mut Ahrs,
        dt: f32,
    ) {
        // todo: This is a good place to apply IMU calibration.

        // todo: Is accel cal required given your bias algo?

        // Calculate angular acceleration. Do this before updating velocities, since we use
        // the prev ones here.
        self.a_pitch = (imu_data.v_pitch - self.v_pitch) / dt;
        self.a_roll = (imu_data.v_roll - self.v_roll) / dt;
        self.a_yaw = (imu_data.v_yaw - self.v_yaw) / dt;

        // Apply filtered gyro and accel readings directly to self.
        self.v_pitch = imu_data.v_pitch;
        self.v_roll = imu_data.v_roll;
        self.v_yaw = imu_data.v_yaw;

        let cal = &ahrs.cal;

        // Update in-place to since we pass this to the attitude-finder.
        let mut imu_data2 = (*imu_data).clone();

        imu_data2.a_x = imu_data.a_x * cal.acc_slope.x - cal.acc_bias.x;
        imu_data2.a_y = imu_data.a_y * cal.acc_slope.y - cal.acc_bias.y;
        imu_data2.a_z = imu_data.a_z * cal.acc_slope.z - cal.acc_bias.z;

        self.a_x = imu_data2.a_x;
        self.a_y = imu_data2.a_y;
        self.a_z = imu_data2.a_z;

        self.attitude = crate::get_attitude(ahrs, &imu_data2, mag_data);

        self.accel_linear = ahrs.linear_acc_estimate;

        if let Some(mag) = mag_data {
            self.mag_x = mag.x;
            self.mag_y = mag.y;
            self.mag_z = mag.z;
        }
    }
}
