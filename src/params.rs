//! This module contains code that stores and updates flight parameters: Attitude, angular
//! rates, altitude etc.
//!
//! Sync with Corvus.

// use crate::{ppks::PositionFused};

use lin_alg::f32::{Quaternion, Vec3};

use crate::{ppks::PositFused, Ahrs, DeviceOrientation, ImuReadings};

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
    /// Linear acceleration, in m/s^2, relative to the AHRS.
    pub accel_linear: Vec3,
    // todo: Vec3 for these?

    // Velocity
    pub v_x: f32,
    pub v_y: f32,
    pub v_z: f32,

    pub v_z_baro: f32,

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
    //
    // pub mag_data: Option<Vec3>,
}

impl Params {
    /// Update params with IMU readings, and attitude. If filtering IMU readings, to so before
    /// running this.
    pub fn update_from_imu_readings(
        &mut self,
        imu_readings: &ImuReadings,
        mag_readings: Option<Vec3>,
        ahrs: &mut Ahrs,
    ) {
        let (accel_data, gyro_data, mag_data) = match ahrs.config.orientation {
            DeviceOrientation::YFwdXRight => {
                let accel_data = Vec3 {
                    x: imu_readings.a_x,
                    y: imu_readings.a_y,
                    z: imu_readings.a_z,
                };

                let gyro_data = Vec3 {
                    x: imu_readings.v_pitch,
                    y: imu_readings.v_roll,
                    z: imu_readings.v_yaw,
                };

                // Invert x and y for mag due to the coordinate system it uses.
                let mag_data = mag_readings.map(|m| Vec3 {
                    x: -m.x, // negative due to our mag's coord system.
                    y: -m.y,
                    z: m.z,
                });

                (accel_data, gyro_data, mag_data)
            }
            DeviceOrientation::YLeftXFwd => {
                let accel_data = Vec3 {
                    x: -imu_readings.a_y,
                    y: imu_readings.a_x,
                    z: imu_readings.a_z,
                };

                let gyro_data = Vec3 {
                    x: -imu_readings.v_roll,
                    y: imu_readings.v_pitch,
                    z: imu_readings.v_yaw,
                };

                // Invert x and y for mag due to the coordinate system it uses.
                let mag_data = mag_readings.map(|m| Vec3 {
                    x: m.y,
                    y: -m.x,
                    z: m.z,
                });

                (accel_data, gyro_data, mag_data)
            }
        };

        ahrs.update(gyro_data, accel_data, mag_data);

        // Now that we've estimated attitude, and applied calibrations, update parameters.
        self.attitude = ahrs.attitude;
        self.accel_linear = ahrs.lin_acc_fused;

        let acc_calibrated = ahrs.acc_calibrated;
        let gyro_calibrated = ahrs.gyro_calibrated;

        self.a_x = acc_calibrated.x;
        self.a_y = acc_calibrated.y;
        self.a_z = acc_calibrated.z;

        // Calculate angular acceleration. Do this before updating velocities, since we use
        // the previous ones here.
        self.a_pitch = (gyro_calibrated.x - self.v_pitch) / ahrs.dt;
        self.a_roll = (gyro_calibrated.y - self.v_roll) / ahrs.dt;
        self.a_yaw = (gyro_calibrated.z - self.v_yaw) / ahrs.dt;

        self.v_pitch = gyro_calibrated.x;
        self.v_roll = gyro_calibrated.y;
        self.v_yaw = gyro_calibrated.z;
    }
}
