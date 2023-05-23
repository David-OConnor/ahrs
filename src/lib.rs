#![no_std]

//! This module contains code related to using sensor fusion to create an attitude platform. It also includes code for interpreting, integrating,
//! and taking the derivatives of sensor readings. Code here is device-agnostic.
//!
//! Some IMUs can integrate with a magnetometer and do some sensor fusion; we use
//! software since it's more general, flexible, and device-agnostic.
//!
//! From HyperShield: "Unscented will not offer any improvement over EKF. The reason for this is
//! that your quadrotor will mainly stay near hover (because of the assumption that gravity
//! will be the dominant acceleration) so the system will appear quasi-linear. If you want to go
//! agile flight then UKF might offer improvements, but then you can't rely on the gravity assumption.
//! There are attitude estimators like the hua filter that circumvents this (it estimates the
//! acceleration in the inertial frame and uses that for prediction instead assuming it's equal
//! to the gravity vector)."
//!
//! Python guide: https://github.com/rlabbe/Kalman-and-Bayesian-Filters-in-Python
//!
//! [Youtube video: Phil's Lab](https://www.youtube.com/watch?v=hQUkiC5o0JI)
//!
//! Important: We define X to be left/right (pitch), Y to be forward/back (roll, and
//! Z to be up/down (yaw)

mod ahrs_fusion;
mod attitude_platform;

pub use ahrs_fusion::Ahrs;

use lin_alg2::f32::{Mat3, Vec3, Quaternion};

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
    pub fn from_buffer(buf: &[u8], gyro_fullscale: f32, accel_fullscale: f32) -> Self {
        // todo: Note: this mapping may be different for diff IMUs, eg if they use a different reading register ordering.
        // todo: Currently hard-set for ICM426xx.

        // Ignore byte 0; it's for the first reg passed during the `write` transfer.
        Self {
            a_x: interpret_accel_or_gyro(i16::from_be_bytes([buf[1], buf[2]]), gyro_fullscale),
            a_y: interpret_accel_or_gyro(i16::from_be_bytes([buf[3], buf[4]]), gyro_fullscale),
            a_z: interpret_accel_or_gyro(i16::from_be_bytes([buf[5], buf[6]]), gyro_fullscale),
            v_pitch: interpret_accel_or_gyro(i16::from_be_bytes([buf[7], buf[8]]), accel_fullscale),
            v_roll: interpret_accel_or_gyro(i16::from_be_bytes([buf[9], buf[10]]), accel_fullscale),
            v_yaw: -interpret_accel_or_gyro(i16::from_be_bytes([buf[11], buf[12]]), accel_fullscale),
        }
    }
}


/// Update the attitude from the AHRS.
pub fn get_attitude(ahrs: &mut Ahrs, imu_readings: &ImuReadings, mag_readings: Option<Vec3>, dt: f32) -> Quaternion {
    // Gyro measurements - not really a vector.
    // In our IMU interpretation, we use direction references that make sense for our aircraft.
    // See `ImuReadings` field descriptions for this. Here, we undo it: The AHRS
    // fusion algorithm expects raw readings. (See the - signs there and here; they correspond)

    // ICM dirs: X right, Y back, Z up
    // Pitch nose up, roll right wing down, yaw CCW\
    // LIS3 mat: X left, Y back, z up

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
    let gyro_data = Vec3::new_zero();

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
            // todo: QC dir and order.
            let mag_data = Vec3 {
                x: m.y,
                y: m.x,
                z: m.z,
            };

            // ahrs.update(gyro_data, accel_data, mag_data, dt);
            // todo:TEmp. Put back version with mag data.
            ahrs.update_no_magnetometer(gyro_data, accel_data, dt);
        }
        None => {
            ahrs.update_no_magnetometer(gyro_data, accel_data, dt);
        }
    }

    ahrs.quaternion
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


