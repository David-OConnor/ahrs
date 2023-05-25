//! Calculate attitude gven a 9-axis IMU.
//!
//!
//! todo: Use magnetic inclination, declination, and strength

// todo: Calibrate mag and acc based on known gravity value, and known magnetic strength at the current
// todo position, if known.

// todo: We are currently moving away from the AHRS Fusion port.

use core::f32::consts::TAU;

use num_traits::float::Float; // abs etc

use lin_alg2::f32::{Quaternion, Vec3};

use crate::{ppks::PositEarthUnits, FORWARD, RIGHT, UP};

use defmt::println;

pub struct Ahrs {
    pub attitude: Quaternion,
    att_from_gyros: Quaternion,
    att_from_acc: Quaternion,
    // att_from_mag: Quaternion,
    linear_acc_estimate: Vec3,
    linear_acc_confidence: f32, // todo?
    /// Time between updates, in seconds.
    dt: f32,
    /// Timestamp, in seconds.
    timestamp: f32,
}



impl Ahrs {
    pub fn new(dt: f32) -> Self {
        Self {
            attitude: Quaternion::new_identity(),
            att_from_gyros: Quaternion::new_identity(),
            att_from_acc: Quaternion::new_identity(),
            linear_acc_estimate: Vec3::new_zero(),
            linear_acc_confidence: 0.,
            dt,
            timestamp: 0.,
        }
    }

    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        let heading_from_prev = 0.; // todo

        let att_acc = att_from_accel(accel_data, heading_from_prev);

        let att_acc_w_lin_removed = att_from_accel(accel_data - self.linear_acc_estimate, heading_from_prev);

        let att_gyro = att_from_gyro(gyro_data, self.attitude, self.dt);

        let diff_acc_gyro = att_acc * att_gyro.inverse();
        let angle_diff_acc_gyro = diff_acc_gyro.angle();

        // todo: Is this up, or down?
        // todo: Inv, or normal?
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP);
        let lin_acc_estimate = att_acc - grav_axis_from_att_gyro * (accel_data.magnitude() - G); // todo QC

        self.linear_acc_estimate = lin_acc_estimate; // todo: DOn't take all of it; fuse with current value.


        // If the magntidue of the acceleration is above this value, we are under linear acceleration,
        // and should ignore the accelerometer.
        let acc_magnitude_thresh_upper = crate::G * 1.2;  // todo setting somewhere
        let acc_magnitude_thresh_lower = crate::G * 0.8;  // todo setting somewhere
        // let acc_magntude_thresh = 7.;  // todo setting somewhere
        // println!("mag: {}", accel_data.magnitude());

        let accel_magnitude = accel_data.magnitude();
        if accel_magnitude > acc_magnitude_thresh_upper || accel_magnitude < acc_magnitude_thresh_lower {
            // We are under linear acceleration; ignore the accelerometer, and use the integrated
            // gyro measurement.
            // (todo: estimate the linear component, remove that, then incorporate acc data)
        } else {
            // If not under much acceleration, re-cage our attitude.
            // todo: Partial, not full.
            self.attitude = att_acc;
        }

        // todo: Ways to identify linear acceleration:
        // - Greater or less than 1G of acceleration, if the accel is calibrated.
        // - Discontinuities or other anomolies when integrating accel-based attitude over time,
        // - or, along those lines, discontinuities etc when fusing with gyro.

        // Identify the angle difference in the vector between the current attitude estimate, and that
        // from the accelerometer alone.

        // self.att_from_gyros = att_gyro;

        // println!("att gyro: {:?}", att_gyro.x);


        match mag_data {
            Some(mag) => {
                let incliantion = -1.09;
                let att_mag = att_from_mag(mag, incliantion);
            }
            None => ()
        }

        self.attitude = att_gyro;
        self.timestamp += self.dt;

    }
}


/// Estimate attitude from accelerometer. This will fail when under
/// linear acceleration. Apply calibration prior to this step.
/// `heading` is in radians
pub fn att_from_accel(accel: Vec3, heading: f32) -> Quaternion {
    let accel_norm = accel.to_normalized();
    let att_without_heading = Quaternion::from_unit_vecs(UP, accel_norm);

    // Remove the final degree of freedom using heading.
    // todo: QC You're rotating around the up or down vec, and not a relative one.
    // let yaw_rotation = Quaternion::from_axis_angle(UP, heading); // todo: down or up?
    let yaw_rotation = Quaternion::from_axis_angle(accel_norm, heading); // todo: down or up?
    yaw_rotation * att_without_heading
}

/// Estimate attitude from magnetometer. This will fail when experiencing magnetic
/// interference, and is noisy in general. Apply calibration prior to this step.
/// Inclination is in radians.
// pub fn att_from_mag(mag: Vec3, posit: &PositEarthUnits) -> Quaternion {
pub fn att_from_mag(mag: Vec3, inclination: f32) -> Quaternion {
    let incl_rot = Quaternion::from_axis_angle(RIGHT, inclination);

    let mag_field_vec = incl_rot.rotate_vec(FORWARD);

    Quaternion::from_unit_vecs(mag_field_vec, mag.to_normalized())
}

/// Estimate attitude from gyroscopes. This will accumulate errors over time.
/// dt is in seconds.
pub fn att_from_gyro(gyro: Vec3, att_prev: Quaternion, dt: f32) -> Quaternion {
    let rot_x = Quaternion::from_axis_angle(RIGHT, gyro.x * dt);
    let rot_y = Quaternion::from_axis_angle(FORWARD, gyro.y * dt);
    let rot_z = Quaternion::from_axis_angle(UP, gyro.z * dt);

    // todo: Rotation order?
    rot_x * rot_y * rot_z * att_prev
}

pub fn get_linear_accel(accel: Vec3, att: Quaternion) -> Vec3 {
    let mut accel_vec_earth_ref = att.rotate_vec(accel.to_normalized());

    let accel_mag = accel.magnitude();

    accel_vec_earth_ref *= accel_mag;
    accel_vec_earth_ref.z -= crate::G;

    accel_vec_earth_ref
}
