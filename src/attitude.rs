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

use crate::{ppks::PositEarthUnits, FORWARD, G, RIGHT, UP};

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
        // let fwd_shifted = self.attitude.rotate_vec(FORWARD);
        // todo: Store accumulated yaw? Find a solution here to not discard the yaw calculated
        // todo from gyros when you realign to the accel.
        let heading_from_prev = 0.; // todo: Temp

        let att_acc = att_from_accel(accel_data, heading_from_prev);

        let att_gyro = att_from_gyro(gyro_data, self.attitude, self.dt);

        let diff_acc_gyro = att_acc * att_gyro.inverse();
        let angle_diff_acc_gyro = diff_acc_gyro.angle();

        // This is the up vector as assessed from the attitude from the gyro. It is equivalent to
        // the accelerometer's normalized vector when linear acceleration is 0.
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP); //Verified correct

        // An alignment of 1 indicates the estimated up direction between the accelerometer
        // and gyro match. A value of 0 means they're perpendicular.
        let acc_gyro_alignment = accel_data.to_normalized().dot(grav_axis_from_att_gyro);

        // Estimate linear acceleration by comparing the accelerometer's normalized vector (indicating the
        // direction it resists gravity) with that estimated from the gyro. This is a proxy for linear
        // acceleration.

        // Some properties this should have:
        // -- If the grav vecs are aligned, we should see a linear accel along it: the accel's value - G * Z.
        // -- If the acc vec is 0, we should see a lin acc at 9.8G along the gyro's G vec
        // -- if the vecs are misaligned, we should see a lin acc IVO the acc's vector, taking grav into account.
        // -- If the vecs are 45 degrees out of alignment, we should see 1G of lateral acceleration.

        // todo: Def wrong: Consider the case of perpendicular alignment. We have a lot of lin acc; not 0!
        // todo: More examining and adjustment of this estimation A/R.
        let lin_acc_estimate = (accel_data - grav_axis_from_att_gyro * G) * acc_gyro_alignment;

        self.linear_acc_estimate = lin_acc_estimate; // todo: DOn't take all of it; fuse with current value.

        let att_acc_w_lin_removed =
            att_from_accel(accel_data - self.linear_acc_estimate, heading_from_prev);

        static mut i: u32 = 0;

        unsafe { i += 1 };

        if unsafe { i } % 1000 == 0 {
            println!("Alignment: {}", acc_gyro_alignment);

            println!(
                "Diff acc gyro: {:?}, gyro grav x{} y{} z{}",
                angle_diff_acc_gyro,
                grav_axis_from_att_gyro.x,
                grav_axis_from_att_gyro.y,
                grav_axis_from_att_gyro.z
            );

            println!(
                "Lin x{} y{} z{}. mag{}",
                self.linear_acc_estimate.x,
                self.linear_acc_estimate.y,
                self.linear_acc_estimate.z,
                lin_acc_estimate.magnitude()
            );
        }

        // If the magntidue of the acceleration is above this value, we are under linear acceleration,
        // and should ignore the accelerometer.
        let acc_magnitude_thresh_upper = G * 1.2; // todo setting somewhere
        let acc_magnitude_thresh_lower = G * 0.8; // todo setting somewhere
                                                  // let acc_magntude_thresh = 7.;  // todo setting somewhere
                                                  // println!("mag: {}", accel_data.magnitude());

        // let accel_magnitude = accel_data.magnitude();

        let lin_acc_thresh = 0.5; // m/s^2
        let total_accel_thresh = 0.2; // m/s^2

        if (accel_data.magnitude() - G).abs() < total_accel_thresh {
            // We guess no linear acc since we're getting close to 1G. Note that
            // this will produce false positives in some cases.
            self.attitude = att_acc;
        } else if lin_acc_estimate.magnitude() < lin_acc_thresh {
            // If not under much acceleration, re-cage our attitude.
            // todo: Partial, not full.
            self.attitude = att_acc;
        } else {
            self.attitude = att_gyro;
            // We are under linear acceleration; ignore the accelerometer, and use the integrated
            // gyro measurement.
            // (todo: estimate the linear component, remove that, then incorporate acc data)
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
            None => (),
        }

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
    // We use negative vectors due to the conventions; I don't have a grasp on it, but this appears
    // required for it to work.
    let rot_x = Quaternion::from_axis_angle(-RIGHT, gyro.x * dt);
    let rot_y = Quaternion::from_axis_angle(-FORWARD, gyro.y * dt);
    let rot_z = Quaternion::from_axis_angle(-UP, gyro.z * dt);

    // todo: Rotation order?
    rot_x * rot_y * rot_z * att_prev
}

pub fn get_linear_accel(accel: Vec3, att: Quaternion) -> Vec3 {
    let mut accel_vec_earth_ref = att.rotate_vec(accel.to_normalized());

    let accel_mag = accel.magnitude();

    accel_vec_earth_ref *= accel_mag;
    accel_vec_earth_ref.z -= G;

    accel_vec_earth_ref
}
