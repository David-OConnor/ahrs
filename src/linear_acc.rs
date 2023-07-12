//! This module contains code related to isolating linear acceleration
//! from gravitational acceleration.
//!
//! Uses:
//! #1: Determining how much faith (and wight) to put into the accelerometer reading for
/// attitude determination.
/// #2: Providing an acc solution that's closer to the true one for this fusing.
/// #3: Removing linear acceleration when computing position from dead-reckoning
///
///  Ways to identify linear acceleration:
/// - Greater or less than 1G of acceleration, if the accel is calibrated.
/// - Discontinuities or other anomolies when integrating accel-based attitude over time,
/// - or, along those lines, discontinuities etc when fusing with gyro.
use crate::{attitude, ppks, Fix, G, UP};

use lin_alg2::f32::{Quaternion, Vec3};

use defmt::println;

/// Estimate linear acceleration by comparing the fused attitude's up direction (based primarily
/// on the gyro in the short term) to that from the accelerometer. This works well for short-duration
/// linear accelerations, but fails for long-term ones, such as an orbit.
pub fn from_gyro(
    accel_data: Vec3,
    // accel_mag: Vec3,
    // att_acc: Quaternion,
    // // att gyro must have heading removed, or acc must have heading added.
    att_gyro: Quaternion,
    acc_len_at_rest: f32,
    lin_acc_bias: Vec3,
) -> (Vec3, Vec3) {
    // This is the up vector as assessed from the attitude from the gyro. It is equivalent to
    // the accelerometer's normalized vector when linear acceleration is 0.
    let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP);

    // Estimate linear acceleration by comparing the accelerometer's normalized vector (indicating the
    // direction it resists gravity) with that estimated from the gyro. This is a proxy for linear
    // acceleration.

    // For this, we postulate that the gyro's attitude is correct, and therefore the force
    // for gravity is that axis's *up* vector, multiplied by G. The difference between the accelerometer
    // readings and this, therefore is linear acceleration.
    // This linear acc estimate is in earth coords.

    // todo: Project to elim heading effects?
    let att_acc_non_norm = attitude::att_from_accel(accel_data.to_normalized());
    let att_diff_rot = att_acc_non_norm * att_gyro.inverse();

    // acc = lin + grav
    // lin = acc - (grav_axis_gyro * G)
    // For the purpose of this calculation, we are assuming the real gravitation axis is
    // that determined by the gyro.
    // This is in the aircraft's frame of reference.
    let lin_acc_estimate = accel_data - (grav_axis_from_att_gyro * acc_len_at_rest);

    let lin_acc_estimate_bias_removed = lin_acc_estimate - lin_acc_bias;

    static mut I: u32 = 0;
    unsafe { I += 1 };

    if unsafe { I } % 1000 == 0 {
        // if false {
        // println!("Ag: {}", _acc_gyro_alignment);
        // println!(
        //     "\n\nLin bias: x{} y{} z{}",
        //     self.cal.linear_acc_bias.x, self.cal.linear_acc_bias.y, self.cal.linear_acc_bias.z,
        // );

        println!(
            "Lin acc: x{} y{} z{}. mag{}",
            lin_acc_estimate_bias_removed.x,
            lin_acc_estimate_bias_removed.y,
            lin_acc_estimate_bias_removed.z,
            lin_acc_estimate_bias_removed.magnitude()
        );

        // print_quat(att_diff_rot, "Att diff rot");
        // println!("Att diff rot angle: {}", att_diff_rot.angle());
        let a = grav_axis_from_att_gyro * acc_len_at_rest;
        // println!("Grav axis gyro x{} y{} z{}", a.x, a.y, a.z);
        println!("Acc x{} y{} z{}", accel_data.x, accel_data.y, accel_data.z);
        //
        // println!(
        //     "Lin: x{} y{} z{}. mag{}",
        //     lin_acc_estimate.x,
        //     lin_acc_estimate.y,
        //     lin_acc_estimate.z,
        //     lin_acc_estimate.magnitude() // lin_acc_estimate.x,
        //                                               // lin_acc_estimate.y,
        //                                               // lin_acc_estimate.z,
        //                                               // lin_acc_estimate.magnitude()
        // );

        // println!(
        //     "Diff acc gyro: {:?}, gyro grav x{} y{} z{}",
        //     angle_diff_acc_gyro,
        //     grav_axis_from_att_gyro.x,
        //     grav_axis_from_att_gyro.y,
        //     grav_axis_from_att_gyro.z
        // );
    }

    // if unsafe { i } % 100 == 0 {
    //     if !update_gyro_from_acc {
    //         println!("Under lin acc");
    //     }
    // }

    (lin_acc_estimate, lin_acc_estimate_bias_removed)
}

/// Estimate linear acceleration by differentiating GNSS-determined velocity.
pub fn from_gnss(fix: &Fix, fix_prev: &Fix) -> Vec3 {
    let d_v = ppks::ned_vel_to_xyz(fix.ned_velocity) - ppks::ned_vel_to_xyz(fix_prev.ned_velocity);
    let d_t = fix.timestamp_s - fix_prev.timestamp_s;

    d_v / d_t
}
