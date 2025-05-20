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
///
use lin_alg::f32::{Quaternion, Vec3};
use num_traits::Float;

use crate::{ppks, ppks::PositVelEarthUnits, Fix, UP};

/// Estimate linear acceleration by comparing the fused attitude's up direction (based primarily
/// on the gyro in the short term) to that from the accelerometer. This works well for short-duration
/// linear accelerations, but fails for long-term ones, such as an orbit.
pub fn from_gyro(accel_data: Vec3, att_gyro: Quaternion, acc_len_at_rest: f32) -> Vec3 {
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

    // let att_acc_non_norm = acc::att_from_accel(accel_data.to_normalized());
    // let att_diff_rot = att_acc_non_norm * att_gyro.inverse();

    // acc = lin + grav
    // lin = acc - (grav_axis_gyro * G)
    // For the purpose of this calculation, we are assuming the real gravitation axis is
    // that determined by the gyro.
    // This is in the aircraft's frame of reference.
    let lin_acc_estimate = accel_data - (grav_axis_from_att_gyro * acc_len_at_rest);

    static mut I: u32 = 0;
    unsafe { I += 1 };

    if unsafe { I } % 1000 == 0 {
        // if false {
        // println!("Ag: {}", _acc_gyro_alignment);
        // println!(
        //     "\n\nLin bias: x{} y{} z{}",
        //     self.cal.linear_acc_bias.x, self.cal.linear_acc_bias.y, self.cal.linear_acc_bias.z,
        // );

        // print_quat(att_diff_rot, "Att diff rot");
        // println!("Att diff rot angle: {}", att_diff_rot.angle());
        let a = grav_axis_from_att_gyro * acc_len_at_rest;
        // println!("Grav axis gyro x{} y{} z{}", a.x, a.y, a.z);

        // println!(
        //     "Diff acc gyro: {:?}, gyro grav x{} y{} z{}",
        //     angle_diff_acc_gyro,
        //     grav_axis_from_att_gyro.x,
        //     grav_axis_from_att_gyro.y,
        //     grav_axis_from_att_gyro.z
        // );
    }

    lin_acc_estimate
}

/// Estimate linear acceleration by differentiating GNSS-determined velocity.
/// Note that there is potentially a feedback loop from needing attitude to do this, and this
/// feeding back into attitude. Hopefully the gyro can keep this from being a problem.
pub(crate) fn from_gnss(fix: &Fix, fix_prev: &Fix, attitude: Quaternion) -> Vec3 {
    let d_v_earth =
        ppks::ned_vel_to_xyz(fix.ned_velocity) - ppks::ned_vel_to_xyz(fix_prev.ned_velocity);

    // IIR seems not appropriate, but you need a moving avg etc.

    let d_v = attitude.inverse().rotate_vec(d_v_earth);

    let d_t = fix.timestamp_s - fix_prev.timestamp_s;

    const EPS: f32 = 0.000001;
    if d_t.abs() < EPS {
        return Vec3::new_zero();
    }

    d_v / d_t
}

/// Estimate linear acceleration by examining short-term ground track turn radius.
pub(crate) fn from_ground_track(posits_recent: &[PositVelEarthUnits; 5], interval_s: f32) -> Vec3 {
    // In meters
    let mut posits_m = [Vec3::new_zero(); 5];

    for (i, posit) in posits_recent.iter().enumerate() {
        posits_m[i] = Vec3 {
            // todo: Acceptable precision using f32?
            x: posit.lon_e8 as f32 / crate::DEG_SCALE_1E8,
            y: posit.lat_e8 as f32 / crate::DEG_SCALE_1E8,
            z: posit.elevation_msl,
        }
    }

    // Ideally here we least-squares fit a circle in 3D space, but I'm not sure how to do this.
    // https://www.mathworks.com/matlabcentral/answers/475212-circle-least-squares-fit-for-3d-data

    // https://stackoverflow.com/questions/15481242/python-optimize-leastsq-fitting-a-circle-to-3d-set-of-points

    //  Coordinates of the barycenter
    let mut total = Vec3::new_zero();
    for p in &posits_m {
        total += *p;
    }
    let mean = total / posits_m.len() as f32;

    // gradient descent minimisation method ###
    // pnts = [[x[k], y[k], z[k]] for k in range(len(x))]
    // meanP = Point(xm, ym, zm) # mean point
    // Ri = [Point(*meanP).distance(Point(*pnts[k])) for k in range(len(pnts))] # radii to the points
    // Rm = math.fsum(Ri) / len(Ri) # mean radius
    // dR = Rm + 10 # difference between mean radii
    // alpha = 0.1
    // c = meanP
    // cArr = []
    // while dR  > eps:
    //     cArr.append(c)
    //     Jx = math.fsum([2 * (x[k] - c[0]) * (Ri[k] - Rm) / Ri[k] for k in range(len(Ri))])
    //     Jy = math.fsum([2 * (y[k] - c[1]) * (Ri[k] - Rm) / Ri[k] for k in range(len(Ri))])
    //     Jz = math.fsum([2 * (z[k] - c[2]) * (Ri[k] - Rm) / Ri[k] for k in range(len(Ri))])
    //     gradJ = [Jx, Jy, Jz] # find gradient
    //     c = [c[k] + alpha * gradJ[k] for k in range(len(c)) if len(c) == len(gradJ)] # find new centre point
    //     Ri = [Point(*c).distance(Point(*pnts[k])) for k in range(len(pnts))] # calculate new radii
    //     RmOld = Rm
    //     Rm = math.fsum(Ri) / len(Ri) # calculate new mean radius
    //     dR = abs(Rm - RmOld) # new difference between mean radii
    //
    // return Point(*c), Rm

    let turn_radius = 1.; //meters
    let speed = 1.;

    // todo: This should points towards the center of the circle; in the plane of rotation.
    let acc_direction = Vec3::new_zero();

    acc_direction * speed * 2. * turn_radius
}
