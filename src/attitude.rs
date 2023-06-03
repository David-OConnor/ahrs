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
    pub linear_acc_estimate: Vec3,
    // linear_acc_confidence: f32, // todo?
    /// linear acceleration, per axis, integrated over time. We use this to
    /// remove biases, under the assumption that this should average out to 0, along
    /// each axis.
    linear_acc_cum: Vec3,
    /// Timestamp, in seconds.
    timestamp: f32,
    // todo: Config vars below. Put in a separate config struct to keep things explicit.
    // todo: And/or a separate state struct that contains things that change during operation.
    /// How far to look back when determining linear acceleration bias from cumulative values.
    lin_bias_timeout: f32,
    /// Time between updates, in seconds.
    dt: f32,
}

impl Ahrs {
    pub fn new(dt: f32) -> Self {
        Self {
            attitude: Quaternion::new_identity(),
            att_from_gyros: Quaternion::new_identity(),
            att_from_acc: Quaternion::new_identity(),
            linear_acc_estimate: Vec3::new_zero(),
            // linear_acc_confidence: 0.,
            linear_acc_cum: Default::default(),
            timestamp: 0.,
            lin_bias_timeout: 10.,
            dt,
        }
    }

    /// Update our AHRS solution given new gyroscope, accelerometer, and mag data.
    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        // todo: A system where any 2/3 of your sources can, in agreement, update the biases
        // todo of the third.

        // todo: Find covariance between various sensors to cancel biases, like positional
        // todo offset of the IMU from center of rotation?
        let accel_norm = accel_data.to_normalized();

        let att_acc = att_from_accel(accel_norm);
        let mut att_gyro = att_from_gyro(gyro_data, self.att_from_gyros, self.dt);

        // Heading from the previous attitude.
        let heading_prev = heading_from_att(self.attitude);

        let mag_heading_weight = 0.1; // todo: FUnction of dt
        let gyro_heading_weight = 1. - mag_heading_weight;

        let mut heading_fused = heading_prev;

        // todo: Separate function to fuse mag.
        match mag_data {
            Some(mag) => {
                let mag_norm = mag.to_normalized();
                // let incliantion = -1.09;
                // let att_mag = att_from_mag(mag_norm, incliantion);
                let heading_mag = heading_from_mag(mag);

                // Fuse heading from gyro with heading from mag.
                // heading_fused = (heading_prev * gyro_heading_weight + heading_mag * mag_heading_weight) / 2.;

                // todo: For now, we are using mag only for the heading.

                if unsafe { i } % 1000 == 0 {
                    println!("mag vec: x{} y{} z{}", mag_norm.x, mag_norm.y, mag_norm.z);

                    println!("Heading mag: {}", heading_mag);
                }
            }
            None => (),
        }

        // To perform this update, we assume the previous attitude has a valid heading.
        // let att_acc_w_heading = find_z_rot(self.att_from_gyros, accel_norm) * att_acc;
        let z_rotation = find_z_rot(heading_fused, accel_norm);
        let att_acc_w_heading = z_rotation * att_acc;

        let (update_gyro_from_acc, lin_acc_estimate) =
            self.handle_linear_acc(accel_data, att_acc_w_heading, att_gyro);

        let att_acc_w_lin_removed = att_from_accel((accel_data - lin_acc_estimate).to_normalized());

        // How much, as a portion of 1., to update the gyro attitude from the accel.
        // 1.0 means replace it.
        // let update_port_acc: f32 = 0.7 * self.dt;
        let update_port_acc: f32 = 0.05;
        let update_port_gyro: f32 = 1. - update_port_acc;

        let mut att_fused = att_gyro;

        if update_gyro_from_acc {
            // https://stackoverflow.com/questions/12374087/average-of-multiple-quaternions
            // There are some sophisticated methods of average quaternions; eg involving
            // eigen values. For now, we cheat by simply averaging their values; this is probably
            // ok if the quaternions are similar. It is more computationally efficient than the
            // proper method, but errors are higher if the quaternions are more different from each other.

            // todo: With lin removed.
            att_fused = Quaternion {
                w: (att_acc_w_heading.w * update_port_acc + att_gyro.w * update_port_gyro),
                x: (att_acc_w_heading.x * update_port_acc + att_gyro.x * update_port_gyro),
                y: (att_acc_w_heading.y * update_port_acc + att_gyro.y * update_port_gyro),
                z: (att_acc_w_heading.z * update_port_acc + att_gyro.z * update_port_gyro),
            }
            .to_normalized();

            att_fused = att_acc_w_heading; // todo: T!
        }

        self.attitude = att_fused;

        self.att_from_gyros = att_gyro;

        self.timestamp += self.dt;

        static mut i: u32 = 0;
        unsafe { i += 1 };
        if unsafe { i } % 1000 == 0 {
            // println!("Alignment: {}", acc_gyro_alignment);

            let euler = self.attitude.to_euler();

            let axis = self.attitude.axis();
            let angle = self.attitude.angle();

            let sign_x = -axis.x.signum();
            let sign_y = -axis.y.signum();
            let sign_z = -axis.z.signum();

            let x_component = (axis.project_to_vec(RIGHT) * angle).magnitude() * sign_x;
            let y_component = (axis.project_to_vec(FORWARD) * angle).magnitude() * sign_y;
            let z_component = (axis.project_to_vec(UP) * angle).magnitude() * sign_z;

            println!(
                "\n\nAxis rots: x{} y{} z{}",
                x_component, y_component, z_component
            );

            println!("Euler: p{} r{} y{}", euler.pitch, euler.roll, euler.yaw);

            println!("Acclen: {}", accel_data.magnitude());

            println!("\nHeading fused: {:?}\n", heading_fused);
        }
    }

    /// Attempt to separate linear from gravitational acceleration.
    /// Returns the estimate mof linear acceleration.
    fn handle_linear_acc(
        &mut self,
        accel_data: Vec3,
        att_acc: Quaternion,
        att_gyro: Quaternion,
    ) -> (bool, Vec3) {
        // todo: Ways to identify linear acceleration:
        // - Greater or less than 1G of acceleration, if the accel is calibrated.
        // - Discontinuities or other anomolies when integrating accel-based attitude over time,
        // - or, along those lines, discontinuities etc when fusing with gyro.

        // Identify the angle difference in the vector between the current attitude estimate, and that
        // from the accelerometer alone.

        // This is the up vector as assessed from the attitude from the gyro. It is equivalent to
        // the accelerometer's normalized vector when linear acceleration is 0.
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP); //Verified correct

        // An alignment of 1 indicates the estimated up direction between the accelerometer
        // and gyro match. A value of 0 means they're perpendicular.
        let acc_gyro_alignment = accel_data.to_normalized().dot(grav_axis_from_att_gyro);

        // Compute the difference, as a quaternion, between the attitude calulated from the accelerometer,
        // and the attitude calculated from the gyroscope. A notable difference implies linear acceleration.
        // todo: We are currently not using these
        // let diff_acc_gyro = att_acc_w_heading * att_gyro.inverse();
        // let angle_diff_acc_gyro = diff_acc_gyro.angle();

        // Estimate linear acceleration by comparing the accelerometer's normalized vector (indicating the
        // direction it resists gravity) with that estimated from the gyro. This is a proxy for linear
        // acceleration.

        // Some properties this should have:
        // -- If the grav vecs are aligned, we should see a linear accel along it: the accel's value - G * Z.
        // -- If the acc vec is 0, we should see a lin acc at 9.8G along the gyro's G vec
        // -- if the vecs are misaligned, we should see a lin acc IVO the acc's vector, taking grav into account.
        // -- If the vecs are 45 degrees out of alignment of equal mag, we should see 1G of lateral acceleration.
        // if they are more than 45 degrees out of alignment, there is >1G of lateral acc.

        // acc = lin + grav
        // lin = acc - (grav_axis_gyro * G)
        // For the purpose of this calculation, we are assuming the real gravitation axis is
        // that determined by the gyro.
        let lin_acc_estimate = accel_data - (grav_axis_from_att_gyro * G);

        // Store our linear acc estimate and accumulator before compensating for bias.
        self.linear_acc_estimate = lin_acc_estimate; // todo: DOn't take all of it; fuse with current value.
                                                     // todo: Be careful about floating point errors over time.
                                                     // todo: Toss extreme values?
                                                     // todo: Lowpass?
        self.linear_acc_cum += lin_acc_estimate * self.dt;

        // Important: This bias assumes acceleration evens out to 0 over time; this may or may not
        // be a valid assumption, under various conditions.
        // todo: Implement your timestamp, to keep bias computation over a set interval of time?
        let lin_acc_bias = if self.timestamp < 0.00001 {
            Vec3::new_zero()
        } else {
            // let _denominator = if self.timestamp > self.lin_bias_timeout {
            //     self.lin_bias_timeout
            // } else {
            //     self.timestamp
            // };
            self.linear_acc_cum / self.timestamp
        };

        let lin_acc_estimate_bias_removed = lin_acc_estimate - lin_acc_bias;

        // If the magntidue of the acceleration is above this value, we are under linear acceleration,
        // and should ignore the accelerometer.
        let acc_magnitude_thresh_upper = G * 1.2; // todo setting somewhere
        let acc_magnitude_thresh_lower = G * 0.8; // todo setting somewhere
                                                  // let accel_magnitude = accel_data.magnitude();

        let lin_acc_thresh = 0.4; // m/s^2
        let total_accel_thresh = 0.2; // m/s^2

        let mut update_gyro_from_acc = false;

        // If it appears there is negligible linear acceleration, update our gyro readings as appropriate.
        if (accel_data.magnitude() - G).abs() < total_accel_thresh {
            // We guess no linear acc since we're getting close to 1G. Note that
            // this will produce false positives in some cases.
            update_gyro_from_acc = true;
        } else if lin_acc_estimate_bias_removed.magnitude() < lin_acc_thresh {
            // If not under much acceleration, re-cage our attitude.
            // todo: Partial, not full.
            update_gyro_from_acc = true;
        }

        static mut i: u32 = 0;
        unsafe { i += 1 };
        if unsafe { i } % 1000 == 0 {
            // println!("Alignment: {}", acc_gyro_alignment);
            println!(
                "Lin bias: x{} y{} z{}",
                lin_acc_bias.x, lin_acc_bias.y, lin_acc_bias.z,
            );

            println!(
                "Lin x{} y{} z{}. mag{}",
                lin_acc_estimate_bias_removed.x,
                lin_acc_estimate_bias_removed.y,
                lin_acc_estimate_bias_removed.z,
                lin_acc_estimate_bias_removed.magnitude()
            );

            //     println!(
            //     "Diff acc gyro: {:?}, gyro grav x{} y{} z{}",
            //     angle_diff_acc_gyro,
            //     grav_axis_from_att_gyro.x,
            //     grav_axis_from_att_gyro.y,
            //     grav_axis_from_att_gyro.z
            // );
        }

        (update_gyro_from_acc, lin_acc_estimate_bias_removed)
    }
}

/// Estimate attitude from accelerometer. This will fail when under
/// linear acceleration. Apply calibration prior to this step.
/// Uses the previous attitude to rotate along the remaining degree of freedom (heading)
pub fn att_from_accel(accel_norm: Vec3) -> Quaternion {
    Quaternion::from_unit_vecs(UP, accel_norm)
}

/// Calculate heading, from an attitude;
fn heading_from_att(att: Quaternion) -> f32 {
    let axis = att.axis();
    let angle = att.angle();

    let sign_z = axis.z.signum();
    (axis.project_to_vec(UP) * angle).magnitude() * sign_z
}

/// Find the rotation around the Z axis associated with an attitude. We use this to apply
/// gyro and mag heading to the remaining degree of freedom on attitude from acceleration.
fn find_z_rot(heading: f32, accel_norm: Vec3) -> Quaternion {
    // Remove the final degree of freedom using heading. Rotate around UP in the earth frame.
    Quaternion::from_axis_angle(accel_norm, heading)
    // Quaternion::from_axis_angle(UP, heading)
}

/// Estimate attitude from magnetometer. This will fail when experiencing magnetic
/// interference, and is noisy in general. Apply calibration prior to this step.
/// Inclination is in radians.
// pub fn att_from_mag(mag: Vec3, posit: &PositEarthUnits) -> Quaternion {
pub fn att_from_mag(mag_norm: Vec3, inclination: f32) -> Quaternion {
    let incl_rot = Quaternion::from_axis_angle(RIGHT, inclination);

    let mag_field_vec = incl_rot.rotate_vec(FORWARD);

    Quaternion::from_unit_vecs(mag_field_vec, mag_norm)
}

/// Calculate heading, in radians, from the magnetometer's X and Y axes.
pub fn heading_from_mag(mag: Vec3) -> f32 {
    -(mag.y.atan2(mag.x))
}

/// Estimate attitude from gyroscopes. This will accumulate errors over time.
/// dt is in seconds.
pub fn att_from_gyro(gyro: Vec3, att_prev: Quaternion, dt: f32) -> Quaternion {
    // We use negative vectors due to the conventions; I don't have a grasp on it, but this appears
    // required for it to work.
    let rot_x = Quaternion::from_axis_angle(-RIGHT, gyro.x * dt);
    let rot_y = Quaternion::from_axis_angle(-FORWARD, gyro.y * dt);
    let rot_z = Quaternion::from_axis_angle(-UP, gyro.z * dt);

    // Rotation order?
    rot_z * rot_x * rot_y * att_prev
}

/// Estimate linear acceleration, given a known (estimated) attitude, and the acceleration vector.
pub fn get_linear_accel(accel: Vec3, att: Quaternion) -> Vec3 {
    let mut accel_vec_earth_ref = att.rotate_vec(accel.to_normalized());

    let accel_mag = accel.magnitude();

    accel_vec_earth_ref *= accel_mag;
    accel_vec_earth_ref.z -= G;

    accel_vec_earth_ref
}

/// Find the vector associated with linear acceleration induced by a position
/// offset in the accelerometer, leading rotations to induce accel.
/// We may subtract this from total accelerometer reading.
/// `posit_offset` is in meters. Rotation is in rad/s
pub fn find_accel_offset(posit_offset: Vec3, gyro_readings: Vec3) -> Vec3 {
    let rotation_x = Quaternion::from_axis_angle(RIGHT, gyro_readings.x);
    let rotation_y = Quaternion::from_axis_angle(FORWARD, gyro_readings.y);
    let rotation_z = Quaternion::from_axis_angle(UP, gyro_readings.z);

    // Order?
    let rotation = rotation_z * rotation_x * rotation_y;

    let rot_axis = rotation.axis();
    let rot_angle = rotation.angle();

    // todo: QC these
    posit_offset.cross(rot_axis * rot_angle)
}
