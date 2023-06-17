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

pub struct AhrsConfig {
    /// How far to look back when determining linear acceleration bias from cumulative values. In seconds.
    pub lin_bias_lookback: f32,
    /// Look back time, in seconds. of mag heading change rate vice gyro heading change rate
    pub mag_diff_lookback: f32,
    /// Difference, in radians/s between mag and gyro heading change. Used to assess if we should
    /// update the gyro heading from the mag.
    pub mag_gyro_diff_thresh: f32,
    /// This value affects how much pitch and roll from the accelerometer are used to
    /// update the gyro. A higher value means more of an update. If this value is 0, the gyro
    /// will not be updated from the acc. If it x dt is 1., the gyro will be synchronized
    /// with the acc. Reasonable values may be around 1. If this x dt is > 1, anomolous behavior
    /// will occur. Higher values lead to more bias from linear acceleration. Low values lead to more
    /// gyro drift. Values from 0.1 to 10 may be optimal.
    ///
    /// This value can be thought of as the 1 / the number of seconds to correct a gyro reading to match the
    /// accelerometer. Keep this in mind re expectations of gyro drift rate.
    pub update_from_acc_amt: f32,
    pub update_port_mag_heading: f32,
    // /// Assume there's minimal linear acceleration if accelerometer magnitude falls between
    // /// G x these values (low, high).
    // /// todo: Alternative: Adjust adjustment-towards-acc weight value based on this, vice
    // /// todo updating or not.
    // pub acc_mag_threshold_no_lin: (f32, f32),
    /// If total acclerometer reading is within this value of G, update gyro from acc.
    pub total_accel_thresh: f32, // m/s^2
    /// If estimated linear acceleration magnitude is greater than this, don't update gyro
    /// from acc.
    pub lin_acc_thresh: f32, // m/s^2
    pub calibration: crate::ImuCalibration,
    /// Time, in seconds, after initialization, to start alignment procedure. This is set up so
    /// powering the device on by plugging it in, etc, doesn't interfere.
    pub start_alignment_time: u8,
    /// Time, in seconds, of the alignment.
    pub alignment_duration: u8,
}

impl Default for AhrsConfig {
    fn default() -> Self {
        Self {
            lin_bias_lookback: 10.,
            mag_diff_lookback: 10.,
            mag_gyro_diff_thresh: 0.01,
            update_from_acc_amt: 0.4,
            update_port_mag_heading: 0.1,
            // acc_mag_threshold_no_lin: (0.8, 1.2),
            total_accel_thresh: 0.2, // m/s^2
            lin_acc_thresh: 0.4,     // m/s^2
            calibration: Default::default(),
            start_alignment_time: 2,
            alignment_duration: 2,
        }
    }
}

pub struct Ahrs {
    pub attitude: Quaternion,
    att_from_gyros: Quaternion,
    att_from_acc: Quaternion,
    /// We use these to stored headings to track magnetometer health over time.
    heading_mag: Option<f32>,
    heading_gyro: f32,
    // att_from_mag: Quaternion,
    pub linear_acc_estimate: Vec3,
    // linear_acc_confidence: f32, // todo?
    /// linear acceleration, per axis, integrated over time. We use this to
    /// remove biases, under the assumption that this should average out to 0, along
    /// each axis.
    linear_acc_cum: Vec3,
    /// Track recent changing in mag heading, per change in gyro heading, in degrees.
    /// We use this to assess magnetometer health, ie if it's being interfered with.
    /// This is `None` if there are no mag reading provided.
    recent_dh_mag__dh_gyro: Option<f32>,
    /// Timestamp, in seconds.
    timestamp: f32,
    pub config: AhrsConfig,
    /// Time between updates, in seconds.
    dt: f32,
    /// We set this var upon the first update; this forces the gyro to take a full update from the
    /// accelerometer. Without this, we maay experience strong disagreement between the gyro and acc
    /// at start, since the gyro initializes to level, regardless of actual aircraft attitude.
    initialized: bool,
}

impl Ahrs {
    pub fn new(dt: f32) -> Self {
        Self {
            attitude: Quaternion::new_identity(),
            att_from_gyros: Quaternion::new_identity(),
            att_from_acc: Quaternion::new_identity(),
            heading_mag: None,
            heading_gyro: 0.,
            linear_acc_estimate: Vec3::new_zero(),
            recent_dh_mag__dh_gyro: None,
            // linear_acc_confidence: 0.,
            linear_acc_cum: Default::default(),
            timestamp: 0.,
            config: AhrsConfig::default(),
            dt,
            initialized: false,
        }
    }

    /// Update our AHRS solution given new gyroscope, accelerometer, and mag data.
    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        // todo: A system where any 2/3 of your sources can, in agreement, update the biases
        // todo of the third.

        // todo: Find covariance between various sensors to cancel biases, like positional
        // todo offset of the IMU from center of rotation?
        let accel_norm = accel_data.to_normalized();

        // Estimate attitude from raw accelerometer and gyro data. Note that
        // The gyro data reguarly receives updates from the acc and mag.
        let att_acc = att_from_accel(accel_norm);

        let mut att_gyro = att_from_gyro(gyro_data, self.att_from_gyros, self.dt);

        let heading_gyro = heading_from_att(att_gyro);

        let z_rotation = find_z_rot(heading_gyro, accel_norm);
        let att_acc_w_heading = z_rotation * att_acc;

        // See comment on the `initialized` field.
        if !self.initialized {
            self.att_from_gyros = att_acc_w_heading;
            att_gyro = self.att_from_gyros;
            self.initialized = true;
        }

        let mut heading_fused = heading_gyro;

        // Fuse with mag data if available.
        match mag_data {
            Some(mag) => {
                let mag_norm = mag.to_normalized();
                // let incliantion = -1.09;
                // let att_mag = att_from_mag(mag_norm, incliantion);
                let heading_mag = heading_from_mag(mag);

                // Assess magnetometer health by its comparison in rate change compared to the gyro.
                match self.heading_mag {
                    Some(heading_mag_prev) => {
                        // todo: Find avg over many readings.
                        // todo: Since even a messed up mag seems to show constant readings
                        // todo when there is no rotation, consider only logging values here if
                        // todo dh/dt exceeds a certain value.
                        self.recent_dh_mag__dh_gyro = Some(
                            (heading_mag - heading_mag_prev) / self.dt
                                - (heading_gyro - self.heading_gyro) / self.dt,
                        );

                        // Fuse heading from gyro with heading from mag.
                        if self.recent_dh_mag__dh_gyro.unwrap().abs()
                            < self.config.mag_gyro_diff_thresh
                        {
                            // heading_fused = (heading_prev * gyro_heading_weight + heading_mag * mag_heading_weight) / 2.;
                        }

                        let update_port_gyro_heading = 1. - self.config.update_port_mag_heading;
                    }
                    None => {
                        self.recent_dh_mag__dh_gyro = None;
                    }
                }

                // todo: For now, we are using mag only for the heading.

                self.heading_mag = Some(heading_mag);

                // if unsafe { i } % 1000 == 0 {
                if false {
                    println!("mag vec: x{} y{} z{}", mag_norm.x, mag_norm.y, mag_norm.z);

                    println!("Heading mag: {}", heading_mag);
                    println!("Heading gyro: {}", heading_gyro);
                }
            }
            None => {
                self.recent_dh_mag__dh_gyro = None;
            }
        }

        // Consider this flow: mag updates gyro heading. Att updates gyro pitch and roll.
        // Gyro turns into fused?

        let (update_gyro_from_acc, lin_acc_estimate) =
            self.handle_linear_acc(accel_data, att_acc_w_heading, att_gyro);

        let att_acc_w_lin_removed = att_from_accel((accel_data - lin_acc_estimate).to_normalized());
        let att_acc_w_lin_removed_and_heading = z_rotation * att_acc_w_lin_removed;

        // let att_acc_w_lin_removed_and_heading = att_acc_w_heading;

        let mut att_fused = att_gyro;

        // todo: Instead of a binary update-or-not, consider weighing the slerp value based
        // todo on how much lin acc we assess, or how much uncertainly in lin acc.
        if update_gyro_from_acc {
            // For now, always update gyro readings from the acc. This may be OK given we've removed ac.
            // Experiment and see, once you have this hooked into corvus connected to preflight for a visual
            // rep.
            // if true {
            // https://stackoverflow.com/questions/12374087/average-of-multiple-quaternions
            // There are some sophisticated methods of average quaternions; eg involving
            // eigen values. For now, we cheat by simply averaging their values; this is probably
            // ok if the quaternions are similar. It is more computationally efficient than the
            // proper method, but errors are higher if the quaternions are more different from each other.

            // let update_port_acc: f32 = self.config.update_port_acc * self.dt;
            // let update_port_gyro: f32 = 1. - self.config.update_port_acc;

            // println!("Updating from acc");
            // todo: With lin removed.

            // Apply a rotation of the gyro solution towards the acc solution, if we think we are not under
            // much linear acceleration.
            // let rot_gyro_to_acc = att_acc_w_heading * att_gyro.inverse();
            let rot_gyro_to_acc = att_acc_w_lin_removed_and_heading * att_gyro.inverse();

            let rot_to_apply_to_gyro = Quaternion::new_identity()
                .slerp(rot_gyro_to_acc, self.config.update_from_acc_amt * self.dt);

            // if unsafe { i } % 1000 == 0 {
            if false {
                let (x_component, y_component, z_component) = att_to_axes(rot_gyro_to_acc);

                // println!(
                //     "\n\nRot gyro to get to acc: x{} y{} z{}",
                //     x_component, y_component, z_component
                // );
                //
                // let (x_component, y_component, z_component) = att_to_axes(rot_to_apply_to_gyro);
                //
                println!(
                    "Rot to apply: x{} y{} z{}",
                    x_component, y_component, z_component
                );
            }

            // todo: Now, you need to interpolate between the identity quat and this one,
            // todo since you're not taking teh whole update

            att_fused = rot_to_apply_to_gyro * att_gyro;

            // let update_port_acc = self.config.update_port_acc;
            // att_fused = Quaternion {
            //     w: (att_acc_w_heading.w * update_port_acc + att_gyro.w * update_port_gyro),
            //     x: (att_acc_w_heading.x * update_port_acc + att_gyro.x * update_port_gyro),
            //     y: (att_acc_w_heading.y * update_port_acc + att_gyro.y * update_port_gyro),
            //     z: (att_acc_w_heading.z * update_port_acc + att_gyro.z * update_port_gyro),
            // }
            //     .to_normalized();

            // Careful; this replaces the gyro data.
            // self.att_from_gyros = att_fused;
            // att_gyro = att_fused;
        }

        // todo note: In your current iteration, att fused and att gyro in state are the same.
        self.attitude = att_fused;
        self.heading_gyro = heading_gyro;

        self.att_from_acc = att_acc;
        self.att_from_gyros = att_fused; // todo: QC if this is what you want.
                                         // self.att_from_gyros = att_gyro;

        self.timestamp += self.dt;

        static mut i: u32 = 0;
        unsafe { i += 1 };

        // if unsafe { i } % 1000 == 0 {
        if false {
            // println!("Alignment: {}", acc_gyro_alignment);

            let euler = self.attitude.to_euler();

            let (x_component, y_component, z_component) = att_to_axes(self.attitude);
            let (x_component_acc, y_component_acc, z_component_acc) =
                att_to_axes(self.att_from_acc);
            let (x_component_acc_h, y_component_acc_h, z_component_acc_h) =
                att_to_axes(att_acc_w_heading);
            let (x_component_gyro, y_component_gyro, z_component_gyro) =
                att_to_axes(self.att_from_gyros);

            println!(
                "\n\nAxis rots fused: x{} y{} z{}",
                x_component, y_component, z_component
            );

            println!(
                "Axis rots acc : x{} y{} z{}",
                x_component_acc, y_component_acc, z_component_acc
            );

            println!(
                "Axis rots acc w hdg: x{} y{} z{}",
                x_component_acc_h, y_component_acc_h, z_component_acc_h
            );

            println!(
                "Axis rots gyro: x{} y{} z{}",
                x_component_gyro, y_component_gyro, z_component_gyro
            );

            println!("Euler: p{} r{} y{}", euler.pitch, euler.roll, euler.yaw);
            println!("Acc: x{} y{} z{}", accel_data.x, accel_data.y, accel_data.z);

            println!("Acclen: {}", accel_data.magnitude());

            println!("\nHeading fused: {:?}\n", heading_fused);

            println!(
                "Hdg diff mag gyro: {}",
                self.recent_dh_mag__dh_gyro.unwrap_or(69.)
            );
        }
    }

    /// Attempt to separate linear from gravitational acceleration.
    /// Returns the estimate mof linear acceleration. This is useful for
    /// #1: Determining how much faith (and wight) to put into the accelerometer reading for
    /// attitude determination.
    /// #2: Providing an acc solution that's closer to the true one for this fusing.
    /// #3: Removing linear acceleration when computing position from dead-reckoning
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
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP); // Verified correct

        // An alignment of 1 indicates the estimated up direction between the accelerometer
        // and gyro match. A value of 0 means they're perpendicular.
        // let _acc_gyro_alignment = accel_data.to_normalized().dot(grav_axis_from_att_gyro);

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

        let accel_mag = accel_data.magnitude();

        // For this, we postulate that the gyro's attitude is correct, and therefore the force
        // for gravity is that axis's *up* vector, multiplied by G. The difference between the accelerometer
        // readings and this, therefore is linear acceleration.

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

        // Bias code below; for now, unused.
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

        // todo: Update biases automatically for a short window after bootup.

        let lin_acc_estimate_bias_removed = lin_acc_estimate - lin_acc_bias;

        // If the magntidue of the acceleration is above this value, we are under linear acceleration,
        // and should ignore the accelerometer.
        // let acc_magnitude_thresh_lower = G * self.config.acc_mag_threshold_no_lin.0;
        // let acc_magnitude_thresh_upper = G * self.config.acc_mag_threshold_no_lin.1;

        let mut update_gyro_from_acc = false;
        // If it appears there is negligible linear acceleration, update our gyro readings as appropriate.
        if (accel_mag - G).abs() < self.config.total_accel_thresh {
            // We guess no linear acc since we're getting close to 1G. Note that
            // this will produce false positives in some cases.
            update_gyro_from_acc = true;
            // } else if lin_acc_estimate_bias_removed.magnitude() < lin_acc_thresh {
        } else if lin_acc_estimate.magnitude() < self.config.lin_acc_thresh {
            // If not under much acceleration, re-cage our attitude.
            update_gyro_from_acc = true;
        }

        // todo: alignment fn
        if self.timestamp > self.config.start_alignment_time as f32
            && self.timestamp
                < (self.config.start_alignment_time + self.config.alignment_duration) as f32
        {
            // Identify and remove linear bias.
        }

        static mut i: u32 = 0;
        unsafe { i += 1 };

        // if unsafe { i } % 1000 == 0 {
        if false {
            // println!("Ag: {}", _acc_gyro_alignment);
            println!(
                "Lin bias: x{} y{} z{}",
                lin_acc_bias.x, lin_acc_bias.y, lin_acc_bias.z,
            );

            println!(
                "Lin x{} y{} z{}. mag{}",
                lin_acc_estimate_bias_removed.x,
                lin_acc_estimate_bias_removed.y,
                lin_acc_estimate_bias_removed.z,
                lin_acc_estimate_bias_removed.magnitude() // lin_acc_estimate.x,
                                                          // lin_acc_estimate.y,
                                                          // lin_acc_estimate.z,
                                                          // lin_acc_estimate.magnitude()
            );

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

        // todo: Temporarily not removing bias.
        (update_gyro_from_acc, lin_acc_estimate_bias_removed)
        // (update_gyro_from_acc, lin_acc_estimate)
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

    let sign_z = -axis.z.signum();
    (axis.project_to_vec(UP) * angle).magnitude() * sign_z
}

/// Find the rotation around the Z axis associated with an attitude. We use this to apply
/// gyro and mag heading to the remaining degree of freedom on attitude from acceleration.
fn find_z_rot(heading: f32, accel_norm: Vec3) -> Quaternion {
    // Remove the final degree of freedom using heading. Rotate around UP in the earth frame.
    Quaternion::from_axis_angle(accel_norm, -heading)
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
    mag.y.atan2(mag.x)
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

/// Convert an attitude to rotations around individual axes.
pub fn att_to_axes(att: Quaternion) -> (f32, f32, f32) {
    let axis = att.axis();
    let angle = att.angle();

    let sign_x = -axis.x.signum();
    let sign_y = -axis.y.signum();
    let sign_z = -axis.z.signum();

    let x_component = (axis.project_to_vec(RIGHT) * angle).magnitude() * sign_x;
    let y_component = (axis.project_to_vec(FORWARD) * angle).magnitude() * sign_y;
    let z_component = (axis.project_to_vec(UP) * angle).magnitude() * sign_z;

    (x_component, y_component, z_component)
}
