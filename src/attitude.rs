//! Calculate attitude gven a 9-axis IMU.
//!
//!
//! todo: Use magnetic inclination, declination, and strength

// todo: Calibrate mag and acc based on known gravity value, and known magnetic strength at the current
// todo position, if known.

// todo: We are currently moving away from the AHRS Fusion port.

use core::{
    f32::consts::TAU,
    sync::atomic::{AtomicUsize, Ordering},
};

use num_traits::float::Float; // abs etc

use lin_alg2::f32::{Mat3, Quaternion, Vec3};

// static MAG_CAL_I: AtomicUsize = AtomicUsize::new(0);

use crate::{
    mag_ellipsoid_fitting::{self, MAG_SAMPLES_PER_CAT, SAMPLE_VERTEX_ANGLE, SAMPLE_VERTICES},
    print_quat, FORWARD, G, RIGHT, UP,
};

use defmt::{println, write};

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
    /// will occur. Higher values lead to more bias from linear acceleration. Low values lead to more
    /// gyro drift. Values from 0.1 to 10 may be optimal.
    ///
    /// This value can be thought of as the 1 / the number of seconds to correct a gyro reading to match the
    /// accelerometer. Keep this in mind re expectations of gyro drift rate. It must be high enough to
    /// compensate for drift (but perhaps not much higher)
    pub update_amt_att_from_acc: f32,
    /// Affects how much gyro biases are updated from the accelerometer-based angular rate
    /// estimation. This should be relatively low, since we don't expect the bias to change much,
    /// and the acc readings are noisy, but stable over time.
    pub update_amt_gyro_bias_from_acc: f32,
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
    // pub calibration: crate::ImuCalibration,
    /// Time, in seconds, after initialization, to start alignment procedure. This is set up so
    /// powering the device on by plugging it in, etc, doesn't interfere.
    pub start_alignment_time: u8,
    /// Time, in seconds, of the alignment.
    pub alignment_duration: u8,
    /// Time in seconds between each magnetometer calibration reading.
    pub mag_cal_timestep: f32,
    /// How much to update the estimated magnetic inclination. See rule of thumb for the other
    /// update rates above. This should be very low, as this should be constant for a given part
    /// of the earth.
    pub update_amt_mag_incl_estimate: f32,
    /// Update the mag incl estimate at this ratio of overall updates.
    pub update_ratio_mag_incl: u16,
    /// Log a mag cal point (or skip logging if that sector already has enough data) every this
    /// many updates. Setting too high may have performance consequences (although probably not a factor),
    /// and cause sector categories to fill with points close to each other. Setting too low many
    /// cause delayed calibration.
    pub update_ratio_mag_cal_log: u16,
    /// This portion of 1 categories must be filled to initiate a calibration.
    pub mag_cal_portion_req: f32,
    /// A value of 1.0 means new mag cals replace the prev. 0.5 means an average.
    pub mag_cal_update_amt: f32,
}

impl Default for AhrsConfig {
    fn default() -> Self {
        Self {
            lin_bias_lookback: 10.,
            mag_diff_lookback: 10.,
            mag_gyro_diff_thresh: 0.01,
            update_amt_att_from_acc: 3.,
            update_amt_gyro_bias_from_acc: 0.01,
            update_port_mag_heading: 0.1,
            // acc_mag_threshold_no_lin: (0.8, 1.2),
            total_accel_thresh: 0.4, // m/s^2
            lin_acc_thresh: 0.3,     // m/s^2
            // calibration: Default::default(),
            start_alignment_time: 2,
            alignment_duration: 2,
            mag_cal_timestep: 0.05,
            update_amt_mag_incl_estimate: 0.05,
            update_ratio_mag_incl: 100,
            update_ratio_mag_cal_log: 160,
            mag_cal_portion_req: 0.85,
            mag_cal_update_amt: 0.3,
        }
    }
}

pub struct AhrsCal {
    /// Subtract this from readings to calibrate.
    pub acc_bias: Vec3,
    pub acc_slope: Vec3,
    /// Subtract this from readings to calibrate.
    pub gyro_bias: Vec3,
    /// linear acceleration, per axis, integrated over time. We use this to
    /// remove biases, under the assumption that this should average out to 0, along
    /// each axis.
    pub linear_acc_cum: Vec3,
    /// Bias, determined from the alignment process.
    pub linear_acc_bias: Vec3,
    pub hard_iron: Vec3,
    pub soft_iron: Mat3,
    /// Magenetometer cal data, per attitude category.
    pub(crate) mag_cal_data_up: [[Vec3; MAG_SAMPLES_PER_CAT]; SAMPLE_VERTICES.len()],
    pub(crate) mag_cal_data_fwd: [[Vec3; MAG_SAMPLES_PER_CAT]; SAMPLE_VERTICES.len()],
    /// Current mag sample index, per attitude category.
    pub(crate) mag_sample_i_up: [usize; SAMPLE_VERTICES.len()],
    pub(crate) mag_sample_i_fwd: [usize; SAMPLE_VERTICES.len()],
    /// Number of samples taken since last calibration, per attitude category.
    pub(crate) mag_sample_count_up: [u8; SAMPLE_VERTICES.len()],
    pub(crate) mag_sample_count_fwd: [u8; SAMPLE_VERTICES.len()],
    /// In m/s^2. Used for determining linear acceleration. This should be close to G.
    acc_len_at_rest: f32,
    /// Used when aligning.
    acc_len_cum: f32,
}

impl Default for AhrsCal {
    fn default() -> Self {
        Self {
            acc_bias: Vec3::new_zero(),
            acc_slope: Vec3::new(1., 1., 1.),
            gyro_bias: Vec3::new_zero(),
            linear_acc_cum: Vec3::new_zero(),
            linear_acc_bias: Vec3::new_zero(),
            // hard_iron: Vec3::new_zero(),
            // Rough, from GPS mag can.
            hard_iron: Vec3::new_zero(),
            soft_iron: Mat3::new_identity(),
            mag_cal_data_up: Default::default(),
            mag_cal_data_fwd: Default::default(),
            mag_sample_i_up: Default::default(),
            mag_sample_i_fwd: Default::default(),
            mag_sample_count_up: Default::default(),
            mag_sample_count_fwd: Default::default(),
            acc_len_cum: 0.,
            acc_len_at_rest: G,
        }
    }
}

impl AhrsCal {
    /// Run this when the device is stationary on a flat surface, with the Z axis up,
    /// to initiate acceleratometer calibration. Updates intercepts only. Readings are in m/s.
    pub fn calibrate_accel(&mut self, acc_data: Vec3) {
        self.acc_bias.x = acc_data.x;
        self.acc_bias.y = acc_data.y;
        self.acc_bias.z = acc_data.z - G;
    }

    fn apply_cal_acc(&self, data: Vec3) -> Vec3 {
        Vec3::new(
            data.x * self.acc_slope.x - self.acc_bias.x,
            data.y * self.acc_slope.y - self.acc_bias.y,
            data.z * self.acc_slope.z - self.acc_bias.z,
        )
    }

    fn apply_cal_gyro(&self, data: Vec3) -> Vec3 {
        Vec3::new(
            // todo: Slope?
            data.x - self.gyro_bias.x,
            data.y - self.gyro_bias.y,
            data.z - self.gyro_bias.z,
        )
    }
}

#[derive(Default)]
pub struct Ahrs {
    pub config: AhrsConfig,
    pub cal: AhrsCal,
    pub attitude: Quaternion,
    pub linear_acc_estimate: Vec3,
    /// We set this var upon the first update; this forces the gyro to take a full update from the
    /// accelerometer. Without this, we maay experience strong disagreement between the gyro and acc
    /// at start, since the gyro initializes to level, regardless of actual aircraft attitude.
    pub initialized: bool,
    pub(crate) att_from_gyros: Quaternion,
    pub(crate) att_from_acc: Quaternion,
    pub(crate) att_from_mag: Option<Quaternion>,
    pub(crate) acc_calibrated: Vec3,
    pub(crate) gyro_calibrated: Vec3,
    pub(crate) mag_calibrated: Option<Vec3>,
    /// We use these to stored headings to track magnetometer health over time.
    pub(crate) heading_mag: Option<f32>,
    pub(crate) heading_gyro: f32,
    // /// Estimatated angular rates from the accelerometer-based attitude.
    // pub(crate) acc_rate_estimate: Vec3,
    // linear_acc_confidence: f32, // todo?
    /// Track recent changing in mag heading, per change in gyro heading, in degrees.
    /// We use this to assess magnetometer health, ie if it's being interfered with.
    /// This is `None` if there are no mag reading provided.
    pub(crate) recent_dh_mag_dh_gyro: Option<f32>,
    /// Use our gyro/acc fused attitude to estimate magnetic inclination.
    pub(crate) mag_inclination_estimate: f32,
    // mag_cal_in_progress: bool,
    /// Time between updates, in seconds.
    pub(crate) dt: f32,
    // /// Radians
    // pub mag_inclination: f32, // todo: Replace with inc estimate above once that's working.
    /// Positive means "east" declination. radians.
    pub(crate) mag_declination: f32,
    /// Long-term difference between accelerometer and gyro readings. Over long periods of time,
    /// the acc-determined angular rate will be reliable due to the fixed reference of gravity.
    /// We use this to zero-out gyro offsets.
    acc_gyro_rate_diff: Vec3,
    /// Timestamp, in seconds.
    timestamp: f32,
}

impl Ahrs {
    pub fn new(dt: f32) -> Self {
        Self {
            dt,
            mag_declination: -0.032, // Raleigh
            mag_inclination_estimate: -1.2,
            ..Default::default()
        }
    }

    /// Update our AHRS solution given new gyroscope, accelerometer, and mag data.
    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        let acc_calibrated = self.cal.apply_cal_acc(accel_data);
        let gyro_calibrated = self.cal.apply_cal_gyro(gyro_data);

        // todo: A system where any 2/3 of your sources can, in agreement, update the biases
        // todo of the third.

        // todo: FIgure out what here should have IIR lowpass filters aplied.

        let accel_norm = acc_calibrated.to_normalized();

        // Estimate attitude from raw accelerometer and gyro data. Note that
        // The gyro data reguarly receives updates from the acc and mag.
        let att_acc = att_from_accel(accel_norm);

        let att_acc_prev = self.att_from_acc;
        self.att_from_acc = att_acc;

        let mut att_gyro = att_from_gyro(gyro_calibrated, self.att_from_gyros, self.dt);

        let heading_gyro = heading_from_att(att_gyro);

        // todo: Remove `heading_gyro` if you end up not using it.
        self.heading_gyro = heading_gyro;

        // See comment on the `initialized` field.
        // We update initialized state at the end of this function, since other steps rely on it.
        if !self.initialized {
            self.att_from_gyros = att_acc;
            att_gyro = self.att_from_gyros;
        }

        let mut heading_fused = heading_gyro;

        // Fuse with mag data if available.
        match mag_data {
            Some(mut mag) => {
                self.handle_mag(mag, heading_gyro, unsafe { I });
            }
            None => {
                self.recent_dh_mag_dh_gyro = None;
            }
        }

        let (update_gyro_from_acc, lin_acc_estimate) =
            self.handle_linear_acc(acc_calibrated, att_gyro);

        let att_acc_w_lin_removed =
            att_from_accel((acc_calibrated - lin_acc_estimate).to_normalized());

        let mut att_fused = att_gyro;

        // todo: Instead of a binary update-or-not, consider weighing the slerp value based
        // todo on how much lin acc we assess, or how much uncertainly in lin acc.

        // These variables here are only used to inspect and debug.
        let mut acc_rate_estimate = Vec3::new_zero();
        let mut acc_gyro_rate_diff = Vec3::new_zero();

        if update_gyro_from_acc {
            // println!("TRUE");
            // Apply a rotation of the gyro solution towards the acc solution, if we think we are not under
            // much linear acceleration.
            // This rotation is heading-invariant: Rotate the gyro *up* towards the acc *up*.
            let gyro_up = att_gyro.rotate_vec(UP);
            let rot_gyro_to_acc = Quaternion::from_unit_vecs(gyro_up, accel_norm);
            // let rot_gyro_to_acc = Quaternion::from_unit_vecs(accel_norm, gyro_up);

            let rot_to_apply_to_gyro = Quaternion::new_identity().slerp(
                rot_gyro_to_acc,
                self.config.update_amt_att_from_acc * self.dt,
            );

            // if unsafe { I } % 1000 == 0 {
            if false {
                print_quat(rot_gyro_to_acc, "Rot to apply");
                println!("rot angle: {}", rot_gyro_to_acc.angle());
            }

            att_fused = rot_to_apply_to_gyro * att_gyro;

            // Note that we are only updating gyro biases if under relatively low linear acceleration.
            (acc_rate_estimate, acc_gyro_rate_diff) =
                self.update_gyro_bias(gyro_data, att_acc, att_acc_prev);
        }

        // todo note: In your current iteration, att fused and att gyro in state are the same.
        self.attitude = att_fused;

        self.att_from_acc = att_acc;
        self.att_from_gyros = att_fused; // todo: QC if this is what you want.
                                         // self.att_from_gyros = att_gyro;

        self.timestamp += self.dt;

        if !self.initialized {
            self.initialized = true;
        }

        static mut I: u32 = 0;
        unsafe { I += 1 };

        if unsafe { I } % 1000 == 0 {
            // if false {
            // println!("Alignment: {}", acc_gyro_alignment);

            // let euler = self.attitude.to_euler();
            // println!("Euler: p{} r{} y{}", euler.pitch, euler.roll, euler.yaw);

            // print_quat(self.attitude, "\n\nAtt fused");

            // print_quat(self.att_from_acc, "Att Acc");

            // Temp: offset at idle appears to be -0.015, -0.01, .004
            println!(
                "\nGyro raw: x{} y{} z{}",
                gyro_data.x, gyro_data.y, gyro_data.z,
            );

            println!(
                "Acc rate: x{} y{} z{}",
                acc_rate_estimate.x, acc_rate_estimate.y, acc_rate_estimate.z,
            );

            println!(
                "Acc gyro rate diff inst: x{} y{} z{}",
                acc_gyro_rate_diff.x, acc_gyro_rate_diff.y, acc_gyro_rate_diff.z,
            );

            println!(
                "Acc gyro rate diff: x{} y{} z{}",
                self.acc_gyro_rate_diff.x, self.acc_gyro_rate_diff.y, self.acc_gyro_rate_diff.z,
            );

            println!(
                "Gyro Cal: x{} y{} z{}\n",
                gyro_calibrated.x, gyro_calibrated.y, gyro_calibrated.z,
            );

            // print_quat(self.att_from_gyros, "Att gyros");

            // println!(
            //     "Acc: x{} y{} z{} mag{}",
            //     acc_calibrated.x,
            //     acc_calibrated.y,
            //     acc_calibrated.z,
            //     acc_calibrated.magnitude()
            // );

            // println!("\nHeading fused: {:?}\n", heading_fused);

            // println!("Heading gyro: {}", heading_gyro);

            println!("Acc len at rest: {:?}", self.cal.acc_len_at_rest);

            println!(
                "Hdg diff mag gyro: {}",
                self.recent_dh_mag_dh_gyro.unwrap_or(69.)
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
        // att_acc: Quaternion,
        // // att gyro must have heading removed, or acc must have heading added.
        att_gyro: Quaternion,
    ) -> (bool, Vec3) {
        //  Ways to identify linear acceleration:
        // - Greater or less than 1G of acceleration, if the accel is calibrated.
        // - Discontinuities or other anomolies when integrating accel-based attitude over time,
        // - or, along those lines, discontinuities etc when fusing with gyro.

        // This is the up vector as assessed from the attitude from the gyro. It is equivalent to
        // the accelerometer's normalized vector when linear acceleration is 0.
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP);

        // Estimate linear acceleration by comparing the accelerometer's normalized vector (indicating the
        // direction it resists gravity) with that estimated from the gyro. This is a proxy for linear
        // acceleration.

        let accel_mag = accel_data.magnitude();

        // For this, we postulate that the gyro's attitude is correct, and therefore the force
        // for gravity is that axis's *up* vector, multiplied by G. The difference between the accelerometer
        // readings and this, therefore is linear acceleration.
        // This linear acc estimate is in earth coords.

        // todo: Project to elim heading effects?
        let att_acc_non_norm = att_from_accel(accel_data.to_normalized());
        let att_diff_rot = att_acc_non_norm * att_gyro.inverse();

        // acc = lin + grav
        // lin = acc - (grav_axis_gyro * G)
        // For the purpose of this calculation, we are assuming the real gravitation axis is
        // that determined by the gyro.
        // This is in the aircraft's frame of reference.
        let lin_acc_estimate = accel_data - (grav_axis_from_att_gyro * self.cal.acc_len_at_rest);

        self.align(lin_acc_estimate, accel_data);

        // Store our linear acc estimate and accumulator before compensating for bias.
        self.linear_acc_estimate = lin_acc_estimate; // todo: DOn't take all of it; fuse with current value.
                                                     // todo: Be careful about floating point errors over time.
                                                     // todo: Toss extreme values?
                                                     // todo: Lowpass?

        let lin_acc_estimate_bias_removed = lin_acc_estimate - self.cal.linear_acc_bias;

        let mut update_gyro_from_acc = false;
        // If it appears there is negligible linear acceleration, update our gyro readings as appropriate.
        if (accel_mag - self.cal.acc_len_at_rest).abs() < self.config.total_accel_thresh {
            // We guess no linear acc since we're getting close to 1G. Note that
            // this will produce false positives in some cases.
            update_gyro_from_acc = true;
        } else if lin_acc_estimate_bias_removed.magnitude() < self.config.lin_acc_thresh {
            // If not under much acceleration, re-cage our attitude.
            update_gyro_from_acc = true;
        }

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
            let a = grav_axis_from_att_gyro * self.cal.acc_len_at_rest;
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

        (update_gyro_from_acc, lin_acc_estimate_bias_removed)
    }

    /// Assumes no linear acceleration. Estimates linear acceleration biases.
    fn align(&mut self, lin_acc_estimate: Vec3, acc_data: Vec3) {
        if self.timestamp > self.config.start_alignment_time as f32
            && self.timestamp
                < (self.config.start_alignment_time + self.config.alignment_duration) as f32
        {
            // Maybe this will help with bias estimates before we set it properly.
            self.cal.acc_len_at_rest = acc_data.magnitude();

            self.cal.linear_acc_cum += lin_acc_estimate * self.dt;
            self.cal.acc_len_cum += acc_data.magnitude() * self.dt;
            // self.cal.linear_acc_cum += lin_acc_estimate;
        }

        // todo: Rework this.
        if self.cal.linear_acc_bias == Vec3::new_zero()
            && self.timestamp
                > (self.config.start_alignment_time + self.config.alignment_duration) as f32
        {
            self.cal.linear_acc_bias =
                self.cal.linear_acc_cum / self.config.alignment_duration as f32;

            self.cal.acc_len_at_rest = self.cal.acc_len_cum / self.config.alignment_duration as f32;

            println!("\n\nAlignment complete \n\n");
        }
    }

    /// Update gyro bias from accelerometer-determined angular rate.
    /// Returns the rate estimate from accelerometer, and the instantaneous difference from
    /// gyro rate, for use in inspection and debugging.
    fn update_gyro_bias(
        &mut self,
        gyro_raw: Vec3,
        att_acc: Quaternion,
        att_from_acc_prev: Quaternion,
    ) -> (Vec3, Vec3) {
        let d_att_acc = att_acc * att_from_acc_prev.inverse();
        let d_att_axes = d_att_acc.to_axes();

        // We notice that, at least when there is little linear acc, this value bounces around
        // the correct one. It generally has a higher error magnitude than the gyros, but perhaps
        // it's average is better.
        // let acc_rate_estimate = Vec3::new(d_att_acc.x, d_att_acc.y, d_att_acc.z) / self.dt;
        let acc_rate_estimate = Vec3::new(d_att_axes.0, d_att_axes.1, d_att_axes.2) / self.dt;

        let acc_gyro_rate_diff = acc_rate_estimate - gyro_raw;

        let update_amt = self.config.update_amt_gyro_bias_from_acc * self.dt;
        let update_amt_inv = 1. - update_amt;

        self.acc_gyro_rate_diff =
            self.acc_gyro_rate_diff * update_amt_inv + acc_gyro_rate_diff * update_amt;

        self.cal.gyro_bias = -self.acc_gyro_rate_diff; // todo temp: Something more sophisticated?

        (acc_rate_estimate, acc_gyro_rate_diff)
    }
}

/// Estimate attitude from accelerometer. This will fail when under
/// linear acceleration. Apply calibration prior to this step.
/// Uses the previous attitude to rotate along the remaining degree of freedom (heading)
pub fn att_from_accel(accel_norm: Vec3) -> Quaternion {
    Quaternion::from_unit_vecs(UP, accel_norm)
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

/// Calculate heading, from an attitude;
fn heading_from_att(att: Quaternion) -> f32 {
    let axis = att.axis();
    let angle = att.angle();

    let sign_z = -axis.z.signum();
    (axis.project_to_vec(UP) * angle).magnitude() * sign_z
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

// todo: For mag cal, you can use acc/gyro-fused attitude to determine which mag sample points to keep,
// todo ie to keep points evenly-spaced. Try that this evening. T

/// Estimate hard and soft iron offsets from sample points
/// todo: Run this periodically/regularly, instead of a separate mag cal procedure.
fn mag_offsets_from_points(sample_pts: &[Vec3]) -> (Mat3, Vec3) {
    // Least-squares approach to model the ellipse.
    // http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html

    // Cruder approach using max and min values: https://www.appelsiini.net/2018/calibrate-magnetometer/

    // For now, our hard-iron procedure uses max and min values. We should at least make some effort
    // to remove outliers.

    // see also: http://www.juddzone.com/ALGORITHMS/least_squares_precision_3D_ellipsoid.html

    let mut x_min = 99999.;
    let mut x_max = -99999.;
    let mut y_min = 99999.;
    let mut y_max = -99999.;
    let mut z_min = 99999.;
    let mut z_max = -99999.;

    for pt in sample_pts {
        if pt.x < x_min {
            x_min = pt.x;
        } else if pt.x > x_max {
            x_max = pt.x;
        }
        if pt.y < y_min {
            y_min = pt.y;
        } else if pt.y > y_max {
            y_max = pt.y;
        }

        if pt.z < z_min {
            z_min = pt.z;
        } else if pt.z > z_max {
            z_max = pt.z;
        }
    }

    let hard_iron = Vec3::new(
        (x_max - x_min) / 2.,
        (y_max - y_min) / 2.,
        (z_max - z_min) / 2.,
    );

    let avg_delta = (hard_iron.x + hard_iron.y + hard_iron.z) / 3.;

    let scale_x = avg_delta / hard_iron.x;
    let scale_y = avg_delta / hard_iron.y;
    let scale_z = avg_delta / hard_iron.z;

    // todo temp
    let soft_iron = Mat3::new_identity();

    (soft_iron, hard_iron)
}
