//! Calculate attitude gven a 9-axis IMU.
//!
//!
//! todo: Use magnetic inclination, declination, and strength

// todo: Calibrate mag and acc based on known gravity value, and known magnetic strength at the current
// todo position, if known.

// todo: We are currently moving away from the AHRS Fusion port.

use num_traits::float::Float; // abs etc

use lin_alg2::f32::{Mat3, Quaternion, Vec3};

// static MAG_CAL_I: AtomicUsize = AtomicUsize::new(0);

use defmt::println;

use crate::{
    blend, linear_acc,
    mag_ellipsoid_fitting::{MAG_SAMPLES_PER_CAT, SAMPLE_VERTICES},
    print_quat, DeviceOrientation, Fix, FORWARD, G, RIGHT, UP,
};

pub struct AhrsConfig {
    /// How far to look back when determining linear acceleration bias from cumulative values. In seconds.
    pub lin_bias_lookback: f32,
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
    /// This should be lower than the amount we update from acc, but enough to overcome
    /// drift, and not be overwhelmed by the AC's update.
    pub update_amt_att_from_mag: f32,
    /// Affects how much gyro biases are updated from the accelerometer-based angular rate
    /// estimation. This should be relatively low, since we don't expect the bias to change much,
    /// and the acc readings are noisy, but stable over time.
    pub update_amt_gyro_bias_from_acc: f32,
    /// If total acclerometer reading is within this value of G, update gyro from acc.
    pub total_accel_thresh: f32, // m/s^2
    /// Similar function as for acc, but comparing to 1. (Our ideal magnetometer magnitude
    /// after calibration must be between 1 + and 1 - this to be considered valid.)
    pub total_mag_thresh: f32,
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
    // A value of 1.0 means new mag cals replace the prev. 0.5 means an average.
    // pub update_amt_mag_cal: f32,
    /// In seconds. If the most recent fix is older than this, don't use it.
    pub max_fix_age_lin_acc: f32,
    pub orientation: DeviceOrientation,
}

impl Default for AhrsConfig {
    fn default() -> Self {
        Self {
            lin_bias_lookback: 10.,
            // mag_gyro_diff_thresh: 0.01,
            update_amt_att_from_acc: 3.,
            update_amt_att_from_mag: 2.,
            update_amt_gyro_bias_from_acc: 0.10,
            total_accel_thresh: 1.0, // m/s^2
            total_mag_thresh: 0.3,   // rel to 1
            lin_acc_thresh: 0.3,     // m/s^2
            start_alignment_time: 2,
            alignment_duration: 2,
            mag_cal_timestep: 0.05,
            update_amt_mag_incl_estimate: 0.05,
            update_ratio_mag_incl: 100,
            update_ratio_mag_cal_log: 160,
            mag_cal_portion_req: 0.65,
            // update_amt_mag_cal: 0.5,
            max_fix_age_lin_acc: 0.5,
            orientation: Default::default(),
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
    pub mag_cal_updated: bool,
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
    pub(crate) acc_len_at_rest: f32,
    /// Used when aligning.
    acc_len_cum: f32,
    gyro_bias_eval_cum: Vec3,
    gyro_bias_eval_num_readings: u32,
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
            mag_cal_updated: false,
            mag_cal_data_up: Default::default(),
            mag_cal_data_fwd: Default::default(),
            mag_sample_i_up: Default::default(),
            mag_sample_i_fwd: Default::default(),
            mag_sample_count_up: Default::default(),
            mag_sample_count_fwd: Default::default(),
            acc_len_cum: 0.,
            acc_len_at_rest: G,
            gyro_bias_eval_cum: Vec3::new_zero(),
            gyro_bias_eval_num_readings: 0,
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
    /// Linear acceleration, as estimated by comparing the fused attitude's "UP"
    /// with that of the accelerometer.
    pub(crate) lin_acc_gyro: Vec3,
    /// Linear acceleration as determined by the GNSS. Note that we get these updates much more
    /// seldom than INS updates; store until ready to use
    pub(crate) lin_acc_gnss: Option<Vec3>,
    pub lin_acc_fused: Vec3,
    /// We set this var upon the first update; this forces the gyro to take a full update from the
    /// accelerometer. Without this, we maay experience strong disagreement between the gyro and acc
    /// at start, since the gyro initializes to level, regardless of actual aircraft attitude.
    pub initialized: bool,
    pub(crate) att_from_acc: Quaternion,
    pub(crate) att_from_mag: Option<Quaternion>,
    pub(crate) acc_calibrated: Vec3,
    pub(crate) gyro_calibrated: Vec3,
    pub(crate) mag_calibrated: Option<Vec3>,
    /// We use these to stored headings to track magnetometer health over time.
    pub(crate) heading_mag: Option<f32>,
    // pub(crate) heading_gyro: f32,
    // /// Estimatated angular rates from the accelerometer-based attitude.
    // pub(crate) acc_rate_estimate: Vec3,
    // linear_acc_confidence: f32, // todo?
    // /// Track recent changing in mag heading, per change in gyro heading, in degrees.
    // /// We use this to assess magnetometer health, ie if it's being interfered with.
    // /// This is `None` if there are no mag reading provided.
    // pub(crate) recent_dh_mag_dh_gyro: Option<f32>,
    /// Use our gyro/acc fused attitude to estimate magnetic inclination.
    pub(crate) mag_inclination_estimate: f32,
    // mag_cal_in_progress: bool,
    /// Time between updates, in seconds.
    pub(crate) dt: f32,
    // /// Radians
    // pub mag_inclination: f32, // todo: Replace with inc estimate above once that's working.
    /// Positive means "east" declination. radians.
    pub(crate) mag_declination: f32,
    /// We use this location for updating linear acceleration using GNSS
    /// todo: You may need a list of several to use the arc-radius approach.
    pub(crate) fix_prev: Option<Fix>,
    /// Long-term difference between accelerometer and gyro readings. Over long periods of time,
    /// the acc-determined angular rate will be reliable due to the fixed reference of gravity.
    /// We use this to zero-out gyro offsets.
    /// We use this to track updates
    pub(crate) num_updates: u32,
    /// Timestamp, in seconds.
    pub(crate) timestamp: f32,
    /// Difference from 1.0 of the magnetometer's magnetude. We use this
    /// to weigh
    pub(crate) recent_mag_variance: f32,
    acc_gyro_rate_diff: Vec3,
    gyro_bias_complete: bool,
    last_fix_timestamp: f32, // seconds, using this struct's timestamp.
}

impl Ahrs {
    pub fn new(dt: f32, orientation: DeviceOrientation) -> Self {
        Self {
            dt,
            mag_declination: -0.032, // Raleigh
            // mag_declination: 6.28 / 8., // todo: testing
            mag_inclination_estimate: 1.2,
            config: AhrsConfig {
                orientation,
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn aligning(&self) -> bool {
        self.timestamp > self.config.start_alignment_time as f32
            && self.timestamp
                < (self.config.start_alignment_time + self.config.alignment_duration) as f32
    }

    /// Update our AHRS solution given new gyroscope, accelerometer, and mag data.
    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        let gyro_calibrated = self.cal.apply_cal_gyro(gyro_data);
        self.gyro_calibrated = gyro_calibrated;

        let mut att_fused = att_from_gyro(gyro_calibrated, self.attitude, self.dt);

        if self.num_updates % ((1. / self.dt) as u32) == 0 {
            print_quat(att_fused, "Att Gyro");
        }

        self.handle_acc(accel_data, &mut att_fused);
        // todo: Temporarily only using acc for lin accel estimate
        self.lin_acc_fused = self.lin_acc_gyro;

        // todo: FIgure out what here should have IIR lowpass filters aplied.

        // Fuse with mag data if available.
        if let Some(mag) = mag_data {
            self.handle_mag(mag, &mut att_fused);
        }

        // These variables here are only used to inspect and debug.
        // let mut acc_rate_estimate = Vec3::new_zero();
        // let mut acc_gyro_rate_diff = Vec3::new_zero();

        // Note that we are only updating gyro biases if under relatively low linear acceleration. We use
        // a lower threshold than the gyro-updates algorithm above, since we don't need to update biases often,
        // so can afford to be pickier.
        // if lin_acc_estimate.magnitude() < self.config.lin_acc_thresh / 2. {
        //     (acc_rate_estimate, acc_gyro_rate_diff) =
        //         self.update_gyro_bias(gyro_data, att_acc, att_acc_prev);
        // }

        // todo: Put this in a fn. (gyro bias estimation at init/rest)
        // todo: Try a simpler gyro bias update here
        // todo: Config vals for these.
        // The lower bound here is to avoid oscillations from physically plugging in the device.
        let min_bias_time = 1.;
        let max_bias_time = 6.;
        let max_bias_val = 0.05; // todo: Adjust this A/R.

        if self.timestamp > min_bias_time
            && self.timestamp < max_bias_time
            && gyro_data.magnitude() < max_bias_val
        {
            self.cal.gyro_bias_eval_cum += gyro_data;
            self.cal.gyro_bias_eval_num_readings += 1;
        } else if self.timestamp > max_bias_time {
            if self.cal.gyro_bias_eval_num_readings > 0 && !self.gyro_bias_complete {
                self.cal.gyro_bias =
                    self.cal.gyro_bias_eval_cum / self.cal.gyro_bias_eval_num_readings as f32;

                self.gyro_bias_complete = true;
            }
        }

        // Time out the GNSS-based lin acc measurements as required
        if self.timestamp - self.last_fix_timestamp > self.config.max_fix_age_lin_acc {
            self.lin_acc_gnss = None;
        }

        // Normalize each update, to prevent normalization errors from accumulating.
        self.attitude = att_fused.to_normalized();

        self.timestamp += self.dt;

        if !self.initialized {
            self.initialized = true;
        }

        if self.num_updates % ((1. / self.dt) as u32) == 0 {
            //     if false {
            // println!("Alignment: {}", acc_gyro_alignment);

            print_quat(self.attitude, "\n\nAtt fused");
            print_quat(self.att_from_acc, "Att acc");

            println!(
                "Gyro raw: x{} y{} z{}. Cal: x{} y{} z{}",
                gyro_data.x,
                gyro_data.y,
                gyro_data.z,
                gyro_calibrated.x,
                gyro_calibrated.y,
                gyro_calibrated.z,
            );

            // println!(
            //     "Acc rate: x{} y{} z{}",
            //     acc_rate_estimate.x, acc_rate_estimate.y, acc_rate_estimate.z,
            // );

            // println!(
            //     "Acc gyro rate diff inst: x{} y{} z{}",
            //     acc_gyro_rate_diff.x, acc_gyro_rate_diff.y, acc_gyro_rate_diff.z,
            // );

            // println!(
            //     "Acc gyro rate diff: x{} y{} z{}",
            //     self.acc_gyro_rate_diff.x, self.acc_gyro_rate_diff.y, self.acc_gyro_rate_diff.z,
            // );

            // println!(
            //     "Acc: x{} y{} z{} mag{}",
            //     acc_calibrated.x,
            //     acc_calibrated.y,
            //     acc_calibrated.z,
            //     acc_calibrated.magnitude()
            // );

            println!(
                "Lin acc gyro x{} y{} z{}",
                self.lin_acc_gyro.x, self.lin_acc_gyro.y, self.lin_acc_gyro.z
            );

            if let Some(la) = self.lin_acc_gnss {
                println!("Lin acc GNSS: x{} y{} z{}", la.x, la.y, la.z);
            }

            print_quat(self.att_from_acc, "ATT from acc");

            // println!("\nHeading fused: {:?}\n", heading_fused);

            // println!("Heading gyro: {}", heading_gyro);

            // println!(
            //     "Hdg diff mag gyro: {}",
            //     self.recent_dh_mag_dh_gyro.unwrap_or(69.)
            // );
        }

        self.num_updates += 1;
    }

    /// Assumes no linear acceleration. Estimates linear acceleration biases.
    fn align(&mut self, lin_acc_estimate: Vec3, acc_data: Vec3) {
        // todo: Rework this.

        if self.aligning() {
            // Maybe this will help with bias estimates before we set it properly.
            self.cal.acc_len_at_rest = acc_data.magnitude();

            self.cal.linear_acc_cum += lin_acc_estimate * self.dt;
            self.cal.acc_len_cum += acc_data.magnitude() * self.dt;
            // self.cal.linear_acc_cum += lin_acc_estimate;
        }

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
    fn _update_gyro_bias(
        &mut self,
        gyro_raw: Vec3,
        att_acc: Quaternion,
        att_from_acc_prev: Quaternion,
    ) -> (Vec3, Vec3) {
        let d_att_acc = att_acc / att_from_acc_prev;
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

    /// Update the linear acc estimates and related state from a GNSS fix.
    pub fn update_from_fix(&mut self, fix: &Fix) {
        self.last_fix_timestamp = self.timestamp;

        // println!("FIX here: {}", fix.timestamp_s);
        if let Some(fix_prev) = &self.fix_prev {
            if fix.timestamp_s - fix_prev.timestamp_s < self.config.max_fix_age_lin_acc {
                self.lin_acc_gnss = Some(linear_acc::from_gnss(fix, fix_prev, self.attitude));
            }
        }

        self.fix_prev = Some(fix.clone());
        // todo: Perhaps
        // let lin_acc_gnd_track = linear_acc::from_ground_track();
    }
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

/// gnss_acc: in m/s; differed from subsequent readings.
pub fn heading_from_gnss_acc(gnss_acc_nse: Vec3, acc_lin_imu: Vec3) -> f32 {
    // Note that we only care about the relative directions; magnitude of acceleration isn't required.

    // let acc_lin_nse = Quaternion::from_unit_vecs(
    //     acc_lin_imu.to_normalized(),
    //     gnss_acc_nse.to_normalized()
    // );

    // todo: Does not appear to be working
    let att_from_gnss =
        Quaternion::from_unit_vecs(gnss_acc_nse.to_normalized(), acc_lin_imu.to_normalized());

    heading_from_att(att_from_gnss)
}

/// Used to nudge the gyro from an absolute, but incomplete attitude reference, eg accelerometer
/// or magnetometer. `ref_vec` is a vector that we know is correct from the updating source.
/// This rotation is heading-invariant: Eg Rotate the gyro *up* towards the acc *up*.
pub(crate) fn make_nudge(
    attitude: Quaternion,
    sensor_norm: Vec3,
    anchor: Vec3,
    update_amt: f32,
) -> Quaternion {
    let anchor_ac_frame = attitude.rotate_vec(anchor);

    let rotation = Quaternion::from_unit_vecs(anchor_ac_frame, sensor_norm);

    Quaternion::new_identity().slerp(rotation, update_amt)
}
