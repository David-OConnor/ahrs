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
    mag_ellipsoid_fitting::{
        self, MAG_SAMPLES_PER_CAT, SAMPLE_VERTEX_ANGLE, SAMPLE_VERTICES, TOTAL_MAG_SAMPLE_PTS,
    },
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
}

impl Default for AhrsConfig {
    fn default() -> Self {
        Self {
            lin_bias_lookback: 10.,
            mag_diff_lookback: 10.,
            mag_gyro_diff_thresh: 0.01,
            update_from_acc_amt: 2.,
            update_port_mag_heading: 0.1,
            // acc_mag_threshold_no_lin: (0.8, 1.2),
            total_accel_thresh: 0.4, // m/s^2
            lin_acc_thresh: 0.4,     // m/s^2
            // calibration: Default::default(),
            start_alignment_time: 2,
            alignment_duration: 2,
            mag_cal_timestep: 0.05,
            update_amt_mag_incl_estimate: 0.05,
            update_ratio_mag_incl: 100,
            update_ratio_mag_cal_log: 50,
        }
    }
}

pub struct AhrsCal {
    pub acc_bias: Vec3,
    pub acc_slope: Vec3,
    /// linear acceleration, per axis, integrated over time. We use this to
    /// remove biases, under the assumption that this should average out to 0, along
    /// each axis.
    pub linear_acc_cum: Vec3,
    /// Bias, determined from the alignment process.
    pub linear_acc_bias: Vec3,
    pub hard_iron: Vec3,
    pub soft_iron: Mat3,
    /// Magenetometer cal data, per attitude category.
    pub mag_cal_data: [[Vec3; MAG_SAMPLES_PER_CAT]; SAMPLE_VERTICES.len()],
    /// Current mag sample index, per attitude category.
    mag_sample_i: [usize; SAMPLE_VERTICES.len()],
    /// Number of samples taken since last calibration, per attitude category.
    mag_sample_count: [u8; SAMPLE_VERTICES.len()],
    /// In m/s^2. Used for determining linear acceleration. This should be close to G.
    acc_len_at_rest: f32,
    /// Used when aligning.
    acc_len_cum: f32,
    // mag_cal_state: MagCalState,
}

impl Default for AhrsCal {
    fn default() -> Self {
        Self {
            acc_bias: Vec3::new_zero(),
            acc_slope: Vec3::new(1., 1., 1.),
            linear_acc_cum: Vec3::new_zero(),
            linear_acc_bias: Vec3::new_zero(),
            // hard_iron: Vec3::new_zero(),
            // Rough, from GPS mag can.
            hard_iron: Vec3::new_zero(),
            soft_iron: Mat3::new_identity(),
            mag_cal_data: [[Vec3::new_zero(); MAG_SAMPLES_PER_CAT]; SAMPLE_VERTICES.len()],
            mag_sample_i: Default::default(),
            mag_sample_count: Default::default(),
            acc_len_cum: 0.,
            acc_len_at_rest: G,
            // mag_cal_state: Default::default(),
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

    /// Apply the hard and soft iron offsets to our readings.
    pub fn apply_mag_cal(&self, mag_data: Vec3) -> Vec3 {
        // todo: Why this clone?
        self.soft_iron.clone() * (mag_data - self.hard_iron)
    }

    // /// Estimate hard iron offset from the filled data buffer.
    // /// todo: How should we do this? Best-fit ellipse, then find its center?
    // /// todo: To keep your buffers from running you out of flash, consider a few independent estimates,
    // /// todo then averaging (etc) them.
    // pub fn estimate_hard_iron(&mut self) {
    //     // self.hard_iron =
    // }

    /// Run this once every x updates.
    /// Stores a calibration point into a directional category. The aim is to evenly mix samples taked
    /// at ~evenly-spaced attitudes, despite attitudes not being distributed evenly over time.
    /// This collects readings, then performs calibration once enough are taken in each category.
    pub fn log_mag_cal_pt(&mut self, att: Quaternion, mag_raw: Vec3) {
        if mag_raw == Vec3::new_zero() {
            return; // todo: Not sure why we occasionally (and at start) get this.
        }

        // Indexed by the sample vec vertices.
        let mut sample_category = 0;
        for (cat, sample) in SAMPLE_VERTICES.iter().enumerate() {
            // `UP` is arbitrary; any Vec will do as long as we're consistent.

            // todo: This still leaves you with an ambiguity for teh Z readings. You might need 2
            // todo categories of different rotation axes.

            let up_rotated = att.rotate_vec(UP);
            if (up_rotated.dot(*sample)).acos() < SAMPLE_VERTEX_ANGLE {
                sample_category = cat;
                break;
            }
        }
        self.mag_cal_data[sample_category][self.mag_sample_i[sample_category]] = mag_raw;

        // This loops around; we always update the oldest value in each category
        self.mag_sample_i[sample_category] =
            (self.mag_sample_i[sample_category] + 1) % MAG_SAMPLES_PER_CAT;

        if self.mag_sample_count[sample_category] < MAG_SAMPLES_PER_CAT as u8 {
            self.mag_sample_count[sample_category] += 1;
        }

        // To display status.
        if false {
            let mut num_pts_left = 0;
            for cat in self.mag_sample_count {
                num_pts_left += MAG_SAMPLES_PER_CAT as u8 - cat;
            }
            println!("Logging mag pt. Num left: {}", num_pts_left);
        }
    }

    /// Update mag calibration based on sample points.
    /// Only run this when there is recent data of each attitude category of points.
    pub fn update_mag_cal(&mut self) {
        // Flatten our sample data organized by points; that format is only used to ensure a
        // good selection of points, ie fairly weighted on all sides.
        let mut sample_pts = [Vec3::new_zero(); TOTAL_MAG_SAMPLE_PTS];

        let mut i = 0;
        println!("\n\n\n[");
        for cat in self.mag_cal_data {
            for point in cat {
                sample_pts[i] = point;
                println!("({}, {}, {}),", point.x, point.y, point.z);
                i += 1;
            }
        }
        println!("]");

        // todo: Put back once working.

        let poly_terms = mag_ellipsoid_fitting::ls_ellipsoid(&sample_pts);
        let (hard_iron, soft_iron) = mag_ellipsoid_fitting::poly_to_params_3d(&poly_terms);
        // // todo: axes and inve; what are they? Check the web page.
        // self.hard_iron = center;
        // self.soft_iron = inve;

        println!(
            "Hard iron: x{} y{} z{}",
            hard_iron.x, hard_iron.y, hard_iron.z
        );

        println!(
            "Soft iron diag: {} {} {} {} {} {}",
            soft_iron.data[0],
            soft_iron.data[1],
            soft_iron.data[2],
            soft_iron.data[4],
            soft_iron.data[5],
            soft_iron.data[8]
        );

        // Reset our sample data.
        self.mag_cal_data = Default::default();
        self.mag_sample_i = Default::default();
        self.mag_sample_count = Default::default();
    }
}

pub struct Ahrs {
        pub config: AhrsConfig,
        pub cal: AhrsCal,
    pub attitude: Quaternion,
    att_from_gyros: Quaternion,
    att_from_acc: Quaternion,
    att_from_mag: Option<Quaternion>,
    /// We use these to stored headings to track magnetometer health over time.
    heading_mag: Option<f32>,
    heading_gyro: f32,
    pub linear_acc_estimate: Vec3,
    // linear_acc_confidence: f32, // todo?
    /// Track recent changing in mag heading, per change in gyro heading, in degrees.
    /// We use this to assess magnetometer health, ie if it's being interfered with.
    /// This is `None` if there are no mag reading provided.
    recent_dh_mag_dh_gyro: Option<f32>,
    /// Use our gyro/acc fused attitude to estimate magnetic inclination.
    mag_inclination_estimate: f32,
    // mag_cal_in_progress: bool,
    /// Timestamp, in seconds.
    timestamp: f32,
    /// Time between updates, in seconds.
    dt: f32,
    // /// Radians
    // pub mag_inclination: f32, // todo: Replace with inc estimate above once that's working.
    /// Positive means "east" declination. radians.
    pub mag_declination: f32,
    /// We set this var upon the first update; this forces the gyro to take a full update from the
    /// accelerometer. Without this, we maay experience strong disagreement between the gyro and acc
    /// at start, since the gyro initializes to level, regardless of actual aircraft attitude.
    initialized: bool,
}

impl Ahrs {
    pub fn new(dt: f32) -> Self {
        Self {
                config: AhrsConfig::default(),
            cal: AhrsCal::default(),
            attitude: Quaternion::new_identity(),
            att_from_gyros: Quaternion::new_identity(),
            att_from_acc: Quaternion::new_identity(),
            att_from_mag: None,
            heading_mag: None,
            heading_gyro: 0.,
            linear_acc_estimate: Vec3::new_zero(),
            recent_dh_mag_dh_gyro: None,
            mag_inclination_estimate: -1.2,
            // mag_cal_in_progress: true, // todo temp true
            timestamp: 0.,
            dt,
            // mag_inclination: -1.2,    // Rough
            mag_declination: -0.032, // Raleigh
            initialized: false,
        }
    }

    /// Update our AHRS solution given new gyroscope, accelerometer, and mag data.
    pub fn update(&mut self, gyro_data: Vec3, accel_data: Vec3, mag_data: Option<Vec3>) {
        // todo: A system where any 2/3 of your sources can, in agreement, update the biases
        // todo of the third.
        let accel_norm = accel_data.to_normalized();

        // Estimate attitude from raw accelerometer and gyro data. Note that
        // The gyro data reguarly receives updates from the acc and mag.
        let att_acc = att_from_accel(accel_norm);
        self.att_from_acc = att_acc;

        let mut att_gyro = att_from_gyro(gyro_data, self.att_from_gyros, self.dt);

        let heading_gyro = heading_from_att(att_gyro);

        // todo: Remove `heading_gyro` if you end up not using it.
        self.heading_gyro = heading_gyro;

        // todo: Sort these att acc w heading and gyro without ones.

        // Remove the heading component from the gyroscope. This allows us to compare
        // it to the accelerometer to estimate linear acceleration.
        // let att_gyro_without_heading = Quaternion::from_axis_angle(UP, heading_gyro) * att_gyro;

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

        // Consider this flow: mag updates gyro heading. Att updates gyro pitch and roll.
        // Gyro turns into fused?

        // todo: Which approach: Apply heading to acc, or remove from gyro?
        let (update_gyro_from_acc, lin_acc_estimate) =
            // self.handle_linear_acc(accel_data,  att_gyro_without_heading);
            self.handle_linear_acc(accel_data,  att_gyro);

        let att_acc_w_lin_removed = att_from_accel((accel_data - lin_acc_estimate).to_normalized());

        let mut att_fused = att_gyro;

        // todo: Instead of a binary update-or-not, consider weighing the slerp value based
        // todo on how much lin acc we assess, or how much uncertainly in lin acc.
        if update_gyro_from_acc {
            // println!("TRUE");
            // Apply a rotation of the gyro solution towards the acc solution, if we think we are not under
            // much linear acceleration.
            // This rotation is heading-invariant: Rotate the gyro *up* towards the acc *up*.
            let gyro_up = att_gyro.rotate_vec(UP);
            let rot_gyro_to_acc = Quaternion::from_unit_vecs(gyro_up, accel_norm);
            // let rot_gyro_to_acc = Quaternion::from_unit_vecs(accel_norm, gyro_up);

            let rot_to_apply_to_gyro = Quaternion::new_identity()
                .slerp(rot_gyro_to_acc, self.config.update_from_acc_amt * self.dt);

            if unsafe { I } % 1000 == 0 {
            // if false {
                print_quat(rot_gyro_to_acc, "Rot to apply");
                println!("rot angle: {}", rot_gyro_to_acc.angle());
            }

            att_fused = rot_to_apply_to_gyro * att_gyro;

            // Careful; this replaces the gyro data.
            // self.att_from_gyros = att_fused;
            // att_gyro = att_fused;
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

            print_quat(self.attitude, "Att fused");

            print_quat(self.att_from_acc, "Att Acc");

            // print_quat(self.att_from_gyros, "Att gyros");

            println!(
                "Acc: x{} y{} z{} mag{}",
                accel_data.x,
                accel_data.y,
                accel_data.z,
                accel_data.magnitude()
            );

            // todo: Temp print to make sure gyro without heading has indeed no heading, and pitch and roll
            // todo are correct for it.
            // print_quat(att_gyro_without_heading, "Gyro without heading");

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
        // todo: Ways to identify linear acceleration:
        // - Greater or less than 1G of acceleration, if the accel is calibrated.
        // - Discontinuities or other anomolies when integrating accel-based attitude over time,
        // - or, along those lines, discontinuities etc when fusing with gyro.

        // Identify the angle difference in the vector between the current attitude estimate, and that
        // from the accelerometer alone.

        // This is the up vector as assessed from the attitude from the gyro. It is equivalent to
        // the accelerometer's normalized vector when linear acceleration is 0.
        let grav_axis_from_att_gyro = att_gyro.rotate_vec(UP);

        // An alignment of 1 indicates the estimated up direction between the accelerometer
        // and gyro match. A value of 0 means they're perpendicular.
        // let _acc_gyro_alignment = accel_data.to_normalized().dot(grav_axis_from_att_gyro);

        // Compute the difference, as a quaternion, between the attitude calulated from the accelerometer,
        // and the attitude calculated from the gyroscope. A notable difference implies linear acceleration.
        // todo: We are currently not using these
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


        // todo: Temp to tS
        // let lin_acc_estimate = accel_data.to_normalized() - grav_axis_from_att_gyro;

        self.align(lin_acc_estimate, accel_data);

        // Store our linear acc estimate and accumulator before compensating for bias.
        self.linear_acc_estimate = lin_acc_estimate; // todo: DOn't take all of it; fuse with current value.
                                                     // todo: Be careful about floating point errors over time.
                                                     // todo: Toss extreme values?
                                                     // todo: Lowpass?

        // todo: Update biases automatically for a short window after bootup.

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
            println!(
                "\n\nLin bias: x{} y{} z{}",
                self.cal.linear_acc_bias.x, self.cal.linear_acc_bias.y, self.cal.linear_acc_bias.z,
            );

            // println!(
            //     "Lin bias rem x{} y{} z{}. mag{}",
            //     lin_acc_estimate_bias_removed.x,
            //     lin_acc_estimate_bias_removed.y,
            //     lin_acc_estimate_bias_removed.z,
            //     lin_acc_estimate_bias_removed.magnitude()
            // );

            print_quat(att_diff_rot, "Att diff rot");
            println!("Att diff rot angle: {}", att_diff_rot.angle());
            let a =  grav_axis_from_att_gyro * self.cal.acc_len_at_rest;
            println!("Grav axis gyro x{} y{} z{}", a.x, a.y, a.z);
            println!("Acc x{} y{} z{}", accel_data.x, accel_data.y, accel_data.z);

            println!(
                "Lin: x{} y{} z{}. mag{}",
                lin_acc_estimate.x,
                lin_acc_estimate.y,
                lin_acc_estimate.z,
                lin_acc_estimate.magnitude() // lin_acc_estimate.x,
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

        (update_gyro_from_acc, lin_acc_estimate_bias_removed)
    }

    fn handle_mag(&mut self, mut mag_raw: Vec3, heading_gyro: f32, i: u32) {
        // // todo: Is this where we want this?
        // if self.mag_cal_in_progress {
        //     // self.cal.mag_cal_state.log_hard_iron(mag);
        //     // return;
        // } else {
        //     let mut mag = self.cal.apply_mag_cal(mag);
        // }
        //
        let mut mag = self.cal.apply_mag_cal(mag_raw);

        const EPS: f32 = 0.0000001;
        if mag.x.abs() < EPS && mag.y.abs() < EPS && mag.z.abs() < EPS {
            // todo: We get this on the first run; not sure why.
            return;
        }

        // todo: Not sure why we have to do this swap.
        // do the swap after applying cal.
        let y = mag.y;
        mag.y = mag.x;
        mag.x = y;

        let mag_norm = mag.to_normalized();

        // We use attitude from the accelerometer only here, to eliminate any
        // Z-axis component
        // todo: This is wrong. You need to find the angle only in the relevant plane.
        let mag_earth_ref = self.att_from_acc.inverse().rotate_vec(mag_norm);
        // let mag_earth_ref = self.att_from_acc.rotate_vec(mag_norm);

        let att_mag = att_from_mag(mag_norm, self.mag_inclination_estimate); // todo: Self.incl estimate
        self.att_from_mag = Some(att_mag);

        // todo: Use your fused/gyro att with heading subtracted for this, or it will be unreliable
        // todo under linear accel.
        // let heading_mag = heading_from_mag(mag_norm, att_acc);
        let heading_mag = heading_from_mag(mag_earth_ref, self.mag_declination);
        // let heading_mag = heading_from_att(att_mag);

        // Assess magnetometer health by its comparison in rate change compared to the gyro.
        match self.heading_mag {
            Some(heading_mag_prev) => {
                // todo: Find avg over many readings.
                // todo: Since even a messed up mag seems to show constant readings
                // todo when there is no rotation, consider only logging values here if
                // todo dh/dt exceeds a certain value.
                self.recent_dh_mag_dh_gyro = Some(
                    (heading_mag - heading_mag_prev) / self.dt
                        - (heading_gyro - self.heading_gyro) / self.dt,
                );

                // Fuse heading from gyro with heading from mag.
                if self.recent_dh_mag_dh_gyro.unwrap().abs() < self.config.mag_gyro_diff_thresh {
                    // heading_fused = (heading_prev * gyro_heading_weight + heading_mag * mag_heading_weight) / 2.;
                }
            }
            None => {
                self.recent_dh_mag_dh_gyro = None;
            }
        }

        // todo: For now, we are using mag only for the heading.

        self.heading_mag = Some(heading_mag);

        // todo: Inclination will be tracked and modified; not set independently every time.
        let mag_on_fwd_plane = mag_earth_ref.project_to_plane(RIGHT);

        // todo: This is currently heavily-dependent on pitch! Likely due to earth ref being wrong?
        // Negative since it's a rotation below the horizon.
        let inclination_estimate = -Quaternion::from_unit_vecs(mag_on_fwd_plane, FORWARD).angle();

        // No need to update the ratio each time.
        if i % self.config.update_ratio_mag_incl as u32 == 0 {
            if self.initialized {
                // Weighted average of current inclination estimate with stored.
                let incl_ratio = self.config.update_amt_mag_incl_estimate
                    * self.dt
                    * self.config.update_ratio_mag_incl as f32;

                self.mag_inclination_estimate = (self.mag_inclination_estimate * (1. - incl_ratio)
                    + inclination_estimate * incl_ratio);

                if i % 1000 == 0 {
                    // if false {
                    println!("Estimated mag incl: {}", inclination_estimate);
                }
            } else {
                // Take the full update on the first run.
                self.mag_inclination_estimate = inclination_estimate;
            }
        }

        // if i % 1000 == 0 {
        if false {
            println!(
                "\n\nMag norm: x{} y{} z{} len{}",
                mag_norm.x,
                mag_norm.y,
                mag_norm.z,
                mag.magnitude()
            );

            println!("Estimated mag incl Cum: {}", self.mag_inclination_estimate);

            println!(
                "Mag earth: x{} y{} z{}",
                mag_earth_ref.x, mag_earth_ref.y, mag_earth_ref.z
            );

            print_quat(att_mag, "Att mag");

            let euler_mag = att_mag.to_euler();
            println!(
                "Euler mag: p{} r{} y{}",
                euler_mag.pitch, euler_mag.roll, euler_mag.yaw,
            );

            println!("Heading mag: {}", heading_mag);
        }

        if i % self.config.update_ratio_mag_cal_log as u32 == 0 {
            self.cal.log_mag_cal_pt(self.attitude, mag_raw);

            // todo: Rework this into something more sophisticated

            let mut ready_to_cal = true;
            for cat_count in self.cal.mag_sample_count {
                if cat_count != MAG_SAMPLES_PER_CAT as u8 {
                    ready_to_cal = false;
                }
            }
            if ready_to_cal {
                self.cal.update_mag_cal();
            }
        }

        // static mut MAG_CAL_PRINTED: bool = false;
        // // todo: use your stored config value, vice this magic 50.
        // if self.mag_cal_in_progress && i % 100 == 0 {
        //     let i = MAG_CAL_I.fetch_add(1, Ordering::Relaxed);
        //
        //     if i == MAG_CAL_DATA_LEN {
        //         self.mag_cal_in_progress = false;
        //     } else {
        //         self.cal.mag_cal_data[i] = (self.attitude, mag);
        //     }
        // }
        //
        // unsafe {
        //     if !self.mag_cal_in_progress && !MAG_CAL_PRINTED {
        //         println!("\n\n[");
        //         for data in self.cal.mag_cal_data {
        //             println!("({}, {}, {}),", data.1.x, data.1.y, data.1.z);
        //         }
        //
        //         println!("\n\n]");
        //
        //         MAG_CAL_PRINTED = true;
        //     }
        // }
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
}

/// Estimate attitude from accelerometer. This will fail when under
/// linear acceleration. Apply calibration prior to this step.
/// Uses the previous attitude to rotate along the remaining degree of freedom (heading)
pub fn att_from_accel(accel_norm: Vec3) -> Quaternion {
    Quaternion::from_unit_vecs(UP, accel_norm)
}

/// Estimate attitude from magnetometer. This will fail when experiencing magnetic
/// interference, and is noisy in general. Apply calibration prior to this step.
/// Inclination is in radians.
/// todo: Automatically determine inclination by comparing att from this function
/// todo with that from our AHRS.
/// todo: Can you use this as a santity check of your ACC/Gyro heading, adn fuse it?
/// todo: The likely play is to apply mag attitude corrections when under linear acceleration
/// todo for long durations.
pub fn att_from_mag(mag_norm: Vec3, inclination: f32) -> Quaternion {
    // `mag_field_vec` points towards magnetic earth, and into the earth IOC inlination.

    let incl_rot = Quaternion::from_axis_angle(RIGHT, inclination);
    let mag_field_vec = incl_rot.rotate_vec(FORWARD);

    // println!("Mag field vec: {} {} {}", mag_field_vec.x, mag_field_vec.y, mag_field_vec.z);

    Quaternion::from_unit_vecs(mag_field_vec, mag_norm)
}

/// Calculate heading, in radians, from the magnetometer's X and Y axes.
// pub fn heading_from_mag(mag_norm: Vec3, att_without_heading: Quaternion) -> f32 {
pub fn heading_from_mag(mag_earth_ref: Vec3, declination: f32) -> f32 {
    // todo: Pass in earth ref instead of recomputing earth ref here.
    // // (mag.y.atan2(mag.x) + TAU/4.) % (TAU / 2.)
    //
    // // todo: Inv, or normal?
    // let mag_earth_ref = att_without_heading.inverse().rotate_vec(mag_norm);
    // let mag_earth_ref = att_without_heading.rotate_vec(mag_norm);

    let mag_heading = (3. * TAU / 4. - mag_earth_ref.x.atan2(mag_earth_ref.y)) % TAU;
    mag_heading + declination

    // return TAU / 4. - mag_norm.x.atan2(mag_norm.y);

    // let mag_on_horizontal_plane = mag_earth_ref.project_to_plane(UP);
    // let rot_to_fwd = Quaternion::from_unit_vecs(mag_on_horizontal_plane, FORWARD);
    // rot_to_fwd.angle()

    // From Honeywell guide
    // todo: Is this equiv to atan2?
    // if mag_earth_ref.y > 0. {
    //     TAU/4. - (mag_earth_ref.x / mag_earth_ref.y).atan()
    // } else {
    //     3. * TAU/4. - (mag_earth_ref.x / mag_earth_ref.y).atan()
    // }
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
