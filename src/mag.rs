//! Code releated to magnetometer reading and calibration.

use core::f32::consts::TAU;

use num_traits::float::Float; // abs etc

// todo: Assess mag health using consistency of mag vector len (calibrated). It should be
// todo close to 1. THrow out readings that aren't close to 1.

use lin_alg2::f32::{Quaternion, Vec3};

use defmt::println;

use crate::{
    attitude::{make_nudge, Ahrs, AhrsCal},
    blend,
    mag_ellipsoid_fitting::{
        self, MAG_SAMPLES_PER_CAT, SAMPLE_VERTEX_ANGLE, SAMPLE_VERTICES, TOTAL_MAG_SAMPLE_PTS,
    },
    ppks::PositVelEarthUnits,
    print_quat,
    util::map_linear,
    FORWARD, RIGHT, UP,
};

impl AhrsCal {
    /// Apply the hard and soft iron offsets to our readings.
    pub fn apply_cal_mag(&self, data: Vec3) -> Vec3 {
        // todo: Why this clone?
        // todo: Temp removing soft iron to TS
        self.soft_iron.clone() * (data - self.hard_iron)
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

        // Indexed by the sample vec vertices. We use 2 perpendicular vector
        // to remove a rotational ambiguity using a single would do.

        for (data, direction, i, count) in [
            (
                &mut self.mag_cal_data_up,
                UP,
                &mut self.mag_sample_i_up,
                &mut self.mag_sample_count_up,
            ),
            (
                &mut self.mag_cal_data_fwd,
                FORWARD,
                &mut self.mag_sample_i_fwd,
                &mut self.mag_sample_count_fwd,
            ),
        ] {
            let mut sample_category = 0;
            for (cat, sample) in SAMPLE_VERTICES.iter().enumerate() {
                // `UP` is arbitrary; any Vec will do as long as we're consistent.

                let rotated = att.rotate_vec(direction);
                if (rotated.dot(*sample)).acos() < SAMPLE_VERTEX_ANGLE {
                    sample_category = cat;
                    break;
                }
            }
            data[sample_category][i[sample_category]] = mag_raw;

            // This loops around; we always update the oldest value in each category
            i[sample_category] = (i[sample_category] + 1) % MAG_SAMPLES_PER_CAT;

            if count[sample_category] < MAG_SAMPLES_PER_CAT as u8 {
                count[sample_category] += 1;
            }
        }
    }

    /// Update mag calibration based on sample points.
    /// Only run this when there is recent data of each attitude category of points.
    pub fn update_mag_cal(&mut self, update_amt: f32) {
        // Flatten our sample data organized by points; that format is only used to ensure a
        // good selection of points, ie fairly weighted on all sides.
        let mut sample_pts = [Vec3::new_zero(); TOTAL_MAG_SAMPLE_PTS];

        const EPS: f32 = 0.0000001;

        println!("\nMag cal complete\n");

        let mut i = 0;
        // println!("\n\n\n[");
        for data in &[self.mag_cal_data_up, self.mag_cal_data_fwd] {
            for cat in data {
                for point in cat {
                    if point.x.abs() < EPS && point.y.abs() < EPS && point.z.abs() < EPS {
                        continue;
                    }

                    sample_pts[i] = *point;
                    // println!("({}, {}, {}),", point.x, point.y, point.z);
                    i += 1;
                }
            }
        }
        // println!("]");

        // todo: Put back once working.

        let poly_terms = mag_ellipsoid_fitting::ls_ellipsoid(&sample_pts);
        let (hard_iron, soft_iron) = mag_ellipsoid_fitting::poly_to_params_3d(&poly_terms);
        // // todo: axes and inve; what are they? Check the web page.
        // self.hard_iron = center;
        // self.soft_iron = inve;

        // println!(
        //     "Hard iron: x{} y{} z{}",
        //     hard_iron.x, hard_iron.y, hard_iron.z
        // );

        // println!(
        //     "Soft iron diag: {} {} {} {} {} {}",
        //     soft_iron.data[0],
        //     soft_iron.data[1],
        //     soft_iron.data[2],
        //     soft_iron.data[4],
        //     soft_iron.data[5],
        //     soft_iron.data[8]
        // );

        let update_amt_inv = 1. - update_amt;

        self.hard_iron = self.hard_iron * update_amt_inv + hard_iron * update_amt;
        self.soft_iron = self.soft_iron.clone() * update_amt_inv + soft_iron * update_amt;

        // Reset our sample counters. We leave the sample data in place, since we can still
        // use the readings in the next calibration. (Since we initiate cal without completely
        // filling it)
        self.mag_sample_i_up = Default::default();
        self.mag_sample_i_fwd = Default::default();
        self.mag_sample_count_up = Default::default();
        self.mag_sample_count_fwd = Default::default();

        self.mag_cal_updated = true; // set to false in firmware once written to flash.
    }
}

impl Ahrs {
    /// We use magnetometer information in two ways. Note that its ambiguity axis is close to
    /// the ambitguity access of the accelerometer, but it provides some heading information, while
    /// the acc provides none. We apply corrections in two ways:
    ///
    /// 1: Create an attitude using the assessed inclination vector, and the magnetometer reading.
    /// 2: Extract heading directly, and apply a correction: This is our primary absolute heading reference.
    pub(crate) fn handle_mag(&mut self, mag_raw: Vec3, att_fused: &mut Quaternion) {
        let mag = self.cal.apply_cal_mag(mag_raw);
        self.mag_calibrated = Some(mag);

        const EPS: f32 = 0.0000001;
        if mag.x.abs() < EPS && mag.y.abs() < EPS && mag.z.abs() < EPS {
            // todo: We get this on the first run; not sure why.
            return;
        }

        let mag_norm = mag.to_normalized();

        let incl_rot = Quaternion::from_axis_angle(RIGHT, -self.mag_inclination_estimate);
        let mag_field_absolute = incl_rot.rotate_vec(FORWARD);

        let att_mag = att_from_mag(mag_norm, mag_field_absolute);

        // todo
        // let att_mag = apply_declination(att_fused, declination);

        self.att_from_mag = Some(att_mag);

        let magnetometer_magnitude = mag.magnitude(); // Say that 10 times fast?

        let update_amt_mag_var = 0.10 * self.dt; // todo: Store this const as a struct param.
                                                 // let mag_variance = (mag.magnitude() - 1.).powi(2);
        let mag_variance = (mag.magnitude() - 1.).abs();
        self.recent_mag_variance =
            blend(self.recent_mag_variance, mag_variance, update_amt_mag_var);

        // Symmetry with acc here; similar logic etc.
        // todo: Move this update_gyro_from_acc logic elsewhere, like a dedicated fn; or, rework it.
        let mut update_gyro_from_mag = false;
        // If it appears there is negligible non-calibrated-away interference, update our
        // gyro readings as appropriate.
        if (magnetometer_magnitude - 1.).abs() < self.config.total_mag_thresh
            && self.mag_inclination_estimate < self.config.mag_incl_max
        {
            update_gyro_from_mag = true;
        }

        if update_gyro_from_mag {
            let rot_correction_from_att = make_nudge(
                self.attitude,
                mag_norm,
                mag_field_absolute,
                self.config.update_amt_att_from_mag * self.dt,
            );


            // We extract the pure heading from the magnetometer, since that is what we primarily don't get
            // from the accelerometer, and we only get a small part from the magnetometer directly.

            // todo: We don't want just acc here, due to linear-acc problems, but do it for now to test.
            let tilt_compensated = self.att_from_acc.inverse().rotate_vec(mag_norm);
            let hdg = tilt_compensated.x.atan2(tilt_compensated.y);

            let heading_rotation = Quaternion::from_unit_vecs(
                self.attitude.rotate_vec(FORWARD).project_to_plane(UP),
                Vec3::new(hdg.sin(), hdg.cos(), 0.),
            );

            // todo: make sure you apply an instantaneous correction on init.
            let rot_correction_from_heading = Quaternion::new_identity().slerp(heading_rotation, 
                self.config.update_amt_hdg_from_mag * self.dt);

            *att_fused = rot_correction_from_att * rot_correction_from_heading * *att_fused;

            if self.num_updates % ((1. / self.dt) as u32) == 0 {
            // if false {
                println!("HDG: {:?}", hdg);
                print_quat(heading_rotation, "HDG ROT");
                println!(" FWD mag: x{} y{}", hdg.sin(), hdg.cos());

                let a = att_fused.rotate_vec(FORWARD).project_to_plane(UP);
                println!(" FWD att: x{} y{} z: {}", a.x, a.y, a.z);
                // print_quat(rot_correction_from_heading * self.att_from_mag.unwrap(), "Att mag w hdg");
            }
        }

        self.update_mag_incl(mag_norm);

        if self.num_updates % ((1. / self.dt) as u32) == 0 {
            // if false {
            println!(
                "\n\nMag raw: x{} y{} z{} len{}",
                mag_raw.x,
                mag_raw.y,
                mag_raw.z,
                mag_raw.magnitude()
            );

            println!(
                "Mag: x{} y{} z{} len{}",
                mag.x,
                mag.y,
                mag.z,
                mag.magnitude()
            );

            println!("Estimated mag incl: {}", self.mag_inclination_estimate);
            println!("Recent mag var: {}", self.recent_mag_variance);

            print_quat(att_mag, "Att mag");
        }

        if self.num_updates % self.config.update_ratio_mag_cal_log as u32 == 0 {
            self.update_mag_cal(mag_raw);
        }
    }

    /// Estimate magnetic inclination, by calculating the angle between Up in the aircraft, and the magnetometer
    /// vector. The resulting value is in radians, down from the horizon.
    fn update_mag_incl(&mut self, mag_norm: Vec3) {
        const TAU_DIV_4: f32 = TAU / 4.;

        let up = self.att_from_acc.rotate_vec(UP);

        // Angle between up and the mag reading. We subtract tau/4 to get angle below the horizon.
        let inclination_estimate = up.dot(mag_norm).acos() - TAU_DIV_4;

        // No need to update the ratio each time.
        if self.num_updates % self.config.update_ratio_mag_incl as u32 == 0 {
            if self.initialized {
                // Weighted average of current inclination estimate with stored.
                let incl_ratio = self.config.update_amt_mag_incl_estimate
                    * self.dt
                    * self.config.update_ratio_mag_incl as f32;

                self.mag_inclination_estimate = blend(
                    self.mag_inclination_estimate,
                    inclination_estimate,
                    incl_ratio,
                );
            } else {
                // Take the full update on the first run.
                // todo: That doesn't appear to work; showing 0 start. use a sane default.
                self.mag_inclination_estimate = 1.2;
                // self.mag_inclination_estimate = inclination_estimate;
            }
        }
    }

    /// Log a magnetic calibration point. Check if we've logged enough points to apply the calibration,
    /// and do if we do.
    fn update_mag_cal(&mut self, mag_raw: Vec3) {
        self.cal.log_mag_cal_pt(self.attitude, mag_raw);

        // If we've filled up most of our sample slots, divided by directional category, initiate cal.
        // Note that additional samples in a filled cat don't count towards this value.
        let mut samples_taken = 0;

        for cat_count in self.cal.mag_sample_count_up {
            samples_taken += cat_count;
        }
        for cat_count in self.cal.mag_sample_count_fwd {
            samples_taken += cat_count;
        }

        if samples_taken as f32 / TOTAL_MAG_SAMPLE_PTS as f32 >= self.config.mag_cal_portion_req
        {
            // Ie, if the recent magnetometer magnitude is too high, take all or almost all of the update.
            // todo: Store these consts in the main struct?
            let mut update_amt_mag_cal =
                map_linear(self.recent_mag_variance, (0., 0.3), (0.1, 1.));
            if update_amt_mag_cal > 1. {
                update_amt_mag_cal = 1.;
            }
            self.cal.update_mag_cal(update_amt_mag_cal);
        }
    }
}

/// Estimate attitude from magnetometer. This will fail when experiencing magnetic
/// interference, and is noisy in general. Apply calibration prior to this step.
/// The magnetic field vector points points towards magnetic earth, and into the earth IOC inlination.
/// It is points forward and down in our coordinate system.
pub fn att_from_mag(mag_norm: Vec3, mag_field_vec_absolute: Vec3) -> Quaternion {
    Quaternion::from_unit_vecs(mag_field_vec_absolute, mag_norm)
}

/// todo: Put this in a separate module A/R
/// Estimate magnetic inclination, by looking up from a table based on geographic position.
fn declination_from_posit(lat_e8: i64, lon_e8: i64) -> f32 {
    0.
}

// todo: Make this work, then use it.
fn apply_declination(att: Quaternion, declination: f32) -> Quaternion {
    // todo: QC order; trial_error
    let up_rel_earth = att.rotate_vec(UP); // todo: QC this!
    // todo: QC direction.
    Quaternion::from_axis_angle(up_rel_earth, declination).to_normalized() * att
    // Quaternion::from_axis_angle(UP, declination).to_normalized() * att
}