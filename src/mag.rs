//! Code releated to magnetometer reading and calibration.

use core::f32::consts::TAU;

use num_traits::float::Float; // abs etc

use crate::{
    attitude::{Ahrs, AhrsCal},
    mag_ellipsoid_fitting::{
        self, MAG_SAMPLES_PER_CAT, SAMPLE_VERTEX_ANGLE, SAMPLE_VERTICES, TOTAL_MAG_SAMPLE_PTS,
    },
    print_quat, FORWARD, RIGHT, UP,
};

// todo: Assess mag health using consistency of mag vector len (calibrated). It should be
// todo close to 1. THrow out readings that aren't close to 1.

use lin_alg2::f32::{Quaternion, Vec3};

use crate::ppks::PositVelEarthUnits;
use defmt::println;

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

        // To display status.
        if false {
            // todo: You should multiply this by the portion req.
            let mut num_pts_left = TOTAL_MAG_SAMPLE_PTS;
            for cat_count in &self.mag_sample_count_up {
                num_pts_left -= *cat_count as usize;
            }
            for cat_count in &self.mag_sample_count_fwd {
                num_pts_left -= *cat_count as usize;
            }

            println!("Logging mag pt. Num left: {}", num_pts_left);
        }
    }

    /// Update mag calibration based on sample points.
    /// Only run this when there is recent data of each attitude category of points.
    pub fn update_mag_cal(&mut self, update_amt: f32) {
        // Flatten our sample data organized by points; that format is only used to ensure a
        // good selection of points, ie fairly weighted on all sides.
        let mut sample_pts = [Vec3::new_zero(); TOTAL_MAG_SAMPLE_PTS];

        const EPS: f32 = 0.0000001;

        let mut i = 0;
        println!("\n\n\n[");
        for data in &[self.mag_cal_data_up, self.mag_cal_data_fwd] {
            for cat in data {
                for point in cat {
                    if point.x.abs() < EPS && point.y.abs() < EPS && point.z.abs() < EPS {
                        continue;
                    }

                    sample_pts[i] = *point;
                    println!("({}, {}, {}),", point.x, point.y, point.z);
                    i += 1;
                }
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

        let update_amt_inv = 1. - update_amt;

        self.hard_iron = self.hard_iron.clone() * update_amt_inv + hard_iron * update_amt;
        self.soft_iron = self.soft_iron.clone() * update_amt_inv + soft_iron * update_amt;

        // todo: Save to flash, but not that this should probably not be done while airborne
        // todo as flash access can cause the program to hang. (?)

        // Reset our sample counters. We leave the sample data in place, since we can still
        // use the readings in the next calibration. (Since we initiate cal without completely
        // filling it)
        // self.mag_cal_data_up = Default::default();
        // self.mag_cal_data_fwd = Default::default();
        self.mag_sample_i_up = Default::default();
        self.mag_sample_i_fwd = Default::default();
        self.mag_sample_count_up = Default::default();
        self.mag_sample_count_fwd = Default::default();
    }
}

impl Ahrs {
    pub(crate) fn handle_mag(
        &mut self,
        mag_raw: Vec3,
        att_fused: &mut Quaternion,
        heading_fused: f32,
    ) {
        let mag = self.cal.apply_cal_mag(mag_raw);

        const EPS: f32 = 0.0000001;
        if mag.x.abs() < EPS && mag.y.abs() < EPS && mag.z.abs() < EPS {
            // todo: We get this on the first run; not sure why.
            return;
        }

        // todo: Not sure why we have to do this swap.
        // do the swap after applying cal.
        // let y = mag.y;
        // mag.y = mag.x;
        // mag.x = y;

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

        let heading_mag = mag_norm.x.atan2(mag_norm.y); // todo temp

        let heading_diff = -heading_mag - heading_fused;

        let rotation = if self.initialized {
            // Accept the whole magnetometer update.
            Quaternion::from_axis_angle(UP, heading_diff)
        } else {
            // Nudge the heading in the direction of the gyro.
            Quaternion::from_axis_angle(
                UP,
                heading_diff * self.config.update_amt_mag_heading * self.dt,
            )
        };

        // todo: Put back.
        // *att_fused = rotation * *att_fused;

        // Assess magnetometer health by its comparison in rate change compared to the gyro.
        match self.heading_mag {
            Some(heading_mag_prev) => {
                // todo: Find avg over many readings.
                // todo: Since even a messed up mag seems to show constant readings
                // todo when there is no rotation, consider only logging values here if
                // todo dh/dt exceeds a certain value.

                let dmag_dt = (heading_mag - heading_mag_prev) / self.dt;
                let dgyro_dt = (heading_fused - self.heading_gyro) / self.dt;

                self.recent_dh_mag_dh_gyro = Some(dmag_dt - dgyro_dt);

                if self.num_updates % ((1. / self.dt) as u32) == 0 {
                    println!(
                        "Dmag: {} Dgyro: {} diff: {}",
                        dmag_dt, dgyro_dt, self.recent_dh_mag_dh_gyro
                    );
                }

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

        self.estimate_mag_inclination(mag_earth_ref);

        if self.num_updates % ((1. / self.dt) as u32) == 0 {
            // if false {
            //     println!(
            //         "\n\nMag: x{} y{} z{} len{}",
            //         mag_norm.x,
            //         mag_norm.y,
            //         mag_norm.z,
            //         mag.magnitude()
            //     );

            // let xy_norm = (mag.x.powi(2) + mag.y.powi(2)).sqrt();
            // println!("\n\nMag xy: x{} y{}", mag.x / xy_norm, mag.y / xy_norm,);

            println!(
                "\n\nMag raw: x{} y{} z{} len{}",
                mag_raw.x,
                mag_raw.y,
                mag_raw.z,
                mag_raw.magnitude()
            );

            // println!(
            //     "\n\nMag raw norm: x{} y{} z{}",
            //     mag_raw.to_normalized().x,
            //     mag_raw.to_normalized().y,
            //     mag_raw.to_normalized().z,
            // );

            println!(
                "Mag: x{} y{} z{} len{}",
                mag.x,
                mag.y,
                mag.z,
                mag.magnitude()
            );

            // println!("Estimated mag incl Cum: {}", self.mag_inclination_estimate);

            // println!(
            //     "Mag earth: x{} y{} z{}",
            //     mag_earth_ref.x, mag_earth_ref.y, mag_earth_ref.z
            // );

            print_quat(att_mag, "Att mag");

            let euler_mag = att_mag.to_euler();
            println!(
                "Euler mag: p{} r{} y{}",
                euler_mag.pitch, euler_mag.roll, euler_mag.yaw,
            );

            println!("Heading mag: {}", heading_mag);
        }

        if self.num_updates % self.config.update_ratio_mag_cal_log as u32 == 0 {
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
                self.cal.update_mag_cal(self.config.mag_cal_update_amt);
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

    /// Estimate magnetic inclination. The returned value is in readians, down from
    /// the horizon.
    fn estimate_mag_inclination(&mut self, mag_earth_ref: Vec3) {
        // todo: Inclination will be tracked and modified; not set independently every time.
        let mag_on_fwd_plane = mag_earth_ref.project_to_plane(RIGHT);

        // todo: This is currently heavily-dependent on pitch! Likely due to earth ref being wrong?
        // Negative since it's a rotation below the horizon.
        let inclination_estimate = -Quaternion::from_unit_vecs(mag_on_fwd_plane, FORWARD).angle();

        // No need to update the ratio each time.
        if self.num_updates % self.config.update_ratio_mag_incl as u32 == 0 {
            if self.initialized {
                // Weighted average of current inclination estimate with stored.
                let incl_ratio = self.config.update_amt_mag_incl_estimate
                    * self.dt
                    * self.config.update_ratio_mag_incl as f32;

                self.mag_inclination_estimate = (self.mag_inclination_estimate * (1. - incl_ratio)
                    + inclination_estimate * incl_ratio)
            } else {
                // Take the full update on the first run.
                self.mag_inclination_estimate = inclination_estimate;
            }
        }
    }
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

    // let mag_heading = (3. * TAU / 4. - mag_earth_ref.x.atan2(mag_earth_ref.y)) % TAU;
    // let mag_heading = (1. * TAU / 4. - mag_earth_ref.x.atan2(mag_earth_ref.y)) % TAU;
    let mag_heading = (1. * TAU / 4. - mag_earth_ref.x.atan2(mag_earth_ref.y)) % TAU;
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

/// todo: Put this in a separate module A/R
/// Estimate magnetic inclination, by looking up from a table based on geographic position.
fn declination_from_posit(lat_e8: i64, lon_e8: i64) -> f32 {
    0.
}
