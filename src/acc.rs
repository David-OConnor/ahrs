//! This module contains some accelerometer-specific code

use num_traits::float::Float; // abs etc

use lin_alg2::f32::{Quaternion, Vec3};

use defmt::println;

use crate::{
    attitude::{Ahrs, AhrsCal},
    linear_acc, print_quat, UP,
};

impl AhrsCal {
    fn apply_cal_acc(&self, data: Vec3) -> Vec3 {
        Vec3::new(
            data.x * self.acc_slope.x - self.acc_bias.x,
            data.y * self.acc_slope.y - self.acc_bias.y,
            data.z * self.acc_slope.z - self.acc_bias.z,
        )
    }
}

impl Ahrs {
    pub(crate) fn handle_acc(&mut self, acc_raw: Vec3, att_fused: &mut Quaternion) {
        let acc = self.cal.apply_cal_acc(acc_raw);
        self.acc_calibrated = acc;

        let accel_norm = acc.to_normalized();

        // Estimate attitude from raw accelerometer and gyro data. Note that
        // The gyro data reguarly receives updates from the acc and mag.
        let att_acc = att_from_accel(accel_norm);
        self.att_from_acc = att_acc;

        // See comment on the `initialized` field.
        // We update initialized state at the end of this function, since other steps rely on it.
        if !self.initialized {
            *att_fused = att_acc;
        }

        // todo: YOu may wish to apply a lowpass filter to linear acc estimate.
        let lin_acc_estimate = linear_acc::from_gyro(acc, *att_fused, self.cal.acc_len_at_rest);

        if let Some(lin_acc_gnss) = self.lin_acc_gnss {
            if let Some(fix) = &self.fix_prev {
                // todo: Put back once you figure out how to compare current time to this.
                if self.timestamp - fix.timestamp_s > self.config.max_fix_age_lin_acc {
                    self.lin_acc_gnss = None;
                } else if self.num_updates % ((1. / self.dt) as u32) == 0 {
                    println!(
                        "Lin acc GNSS: x{} y{} z{} mag{}",
                        lin_acc_gnss.x,
                        lin_acc_gnss.y,
                        lin_acc_gnss.z,
                        lin_acc_gnss.magnitude()
                    );
                }
                // todo: Here etc, include your fusing with gyro lin acc estimate.
            }
        }

        // let lin_acc_estimate_bias_removed = lin_acc_estimate - self.cal.linear_acc_bias;

        // todo: Rework alignment
        // self.align(lin_acc_estimate, accel_data);

        // self.linear_acc_estimate = lin_acc_estimate_bias_removed;
        self.linear_acc_estimate = lin_acc_estimate;

        let lin_acc_estimate_bias_removed = lin_acc_estimate; // todo: For now we removed bias removal

        // todo: Move this update_gyro_from_acc logic elsewhere, like a dedicated fn; or, rework it.
        let mut update_gyro_from_acc = false;
        // If it appears there is negligible linear acceleration, update our gyro readings as appropriate.
        if (acc.magnitude() - self.cal.acc_len_at_rest).abs() < self.config.total_accel_thresh {
            // We guess no linear acc since we're getting close to 1G. Note that
            // this will produce false positives in some cases.
            update_gyro_from_acc = true;
        }
        // else if lin_acc_estimate_bias_removed.magnitude() < self.config.lin_acc_thresh {
        //     // If not under much acceleration, re-cage our attitude.
        //     update_gyro_from_acc = true;
        // }

        let att_acc_w_lin_removed = att_from_accel((acc - lin_acc_estimate).to_normalized());

        // Make sure we update heading_gyro after mag handling; we use it to diff gyro heading.
        // todo: Remove `heading_gyro` if you end up not using it.
        // self.heading_gyro = heading_fused;

        // todo: Instead of a binary update-or-not, consider weighing the slerp value based
        // todo on how much lin acc we assess, or how much uncertainly in lin acc.

        if update_gyro_from_acc {
            // Apply a rotation of the gyro solution towards the acc solution, if we think we are not under
            // much linear acceleration.
            // This rotation is heading-invariant: Rotate the gyro *up* towards the acc *up*.
            let gyro_up = att_fused.rotate_vec(UP);
            let rot_gyro_to_acc = Quaternion::from_unit_vecs(gyro_up, accel_norm);

            let rot_acc_correction = Quaternion::new_identity().slerp(
                rot_gyro_to_acc,
                self.config.update_amt_att_from_acc * self.dt,
            );

            // if self.num_updates % ((1. / self.dt) as u32) == 0 {
            if false {
                print_quat(rot_gyro_to_acc, "Rot to apply");
                println!("rot angle: {}", rot_gyro_to_acc.angle());
            }

            *att_fused = rot_acc_correction * *att_fused;
        }

        if self.num_updates % ((1. / self.dt) as u32) == 0 {
            println!("Acc cal x{} y{} z{}", acc.x, acc.y, acc.z);

            println!(
                "\nLin acc: x{} y{} z{} mag{}",
                lin_acc_estimate.x,
                lin_acc_estimate.y,
                lin_acc_estimate.z,
                lin_acc_estimate.magnitude(),
            );
        }
    }
}

/// Estimate attitude from accelerometer. This will fail when under
/// linear acceleration. Apply calibration prior to this step.
/// Uses the previous attitude to rotate along the remaining degree of freedom (heading)
pub fn att_from_accel(accel_norm: Vec3) -> Quaternion {
    Quaternion::from_unit_vecs(UP, accel_norm)
}
