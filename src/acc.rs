//! This module contains some accelerometer-specific code

use lin_alg::f32::{Quaternion, Vec3};
use num_traits::float::Float; // abs etc

use crate::{
    attitude::{make_nudge, Ahrs, AhrsCal, NUM_LIN_ACC_CUM_SAMPLES, SAMPLES_BEFORE_ACC_CALC},
    linear_acc, print_quat, UP,
};

#[cfg(feature = "defmt")]
use defmt::println;


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
        let lin_acc_gyro = linear_acc::from_gyro(acc, *att_fused, self.cal.acc_len_at_rest);
        self.lin_acc_gyro = lin_acc_gyro;

        // let lin_acc_estimate_bias_removed = lin_acc_estimate - self.cal.linear_acc_bias;

        // todo: Rework alignment
        // self.align(lin_acc_estimate, accel_data);
        if self.cal.acc_len_count == NUM_LIN_ACC_CUM_SAMPLES + SAMPLES_BEFORE_ACC_CALC {
            self.cal.acc_len_at_rest = self.cal.acc_len_cum / NUM_LIN_ACC_CUM_SAMPLES as f32;
            // println!("\n\nAcc len at rest found: {:?}", self.cal.acc_len_at_rest);
            self.cal.acc_len_count += 1; // so this no longer triggers.
        } else if self.cal.acc_len_count > SAMPLES_BEFORE_ACC_CALC
            && self.cal.acc_len_count < NUM_LIN_ACC_CUM_SAMPLES + SAMPLES_BEFORE_ACC_CALC
        {
            self.align(acc);
        } else if self.cal.acc_len_count <= SAMPLES_BEFORE_ACC_CALC {
            self.cal.acc_len_count += 1;
        }

        // self.linear_acc_estimate = lin_acc_estimate_bias_removed;
        self.lin_acc_fused = lin_acc_gyro;

        let lin_acc_estimate_bias_removed = lin_acc_gyro; // todo: For now we removed bias removal

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

        let att_acc_w_lin_removed = att_from_accel((acc - lin_acc_gyro).to_normalized());

        // Make sure we update heading_gyro after mag handling; we use it to diff gyro heading.
        // todo: Remove `heading_gyro` if you end up not using it.
        // self.heading_gyro = heading_fused;

        // todo: Instead of a binary update-or-not, consider weighing the slerp value based
        // todo on how much lin acc we assess, or how much uncertainly in lin acc.

        if update_gyro_from_acc {
            // Apply a rotation of the gyro solution towards the acc solution, if we think we are not under
            // much linear acceleration.
            let rot_correction = make_nudge(
                self.attitude,
                accel_norm,
                UP,
                self.config.update_amt_att_from_mag * self.dt,
            );

            *att_fused = rot_correction * *att_fused;
        }

        #[cfg(feature = "defmt")]
        // if self.num_updates % ((1. / self.dt) as u32) == 0 {
        if false {
            println!("Acc cal x{} y{} z{}", acc.x, acc.y, acc.z);

            println!(
                "\nLin acc gyro: x{} y{} z{} mag{}",
                lin_acc_gyro.x,
                lin_acc_gyro.y,
                lin_acc_gyro.z,
                lin_acc_gyro.magnitude(),
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
