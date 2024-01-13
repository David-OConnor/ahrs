//! THis module is used to apply digital filters to data.

use cmsis_dsp_api as dsp_api;
use cmsis_dsp_sys as dsp_sys;

use lin_alg2::f32::Vec3;

/// Used to satisfy RTIC resource Send requirements.
pub struct IirInstWrapper {
    pub inner: dsp_sys::arm_biquad_casd_df1_inst_f32,
}

unsafe impl Send for IirInstWrapper {}

const BLOCK_SIZE: u32 = 1;

static mut FILTER_STATE_GNSS_ACCEL_X: [f32; 4] = [0.; 4];
static mut FILTER_STATE_GNSS_ACCEL_Y: [f32; 4] = [0.; 4];
static mut FILTER_STATE_GNSS_ACCEL_Z: [f32; 4] = [0.; 4];

// filter_ = signal.iirfilter(1, 300, btype="lowpass", ftype="bessel", output="sos", fs=1_000)
// coeffs = []
// for row in filter_:
//     coeffs.extend([row[0] / row[3], row[1] / row[3], row[2] / row[3], -row[4] / row[3], -row[5] / row[3]])

#[allow(clippy::excessive_precision)]
static COEFFS_LP_GNSS_ACCEL: [f32; 5] = [
    0.2452372752527856,
    0.2452372752527856,
    0.0,
    0.509525449494429,
    -0.0,
];

// todo: Calibration for IMU: Hardware, software, or both?

/// Store lowpass IIR filter instances, for use with lowpass and notch filters for IMU readings.
pub struct Filters {
    pub accel_gnss_x: IirInstWrapper,
    pub accel_gnss_y: IirInstWrapper,
    pub accel_gnss_z: IirInstWrapper,
}

impl Default for Filters {
    fn default() -> Self {
        let mut result = Self {
            accel_gnss_x: IirInstWrapper {
                inner: dsp_api::biquad_cascade_df1_init_empty_f32(),
            },
            accel_gnss_y: IirInstWrapper {
                inner: dsp_api::biquad_cascade_df1_init_empty_f32(),
            },
            accel_gnss_z: IirInstWrapper {
                inner: dsp_api::biquad_cascade_df1_init_empty_f32(),
            },
        };

        unsafe {
            // todo: Re-initialize fn?
            dsp_api::biquad_cascade_df1_init_f32(
                &mut result.accel_gnss_x.inner,
                &COEFFS_LP_GNSS_ACCEL,
                &mut FILTER_STATE_GNSS_ACCEL_X,
            );
            dsp_api::biquad_cascade_df1_init_f32(
                &mut result.accel_gnss_y.inner,
                &COEFFS_LP_GNSS_ACCEL,
                &mut FILTER_STATE_GNSS_ACCEL_Y,
            );
            dsp_api::biquad_cascade_df1_init_f32(
                &mut result.accel_gnss_z.inner,
                &COEFFS_LP_GNSS_ACCEL,
                &mut FILTER_STATE_GNSS_ACCEL_Z,
            );
        }

        result
    }
}

impl Filters {
    /// Apply the filters to IMU readings, modifying in place. Block size = 1.
    pub fn apply(&mut self, gnss_acc: &mut Vec3) {
        gnss_acc.a_x = filter_one(&mut self.accel_gnss_x, gnss_acc.a_x);
        gnss_acc.a_y = filter_one(&mut self.accel_gnss_y, gnss_acc.a_y);
        gnss_acc.a_z = filter_one(&mut self.accel_gnss_z, gnss_acc.a_z);
    }
}

/// Helper fn to reduce repetition. Applies an IIR filter to a single value.
pub fn filter_one(filter: &mut IirInstWrapper, value: f32) -> f32 {
    let mut out_buf = [0.];
    dsp_api::biquad_cascade_df1_f32(
        &mut filter.inner,
        &[value],
        &mut out_buf,
        BLOCK_SIZE,
    );
    out_buf[0]
}
