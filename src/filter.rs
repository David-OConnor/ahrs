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
        unsafe {
            Self {
                accel_gnss_x: IirInstWrapper { inner: iir_new(&COEFFS_LP_GNSS_ACCEL, &mut FILTER_STATE_GNSS_ACCEL_X, ) },
                accel_gnss_y: IirInstWrapper { inner: iir_new(&COEFFS_LP_GNSS_ACCEL, &mut FILTER_STATE_GNSS_ACCEL_Y, ) },
                accel_gnss_z: IirInstWrapper { inner: iir_new(&COEFFS_LP_GNSS_ACCEL, &mut FILTER_STATE_GNSS_ACCEL_Z, ) },
            }
        }
    }
}

impl Filters {
    /// Apply the filters to IMU readings, modifying in place. Block size = 1.
    pub fn apply(&mut self, gnss_acc: &mut Vec3) {
        gnss_acc.a_x = iir_apply(&mut self.accel_gnss_x, gnss_acc.a_x);
        gnss_acc.a_y = iir_apply(&mut self.accel_gnss_y, gnss_acc.a_y);
        gnss_acc.a_z = iir_apply(&mut self.accel_gnss_z, gnss_acc.a_z);
    }
}

/// Helper fn to reduce repetition. Applies an IIR filter to a single value.
pub fn iir_apply(filter: &mut IirInstWrapper, value: f32) -> f32 {
    dsp_api::iir_apply(&mut filter.inner, value)
}
