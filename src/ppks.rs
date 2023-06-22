//! Present-position keeping system. Fuses GPS, and dead-reckoning.
//!
//! Sync Corvus' implemention from this.

// todo: Move to ahrs module, but we need to sort out how to handle `Fix`.

use num_traits::Float;

use lin_alg2::f32::{Quaternion, Vec3};

use crate::{params::Params, Fix};

const FIX_FUSED_SIZE: usize = 36;

const DEG_SCALE_1E8: f32 = 100_000_000.;

use defmt::println;

/// Convert NED velocity in mm/s to xyz velocity in m/s. (Still in earth frame.
fn ned_vel_to_xyz(ned_vel: [i32; 3]) -> Vec3 {
    Vec3::new(
        ned_vel[1] as f32 / 1_000.,
        ned_vel[0] as f32 / 1_000.,
        -ned_vel[2] as f32 / 1_000.,
    )
}

/// This includes position information based on fusing GNSS and INS data.
/// It doesn't include much of the metadata of a GNSS fix. We use Ublox
/// conventions when able. For example, lat and lon are in deg / 1e7.
/// todo: You may need to copy the exact bi lens of Fix2, and use PackedStruct.
#[derive(Default, Clone)]
pub struct PositFused {
    pub timestamp_us: u64, // microseconds
    pub lat_e8: i64,
    pub lon_e8: i64,
    pub elevation_hae: f32,     // meters
    pub elevation_msl: f32,     // meters
    pub ned_velocity: [f32; 3], // m/s
}

impl PositFused {
    /// `fix` is the most-recently-received fix. `attitude` and `position_inertial` are the most recent ones
    /// as well.
    ///
    /// todo: We need to update position_inertial when we get a new fix, taking some or all
    /// todo of the fix.
    pub fn new(
        fix: &Fix,
        inertial: &mut PositInertial,
        params: &Params,
        dt: f32,
        timestamp: f32,
    ) -> Self {
        let gnss_dr = PositVelEarthUnits::from_fix_dr(fix, timestamp);

        inertial.update(params, dt);

        // We use separate (vice a single float) weight since we also need to apply it to integers.
        let weight_gnss = 3;
        let weight_imu = 10 - weight_gnss;
        let weight_gnss_float = weight_gnss as f32 / 10.;
        let weight_imu_float = 1. - weight_gnss_float;

        // todo: If dt since last fix is longer, weigh INS more.

        let inertial_earth = inertial.get_earth_units();

        let ned_velocity = [
            (fix.ned_velocity[0] as f32 / 1_000.) + inertial.v_y,
            (fix.ned_velocity[1] as f32 / 1_000.) + inertial.v_x,
            (fix.ned_velocity[2] as f32 / 1_000.) + inertial.v_z,
        ];

        let elevation_msl = gnss_dr.elevation_msl * weight_gnss_float
            + inertial_earth.elevation_msl * weight_imu_float;

        // todo: Qc
        let elevation_hae = (fix.elevation_hae as f32 / 1_000.)
            - (fix.elevation_msl - fix.elevation_hae) as f32 / 1_000.;

        Self {
            // todo: Inertial is in mm I think. Maybe have it return lat/lon instead.
            timestamp_us: (timestamp * 1_000_000.) as u64,
            lon_e8: (gnss_dr.lon_e8 * weight_gnss + inertial_earth.lon_e8 * weight_imu) / 10,
            lat_e8: (gnss_dr.lat_e8 * weight_gnss + inertial_earth.lat_e8 * weight_imu) / 10,
            elevation_msl,
            elevation_hae,
            ned_velocity,
        }
    }

    /// We use this as our CAN wire format for fused fixes.
    pub fn to_bytes(&self) -> [u8; FIX_FUSED_SIZE] {
        let mut result = [0; FIX_FUSED_SIZE];

        result[0..8].copy_from_slice(&self.timestamp_us.to_le_bytes());
        result[8..12].copy_from_slice(&self.lat_e8.to_le_bytes());
        result[12..16].copy_from_slice(&self.lon_e8.to_le_bytes());
        result[16..20].copy_from_slice(&self.elevation_hae.to_le_bytes());
        result[20..24].copy_from_slice(&self.elevation_msl.to_le_bytes());

        result[24..28].copy_from_slice(&self.ned_velocity[0].to_le_bytes());
        result[28..32].copy_from_slice(&self.ned_velocity[1].to_le_bytes());
        result[32..36].copy_from_slice(&self.ned_velocity[2].to_le_bytes());

        result
    }
}

/// A simple intertial position. Units are in m and m/s. The `posit` and `velocity` values are in
/// relation to the anchor. They are East, North, up.
#[derive(Default)]
pub struct PositInertial {
    /// This anchor is what we add the other readings to. It's not necessarily dead-recocking,
    /// but the info we need is the same as that struct.
    pub anchor: PositVelEarthUnits,
    // pub anchor_heading: f32,
    /// Position in meters, earth frame. East, North, Up.
    pub posit: Vec3,
    /// Velocity in m/s, earth frame. East, North, Up.
    pub velocity: Vec3,
}

impl PositInertial {
    /// Update position intertially. Attitude is the aircraft's orientation
    /// relative to the earth's surface. Note that we don't use the gyro readings
    /// directly here. (We use accel readings). We use gyro readings as part of
    /// attitude determination.
    ///
    /// We fuse with GNSS and baro elsewhere.
    pub fn update(
        &mut self,
        params: &Params,
        dt: f32, // seconds
    ) {
        // todo: temp hard Hard-codefor testing.
        self.anchor = PositVelEarthUnits {
            lat_e8: 3500000000,
            lon_e8: -7800000000,
            elevation_msl: 0.,
            velocity: Vec3::new(0., 1., 0.),
        };

        // todo: Using this to experiment with new att platform
        let accel_data = Vec3 {
            x: params.a_x,
            y: params.a_y, // negative due to our IMU's coord system.
            z: params.a_z,
        };

        let gyro_data = Vec3 {
            x: params.v_pitch,
            y: params.v_roll,
            z: params.v_yaw,
        };
        let mag_data = Vec3 {
            x: params.mag_x, // negative due to our mag's coord systeparams.mag_.
            y: params.mag_y,
            z: params.mag_z,
        };

        // println!("Accel mag: {}", accel_vec.magnitude());
        // println!("Magnetic mag: {}", mag_vec.magnitude());

        // let accel_mag = 10.084;
        // let mag_mag = 1.0722; // gauss
        // let local_mag_str = 0.47;

        // accel_data *= crate::G / accel_mag;
        // mag_data *= local_mag_str / mag_mag;

        // let accel_lin = crate::attitude::get_linear_accel(accel_data, att_accel);
        let accel_lin = params.accel_linear;

        self.posit += self.velocity * dt;
        self.velocity += accel_lin * dt;

        static mut I: u32 = 0;
        unsafe { I += 1 };

        if unsafe { I } % 2000 == 0 {
            println!(
                "Inertial: x{} y{} z{} -- vx{} vy{} vz{}",
                self.x, self.y, self.z, self.v_x, self.v_y, self.v_z
            );
        }
    }

    /// Update the anchor point, based on a new fix.
    pub fn update_anchor(&mut self, fix: &Fix) {
        // todo: For now, we accept the whole fix. In the future, weight with our inertial position.
        let diff = ned_vel_to_xyz(fix.ned_velocity) - self.velocity;

        // let update_factor = 0.8; // 1.0 means replace the anchor entirely.
        let update_factor = 1.; // 1.0 means replace the anchor entirely.
        let velocity = self.velocity + diff * update_factor;

        self.anchor = PositVelEarthUnits {
            lat_e8: fix.lat as i64 * 10,
            lon_e8: fix.lon as i64 * 10,
            elevation_msl: fix.elevation_msl as f32 / 1_000.,
            velocity,
        };
    }

    pub fn get_earth_units(&self) -> PositVelEarthUnits {
        let (lat, lon) = augment_latlon(self.anchor.lat_e8, self.anchor.lon_e8, self.y, self.x);

        let velocity = self.anchor.velocity + self.velocity;

        PositVelEarthUnits {
            lat_e8: lat,
            lon_e8: lon,
            elevation_msl: self.anchor.elevation_msl / 1_000. + self.z,
            velocity,
        }
    }
}

/// Lat and lon are x 1e8, and alt in m. We use this for several things.
#[derive(Clone, Default)]
pub struct PositVelEarthUnits {
    pub lat_e8: i64,
    pub lon_e8: i64,
    /// meters
    pub elevation_msl: f32,
    /// m/s east, north, up
    pub velocity: Vec3,
}

impl PositVelEarthUnits {
    /// Apply dead-reckoning to the most recent GNSS fix. This is what we fuse with inertial data.
    /// Returns a result with lat, lon, and alt as f32.
    /// time is in seconds.
    /// Returns lat, lon (both 1e7), and alt in mm.
    pub fn from_fix_dr(fix: &Fix, timestamp_s: f32) -> Self {
        let dt = timestamp_s - fix.timestamp_s;

        // These are in mm/s
        let [v_north, v_east, v_down] = fix.ned_velocity;

        let (lat, lon) = augment_latlon(
            (fix.lat as i64) * 10,
            (fix.lon as i64) * 10,
            v_north as f32 * dt / 1_000.,
            v_east as f32 * dt / 1_000.,
        );

        Self {
            lat_e8: lat,
            lon_e8: lon,
            elevation_msl: (fix.elevation_msl as f32 / 1_000.) - (v_down as f32 * dt / 1_000.),
        }
    }
}

/// Modify a latitude and longitude, with a distance to add in meters. Input and output are in degrees x 1e7.
/// This format matches what we receive from the GNSS.
/// Per this article (https://en.wikipedia.org/wiki/Decimal_degrees), this allows for a precision of 11mm
/// at the equator.
/// lat and lon are in degrees x 1e8. to_add is in m.
/// We use this combined function for lat and lon, since latitude is required to convert meters to
/// degrees longitude.
/// Note: This only roughly takes into account the non-spherical character of the earth.
fn augment_latlon(lat_e8: i64, lon_e8: i64, add_lat_m: f32, add_lon_m: f32) -> (i64, i64) {
    const M_PER_DEG_LAT_EQUATOR: f32 = 1_843. * 60.;
    const M_PER_DEG_LAT_POLE: f32 = 1_861. * 60.;

    let port_through_lat = (lat_e8 as f32).abs() / DEG_SCALE_1E8 / 90.; // todo: Use cosine instead?

    let m_per_deg_lat =
        M_PER_DEG_LAT_EQUATOR + port_through_lat * (M_PER_DEG_LAT_POLE - M_PER_DEG_LAT_EQUATOR);

    let result_lat = lat_e8 + (add_lat_m / m_per_deg_lat * DEG_SCALE_1E8) as i64;

    // todo: What's a good avg value for lon at eq?
    let m_per_deg_lon = (result_lat as f32 / DEG_SCALE_1E8).cos() * M_PER_DEG_LAT_EQUATOR;

    let result_lon = lon_e8 + (add_lon_m / m_per_deg_lon * DEG_SCALE_1E8) as i64;

    (result_lat, result_lon)
}
