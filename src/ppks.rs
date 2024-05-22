//! Present-position keeping system. Fuses GPS, and dead-reckoning.
//!
//! Sync Corvus' implemention from this.

// todo: Move to ahrs module, but we need to sort out how to handle `Fix`.

use lin_alg2::f32::{Quaternion, Vec3};
use num_traits::Float;

use crate::{Fix, DEG_SCALE_1E8};

pub const FIX_FUSED_SIZE: usize = 8 * 3 + 4 * 5;

use defmt::println;

/// Convert NED velocity in mm/s to xyz velocity in m/s. (Still in earth frame.)
pub fn ned_vel_to_xyz(ned_vel: [i32; 3]) -> Vec3 {
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
    /// Fuse dead-reckoning position and velocity from the most-recent fix with inertial
    /// position and velocity.
    /// `fix` is the most-recently-received fix.
    pub fn new(
        // todo: fix here vs inertial's anchor? Likely the same.
        fix: &Fix,
        inertial: &mut PositInertial,
        timestamp: u64,
    ) -> Self {
        let gnss_dr = PositVelEarthUnits::from_fix_dr(fix, timestamp);

        // todo: If dt since last fix is longer, weigh INS more.

        // todo: Update heading from Inertial + GNSS accelerations (Change in NSE velocity).

        let inertial_absolute = inertial.combine_with_anchor();

        // For this reference, we update gnss DR with our inertial position and velocity
        let update_amount = 0.5;
        let update_amount = 1.;

        let add_to_gnss = PositVelEarthUnits {
            // todo: This /2 is a rough proxy for update amount on fixed-point.
            lat_e8: (inertial_absolute.lat_e8 - gnss_dr.lat_e8) / 2,
            lon_e8: (inertial_absolute.lon_e8 - gnss_dr.lon_e8) / 2,
            elevation_msl: (inertial_absolute.elevation_msl - gnss_dr.elevation_msl)
                * update_amount,
            velocity: (inertial_absolute.velocity - gnss_dr.velocity) * update_amount,
        };

        // todo: Temp, until inertial works.
        let add_to_gnss = PositVelEarthUnits {
            // todo: This /2 is a rough proxy for update amount on fixed-point.
            lat_e8: 0,
            lon_e8: 0,
            elevation_msl: 0.,
            velocity: Vec3::new_zero(),
        };

        let vel = gnss_dr.velocity + add_to_gnss.velocity;

        let elevation_msl = gnss_dr.elevation_msl + add_to_gnss.elevation_msl;

        let elevation_hae = (fix.elevation_hae as f32 / 1_000.)
            + (elevation_msl - (fix.elevation_msl as f32 / 1_000.));

        Self {
            timestamp_us: timestamp * 1_000,
            lon_e8: gnss_dr.lon_e8 + add_to_gnss.lon_e8,
            lat_e8: gnss_dr.lat_e8 + add_to_gnss.lat_e8,
            elevation_msl,
            elevation_hae,
            ned_velocity: [vel.x, vel.y, vel.z],
        }
    }

    /// We use this as our CAN wire format for fused fixes.
    pub fn to_bytes(&self) -> [u8; FIX_FUSED_SIZE] {
        let mut result = [0; FIX_FUSED_SIZE];

        result[0..8].copy_from_slice(&self.timestamp_us.to_le_bytes());
        result[8..16].copy_from_slice(&self.lat_e8.to_le_bytes());
        result[16..24].copy_from_slice(&self.lon_e8.to_le_bytes());
        result[24..28].copy_from_slice(&self.elevation_hae.to_le_bytes());
        result[28..32].copy_from_slice(&self.elevation_msl.to_le_bytes());

        result[32..36].copy_from_slice(&self.ned_velocity[0].to_le_bytes());
        result[36..40].copy_from_slice(&self.ned_velocity[1].to_le_bytes());
        result[40..44].copy_from_slice(&self.ned_velocity[2].to_le_bytes());

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
        acc_lin: Vec3,
        attitude: Quaternion,
        dt: f32, // seconds
    ) {
        // todo: temp hard Hard-codefor testing.
        // self.anchor = PositVelEarthUnits {
        //     lat_e8: 3500000000,
        //     lon_e8: -7800000000,
        //     elevation_msl: 0.,
        //     velocity: Vec3::new(0., 1., 0.),
        // };

        // todO: For debugging
        // let acc_lin = Vec3::new(
        //     0.,
        //     0.1,
        //     0.
        // );

        self.posit += self.velocity * dt;
        self.velocity += attitude.inverse().rotate_vec(acc_lin) * dt;

        static mut I: u32 = 0;
        unsafe { I += 1 };

        let combined = self.combine_with_anchor();

        // if unsafe { I } % 2000 == 0 {
        if false {
            println!(
                "\n\nInertial: x{} y{} z{} -- vx{} vy{} vz{}",
                self.posit.x,
                self.posit.y,
                self.posit.z,
                self.velocity.x,
                self.velocity.y,
                self.velocity.z
            );

            println!(
                "Combined: lat{} lon{} msl{} -- vx{} vy{} vz{}",
                combined.lat_e8,
                combined.lon_e8,
                combined.elevation_msl,
                combined.velocity.x,
                combined.velocity.y,
                combined.velocity.z
            );
        }
    }

    /// Update the anchor point, based on a new fix.
    pub fn update_anchor(&mut self, fix: &Fix) {
        // todo: For now, we accept the whole fix. In the future, weight with our inertial position.
        // let diff = ned_vel_to_xyz(fix.ned_velocity) - self.velocity;
        //
        // // let update_factor = 0.8; // 1.0 means replace the anchor entirely.
        // let update_factor = 1.; // 1.0 means replace the anchor entirely.

        self.anchor = PositVelEarthUnits {
            lat_e8: fix.lat_e7 as i64 * 10,
            lon_e8: fix.lon_e7 as i64 * 10,
            elevation_msl: fix.elevation_msl as f32 / 1_000.,
            velocity: ned_vel_to_xyz(fix.ned_velocity),
        };

        // For now, completely sync to anchor by resetting our relative position and velocity.
        // todo: In the future, don't do this
        self.posit = Vec3::new_zero();
        self.velocity = Vec3::new_zero();
    }

    /// Combine the relative units stored by this struct's position and velocity with its anchor.
    pub fn combine_with_anchor(&self) -> PositVelEarthUnits {
        let (lat, lon) = augment_latlon(
            self.anchor.lat_e8,
            self.anchor.lon_e8,
            self.posit.y,
            self.posit.x,
        );

        PositVelEarthUnits {
            lat_e8: lat,
            lon_e8: lon,
            elevation_msl: self.anchor.elevation_msl + self.posit.z,
            velocity: self.anchor.velocity + self.velocity,
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
    pub fn from_fix(fix: &Fix) -> Self {
        Self {
            lat_e8: (fix.lat_e7 as i64) * 10,
            lon_e8: (fix.lon_e7 as i64) * 10,
            elevation_msl: fix.elevation_msl as f32 / 1_000.,
            velocity: ned_vel_to_xyz(fix.ned_velocity),
        }
    }

    /// Apply dead-reckoning to the most recent GNSS fix.
    /// Returns a result with lat, lon, and alt as f32.
    /// time is in seconds.
    /// Returns lat, lon (both 1e7), and alt in mm.
    pub fn from_fix_dr(fix: &Fix, timestamp: u64) -> Self {
        // todo: Is this fix precision acceptable?
        let dt = (timestamp as f32 / 1_000.)  - fix.timestamp_s;

        let velocity = ned_vel_to_xyz(fix.ned_velocity);

        let (lat_e8, lon_e8) = augment_latlon(
            (fix.lat_e7 as i64) * 10,
            (fix.lon_e7 as i64) * 10,
            velocity.y * dt,
            velocity.x * dt,
        );

        Self {
            lat_e8,
            lon_e8,
            elevation_msl: (fix.elevation_msl as f32 / 1_000.) + velocity.z * dt,
            velocity,
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
