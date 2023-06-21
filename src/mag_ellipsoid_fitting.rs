#![allow(non_snake_case)]

//! http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
//! Has detailed examples, but unfortunately they rely on lin alg ops like inverse
//! on n-degree matrices.
//!
//! See also: https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c

// todo: Good approach, but you may be limited by embedded here on the matrix inverse. Maybe
// todo loop in nalgebra

use nalgebra as na;
// use ndarray;
use na::{Const, Matrix3, Matrix4, RowVector4, SMatrix, SVector};

use num_traits::Float;

use lin_alg2::f32::{Mat3, Vec3};

use defmt::println;

// todo: In addition to fitting hard and soft iron offsets from ellipsoids,
// todo: You can use the gyro to factor into these.

/// Vectors of attitudes the craft is at while taking elipsoid sample points. Make sure there
/// are several points near each of these prior to calibrating.
/// Dodecahedron vertices.
/// https://math.fandom.com/wiki/Dodecahedron
const PHI: f32 = 1.6180339887; // 1/2 + (5/4).sqrt()

// todo: If this doesn't compile:
// use core::cell::OnceCell;
// static SAMPLE_VECS = OnceCell::new();
// SAMPLE_VECS.set().unwrap();
// todo: Normalize these!
pub const SAMPLE_VERTICES: [Vec3; 20] = [
    Vec3::new(1., 1., 1.),
    Vec3::new(1., 1., -1.),
    Vec3::new(1., -1., 1.),
    Vec3::new(1., -1., -1.),
    //
    Vec3::new(-1., 1., 1.),
    Vec3::new(-1., 1., -1.),
    Vec3::new(-1., -1., 1.),
    Vec3::new(-1., -1., -1.),
    //
    Vec3::new(0., 1. / PHI, PHI),
    Vec3::new(0., 1. / PHI, -PHI),
    Vec3::new(0., -1. / PHI, PHI),
    Vec3::new(0., -1. / PHI, -PHI),
    //
    Vec3::new(1. / PHI, PHI, 0.),
    Vec3::new(1. / PHI, -PHI, 0.),
    Vec3::new(-1. / PHI, PHI, 0.),
    Vec3::new(-1. / PHI, -PHI, 0.),
    //
    Vec3::new(PHI, 0., 1. / PHI),
    Vec3::new(PHI, 0., -1. / PHI),
    Vec3::new(-PHI, 0., 1. / PHI),
    Vec3::new(-PHI, 0., -1. / PHI),
];

// The angular distance between neighboring sample vecs. If a quaternion's *up* vec (We chose this
// arbitrarily; any vec will do) is closer than this to a sample vec, we group it with that sample vec.
pub const SAMPLE_VERTEX_ANGLE: f32 = 0.7297276562166872;

// We need at least this many samples per category before calibrating.
pub const MAG_SAMPLES_PER_CAT: usize = 10;
pub const TOTAL_MAG_SAMPLE_PTS: usize = MAG_SAMPLES_PER_CAT * SAMPLE_VERTICES.len();

/// least squares fit to a 3D-ellipsoid
/// Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1
///
/// Note that sometimes it is expressed as a solution to
///  Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
/// where the last six terms have a factor of 2 in them
/// This is in anticipation of forming a matrix with the polynomial coefficients.
/// Those terms with factors of 2 are all off diagonal elements.  These contribute
/// two terms when multiplied out (symmetric) so would need to be divided by two
pub fn ls_ellipsoid(sample_pts: &[Vec3; TOTAL_MAG_SAMPLE_PTS]) -> [f32; 10] {
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    // np.hstack performs a loop over all samples and creates
    // a row in J for each x,y,z sample:
    // J[ix,0] = x[ix]*x[ix]
    // J[ix,1] = y[ix]*y[ix]
    // etc.
    // let mut J = [[0.; 9]; TOTAL_MAG_SAMPLE_PTS];
    // for i in 0..TOTAL_MAG_SAMPLE_PTS {
    //     let (xi, yi, zi) = (x[i], y[i], z[i]);
    //     J[i] = [xi * xi, yi * yi, zi * zi, xi * yi, xi * zi, yi * zi, xi, yi, zi];
    // }

    // let J_storate = nalgebra::base::ArrayStorage<f32, TOTAL_MAG_SAMPLE_PTS, 9>

    // let mut J = na::DMatrix::from_diagonal_elements(TOTAL_MAG_SAMPLE_PTS, 9, 1.);
    // let j_data = [0.0_f32; TOTAL_MAG_SAMPLE_PTS * 9];
    // let mut J = na::MatrixView::from_slice_with_strides_generic(
    //     &j_data,
    //     9,
    //     TOTAL_MAG_SAMPLE_PTS,
    //     9 * 8,
    //     TOTAL_MAG_SAMPLE_PTS * 8,
    // );

    // let j_storage = ArrayStorage::new(j_data);
    // let mut J: SMatrix<f32, TOTAL_MAG_SAMPLE_PTS, 9> = SMatrix::from_array_storage(j_storage);
    //  Const<4>, Const<TOTAL_MAG_SAMPLE_PTS> etc?
    let mut J: SMatrix<f32, { TOTAL_MAG_SAMPLE_PTS }, 9> = SMatrix::zeros();

    for i in 0..TOTAL_MAG_SAMPLE_PTS {
        let (xi, yi, zi) = (sample_pts[i].x, sample_pts[i].y, sample_pts[i].z);
        J.set_row(
            i,
            &na::RowSVector::from([
                xi * xi,
                yi * yi,
                zi * zi,
                xi * yi,
                xi * zi,
                yi * zi,
                xi,
                yi,
                zi,
            ]),
        );
    }

    // column of ones
    let K: SVector<f32, { TOTAL_MAG_SAMPLE_PTS }> = SVector::from([1.; TOTAL_MAG_SAMPLE_PTS]);
    // let K = [1.; TOTAL_MAG_SAMPLE_PTS];

    let JT = J.transpose();
    let JTJ = JT * J;
    let InvJTJ = JTJ.try_inverse().unwrap();
    let ABC = InvJTJ * (JT * K);

    // Rearrange, move the 1 to the other side
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    //    or
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    //  where J = -1
    // let mut result = [0.; 10];
    // result[0..9].copy_from_slice(ABC[0..9]);
    // result[9] = -1.; // J term
    let result = [
        ABC[0], ABC[1], ABC[2], ABC[3], ABC[4], ABC[5], ABC[6], ABC[7], ABC[8], -1.,
    ];

    result
}

/// http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
/// convert the polynomial form of the 3D-ellipsoid to parameters
/// center, axes, and transformation matrix
/// vec is the vector whose elements are the polynomial
/// coefficients A..J
/// returns (center, axes, rotation matrix)
///
/// Algebraic form: X.T * Amat * X --> polynomial form
pub fn poly_to_params_3d(coeffs: &[f32; 10]) -> (Vec3, Vec3, Mat3) {
    // #[rustfmt::skip]
    // let amat = Mat4::new([
    //     vec[0], vec[3]/2., vec[4]/2., vec[6]/2.,
    //     vec[3]/2., vec[1], vec[5]/2., vec[7] / 2.,
    //     vec[4]/2., vec[5]/2., vec[2], vec[8] / 2.,
    //     vec[6]/2., vec[7]/2., vec[8]/2., vec[9],
    // ]);

    let c = &coeffs;
    let Amat = Matrix4::new(
        c[0],
        c[3] / 2.0,
        c[4] / 2.0,
        c[6] / 2.0,
        c[3] / 2.0,
        c[1],
        c[5] / 2.0,
        c[7] / 2.0,
        c[4] / 2.0,
        c[5] / 2.0,
        c[2],
        c[8] / 2.0,
        c[6] / 2.0,
        c[7] / 2.0,
        c[8] / 2.0,
        c[9],
    );

    // println!(""\nAlgebraic form of polynomial\n',Amat");"

    // See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    // equation 20 for the following method for finding the center
    // let A3 = Amat[0:3,0:3]
    // #[rustfmt::skip]
    //     let a3 = Mat3::new([
    //     vec[0],
    //     vec[3] / 2.,
    //     vec[4] / 2.,
    //     vec[3] / 2.,
    //     vec[1],
    //     vec[5] / 2.,
    //     vec[4] / 2.,
    //     vec[5],
    //     vec[2],
    // ]);
    //

    let A3 = Matrix3::new(
        c[0],
        c[3] / 2.0,
        c[4] / 2.0,
        c[3] / 2.0,
        c[1],
        c[5] / 2.0,
        c[4] / 2.0,
        c[5] / 2.0,
        c[2],
    );

    let A3_inv = A3.try_inverse().unwrap();

    // let ofs = Vec3::new(c[6], c[7], c[8]) / 2.;
    let ofs = SVector::from([c[6], c[7], c[8]]) / 2.;
    // let center = -(A3_inv.dot(ofs));
    let center = -(A3_inv * ofs);

    // println!("Center: {:?}", center);

    // Center the ellipsoid at the origin
    // let mut tofs = Mat4::new_identity();

    let mut Tofs = Matrix4::identity();
    // Tofs.data[12] = center.x;
    // Tofs.data[13] = center.y;
    // Tofs.data[14] = center.z;
    Tofs.set_row(3, &RowVector4::new(center[0], center[1], center[2], 1.));

    let R = Tofs * (Amat * &Tofs.transpose());
    // if printMe: print '\nAlgebraic form translated to center\n',R,'\n'

    // let rd = &R.data;
    // #[rustfmt::skip]
    // let R3 = Mat3::new([
    //     rd[0],rd[1], rd[2],
    //     rd[4], rd[5], rd[6],
    //     rd[8], rd[9], rd[10]
    // ]);
    let R3 = Matrix3::new(
        R[(0, 0)],
        R[(0, 1)],
        R[(0, 2)],
        R[(1, 0)],
        R[(1, 1)],
        R[(1, 2)],
        R[(2, 0)],
        R[(2, 1)],
        R[(2, 2)],
    );

    let R3_test = R3 / R3[(0, 0)];
    // print 'normed \n',R3test

    let s1 = -R[(3, 3)];
    let R3S = R3 / s1;

    let eigen = R3S.symmetric_eigen();
    let (eigen_vals, eigen_vecs) = (eigen.eigenvalues, eigen.eigenvectors);

    // let recip = 1.0_f32 / eigen_vals.abs();
    // let recip = [
    //     1. / eigen_vals[0].abs(),
    //     1. / eigen_vals[1].abs(),
    //     1. / eigen_vals[2].abs(),
    // ];
    let axes = Vec3::new(
        (1. / eigen_vals[0].abs()).sqrt(),
        (1. / eigen_vals[1].abs()).sqrt(),
        (1. / eigen_vals[2].abs()).sqrt(),
    );

    // let axes = recip.sqrt();
    // if printMe: print '\nAxes are\n',axes  ,'\n'

    let inve = eigen_vecs.try_inverse().unwrap(); // inverse is actually the transpose here
                                                  // if printMe: print '\nRotation matrix\n',inve

    let center = Vec3::new(center[0], center[1], center[2]);

    // Our constructor is column-major. nalgebra indices are row-first.
    // let inve = Mat3::new([
    //     inve[0][0], inve[1][0], inve[2][0], inve[0][1], inve[1][1], inve[2][1], inve[0][2],
    //     inve[1][2], inve[2][2],
    // ]);
    //
    let inve = Mat3::new_identity(); // todo temp; above not working yet

    (center, axes, inve)
}

//
// /// http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
// fn print_ans_3d(center,axes,R,xin,yin,zin,verbose) {
//
//     println!("\nCenter at  %10.4f,%10.4f,%10.4f" % (center[0],center[1],center[2]));
//     println!("Axes gains %10.4f,%10.4f,%10.4f " % (axes[0],axes[1],axes[2]));
//     println!("Rotation Matrix\n%10.5f,%10.5f,%10.5f\n%10.5f,%10.5f,%10.5f\n%10.5f,%10.5f,%10.5f" % (
//       R[0,0],R[0,1],R[0,2],R[1,0],R[1,1],R[1,2],R[2,0],R[2,1],R[2,2]));
//
//
//     // Check solution
//     // Convert to unit sphere centered at origin
//     //  1) Subtract off center
//     //  2) Rotate points so bulges are aligned with axes (no xy,xz,yz terms)
//     //  3) Scale the points by the inverse of the axes gains
//     //  4) Back rotate
//     // Rotations and gains are collected into single transformation matrix M
//
//     // subtract the offset so ellipsoid is centered at origin
//     xc=xin-center[0]
//     yc=yin-center[1]
//     zc=zin-center[2]
//
//     // create transformation matrix
//     L = np.diag([1/axes[0],1/axes[1],1/axes[2]])
//     M=np.dot(R.T,np.dot(L,R))
//     print '\nTransformation Matrix\n',M
//
//     // apply the transformation matrix
//     [xm,ym,zm]=np.dot(M,[xc,yc,zc])
//     # Calculate distance from origin for each point (ideal = 1.0)
//     rm = np.sqrt(xm*xm + ym*ym + zm*zm)
//
//     println!("\nAverage Radius  %10.4f (truth is 1.0)" % (np.mean(rm)));
//     println!("Stdev of Radius %10.4f\n " % (np.std(rm)));
// }