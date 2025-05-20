#![allow(non_snake_case)]

//! http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
//! http://www.juddzone.com/ALGORITHMS/least_squares_precision_3D_ellipsoid.html
//! Has detailed examples, but unfortunately they rely on lin alg ops like inverse
//! on n-degree matrices.
//!
//! See also: https://github.com/PaulStoffregen/MotionCal/blob/master/magcal.c
//!
//! the `nalgebra` crate could simplify our matrix operations, but it increases binary size
//! too much.

// todo: Good approach, but you may be limited by embedded here on the matrix inverse. Maybe
// todo loop in nalgebra

use lin_alg::f32::{Mat3, Vec3};
use na::{Matrix3, Matrix4, RowVector4, SMatrix, SVector};
use nalgebra as na;
use num_traits::Float;

// todo: In addition to fitting hard and soft iron offsets from ellipsoids,
// todo: You can use the gyro to factor into these.

/// Vectors of attitudes the craft is at while taking elipsoid sample points. Make sure there
/// are several points near each of these prior to calibrating.
/// Dodecahedron vertices.
/// https://math.fandom.com/wiki/Dodecahedron
const PHI: f32 = 1.618_034; // 1/2 + (5/4).sqrt()

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
pub const SAMPLE_VERTEX_ANGLE: f32 = 0.729_727_7;

// We need at least this many samples per category before calibrating.
// Note: We'd ideally like this to be higher, but are experiencing stack overflows eg at 5.
pub const MAG_SAMPLES_PER_CAT: usize = 2;

// // times 2, since we're comparing to perpendicular vectors to each point.
pub const TOTAL_MAG_SAMPLE_PTS: usize = MAG_SAMPLES_PER_CAT * SAMPLE_VERTICES.len() * 2;

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

    let mut J: SMatrix<f32, { TOTAL_MAG_SAMPLE_PTS }, 9> = SMatrix::zeros();

    for i in 0..TOTAL_MAG_SAMPLE_PTS {
        let (xi, yi, zi) = (sample_pts[i].x, sample_pts[i].y, sample_pts[i].z);

        // A sample pt value of 0 means no data on that point.
        const EPS: f32 = 0.0000001;
        if xi.abs() < EPS && yi.abs() < EPS && zi.abs() < EPS {
            continue;
        }

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
    let JTJ = JT * J; // 9x9 matrix.
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
    [
        ABC[0], ABC[1], ABC[2], ABC[3], ABC[4], ABC[5], ABC[6], ABC[7], ABC[8], -1.,
    ]
}

/// http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
/// convert the polynomial form of the 3D-ellipsoid to parameters
/// center, axes, and transformation matrix
/// vec is the vector whose elements are the polynomial
/// coefficients A..J
/// returns (center, axes, rotation matrix)
///
/// Algebraic form: X.T * Amat * X --> polynomial form
pub fn poly_to_params_3d(coeffs: &[f32; 10]) -> (Vec3, Mat3) {
    // #[rustfmt::skip]
    // let amat = Mat4::new([
    //     vec[0], vec[3]/2., vec[4]/2., vec[6]/2.,
    //     vec[3]/2., vec[1], vec[5]/2., vec[7] / 2.,
    //     vec[4]/2., vec[5]/2., vec[2], vec[8] / 2.,
    //     vec[6]/2., vec[7]/2., vec[8]/2., vec[9],
    // ]);

    let c = &coeffs;
    #[rustfmt::skip]
        let Amat = Matrix4::new(
        c[0],c[3]/2.,c[4]/2.,c[6]/2.,
        c[3]/2.,c[1],c[5]/2.,c[7]/2.,
        c[4]/2.,c[5]/2.,c[2],c[8]/2.,
        c[6]/2.,c[7]/2.,c[8]/2.,c[9],
    );

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

    #[rustfmt::skip]
        let A3 = Matrix3::new(
        c[0],c[3]/2.,c[4]/2.,
        c[3]/2.,c[1],c[5]/2.,
        c[4]/2.,c[5]/2.,c[2],
    );

    let A3_inv = A3.try_inverse().unwrap_or(Matrix3::identity());

    // let ofs = Vec3::new(c[6], c[7], c[8]) / 2.;
    let ofs = SVector::from([c[6], c[7], c[8]]) / 2.;
    let center = -(A3_inv * ofs);

    // Center the ellipsoid at the origin
    // let mut tofs = Mat4::new_identity();

    let mut Tofs = Matrix4::identity();
    Tofs.set_row(3, &RowVector4::new(center[0], center[1], center[2], 1.));

    let R = Tofs * (Amat * Tofs.transpose());
    // if printMe: print '\nAlgebraic form translated to center\n',R,'\n'

    // let rd = &R.data;

    #[rustfmt::skip]
        let R3 = Matrix3::new(
        R[(0, 0)],R[(0, 1)],R[(0, 2)],
        R[(1, 0)],R[(1, 1)],R[(1, 2)],
        R[(2, 0)],R[(2, 1)],R[(2, 2)],
    );

    let s1 = -R[(3, 3)];
    let R3S = R3 / s1;

    let eigen = R3S.symmetric_eigen();
    let (eigen_vals, eigen_vecs) = (eigen.eigenvalues, eigen.eigenvectors);

    let axes = Vec3::new(
        (1. / eigen_vals[0].abs()).sqrt(),
        (1. / eigen_vals[1].abs()).sqrt(),
        (1. / eigen_vals[2].abs()).sqrt(),
    );

    // Inverse is actually the transpose here
    let Rot = eigen_vecs.try_inverse().unwrap_or(Matrix3::identity());

    // L = np.diag([1/axes[0],1/axes[1],1/axes[2]])
    // M=np.dot(R.T,np.dot(L,R))

    let L = Matrix3::new(
        1. / axes.x,
        0.,
        0.,
        0.,
        1. / axes.y,
        0.,
        0.,
        0.,
        1. / axes.z,
    );

    let M = Rot.transpose() * (L * Rot);

    // Our constructor is column-major. nalgebra indices are row-first.
    #[rustfmt::skip]
        let M = Mat3::new([
        M[(0, 0)], M[(1, 0)], M[(2, 0)],
        M[(0, 1)], M[(1, 1)], M[(2, 1)],
        M[(0, 2)], M[(1, 2)], M[(2, 2)],
    ]);

    let center = Vec3::new(center[0], center[1], center[2]);

    // ie hard-iron and soft-iron
    (center, M)
}

// "Iterative" appraoch below

// todo: If this approach works and is performant enough, remove the nalg dependency and above code.

// fn mag_dot_acc_err(mag: Vec3,acc: Vec3, mdg: Vec3, params: &[f32; 12]) -> Vec3 {
//     // offset and transformation matrix from parameters
//     let ofs = Vec3::new(params[0], params[1], params[2]);
//
//     // let mat = np.reshape(params[3..12],(3,3));
//     // Col-major.
//     #[rustfmt::skip]
//     let mat = Mat3::new([
//         params[3], params[6], params[9],
//         params[4], params[7], params[10],
//         params[5], params[8], params[11],
//     ]);
//
//     // subtract offset, then apply transformation matrix
//     let mc= mag - ofs;
//
//     let mm = mat.dot(mc);
//
//     // calculate dot product from corrected mags
//     let mdg1 = mm.dot(acc);
//
//     mdg-mdg1
//
// }
//
// fn analytic_partial_row(mag: Vec3,acc: Vec3, target: Vec3, params: &[f32; 12]) -> (Vec3, [f32; 12]){
//     let err0 = magDotAccErr(mag, acc, target, params);
//     // ll = len(params)
//     let mut slopeArr = [0.; 12];
//
//     slopeArr[0] = -(params[3] * acc.x + params[4] * acc.y + params[5] * acc.z);
//     slopeArr[1] = -(params[6] * acc.x + params[7] * acc.y + params[8] * acc.z);
//     slopeArr[2] = -(params[9] * acc.x + params[10] * acc.y + params[11] * acc.z);
//
//     slopeArr[3] = (mag.x - params[0]) * acc.x;
//     slopeArr[4] = (mag.y - params[1]) * acc.x;
//     slopeArr[5] = (mag.z - params[2]) * acc.x;
//
//     slopeArr[6] = (mag.x - params[0]) * acc.y;
//     slopeArr[7] = (mag.y - params[1]) * acc.y;
//     slopeArr[8] = (mag.z - params[2]) * acc.y;
//
//     slopeArr[9] = (mag.x - params[0]) * acc.z;
//     slopeArr[10] = (mag.y - params[1]) * acc.z;
//     slopeArr[11] = (mag.z - params[2]) * acc.z;
//
//     (err0, slopeArr)
// }
//
// fn numericPartialRow(mag: Vec3,acc: Vec3, target: Vec3, params: &mut [f32; 12], step: usize) -> (Vec3, [f32; 12]) {
//     let err0 = mag_dot_acc_err(mag, acc, target, params);
//
//     // let ll = params.len();
//     const ll: usize = 12;
//     let mut slopeArr = [0.; ll];
//
//     for ix in 0..ll {
//         params[ix] = params[ix] + step[ix];
//         errA = mag_dot_acc_err(mag, acc, target, params);
//
//         params[ix] = params[ix] - 2.0 * step[ix];
//         errB = mag_dot_acc_err(mag, acc, target, params);
//
//         params[ix] = params[ix] + step[ix];
//         slope = (errB - errA) / (2.0 * step[ix]);
//
//         slopeArr[ix] = slope;
//     }
//
//     (err0, slopeArr)
// }
//
//
// fn ellipsoid_iterate(mag: Vec3,accel: Vec3, verbose: bool) -> ([f32; 12], Vec3) {
//
//    let magCorrected= mag.clone();
//    // Obtain an estimate of the center and radius
//    // For an ellipse we estimate the radius to be the average distance
//    // of all data points to the center
//    let (centerE,magR,magSTD)= ellipsoid_estimate2(mag,verbose);
//
//    // Work with normalized data
//    let magScaled=mag/magR;
//    let centerScaled = centerE/magR;
//
//    let accNorm = accel.to_normalized();
//
//    let params = [0.; 12];
//    // Use the estimated offsets, but our transformation matrix is unity
//    params[0..3].clone_from_slice(&[centerScaled.x, centerScaled.y, centerScaled.z]);
//
//     let mat = Mat3::new_identity();
//
//     // todo: Probably needs to be recast as row-order to make this work.
//    params[3..12].copy_from_slice(mat3.data);
//
//    // initial dot based on centered mag, scaled with average radius
//    let magCorrected=applyParams12(magScaled,params);
//    let (avgDot,stdDot)=mgDot(magCorrected,accNorm);
//
//    let nSamples= magScaled.len();
//    let sigma = errorEstimate(magScaled,accNorm,avgDot,params);
//    if verbose {
//        println!("Initial sigma: {}", sigma);
//    }
//
//    // pre allocate the data.  We do not actually need the entire
//    // D matrix if we calculate DTD (a 12x12 matrix) within the sample loop
//    // Also DTE (dimension 12) can be calculated on the fly.
//
//     // let D = [0.]
//    let D = np.zeros([nSamples,12]);
//    let E= [0.; nSAMPLES];
//
//    // If numeric derivatives are used, this step size works with normalized data.
//    let step = [1./5_000.; 12];
//
//    // Fixed number of iterations for testing.  In production you check for convergence
//
//    let nLoops= 5;
//
//    for iloop in 0..nLoops {
//       // Numeric or analytic partials each give the same answer
//       for ix in 0..nSamples {
//           // (f0,pdiff)=numericPartialRow(magScaled[ix],accNorm[ix],avgDot,params,step,1)
//           (f0, pdiff) = analyticPartialRow(magScaled[ix], accNorm[ix], avgDot, params)
//           E[ix] = f0;
//           D[ix] = pdiff;
//       }
//       // Use the pseudo-inverse
//       DT=D.T;
//       DTD=np.dot(DT,D);
//       DTE=np.dot(DT,E);
//       invDTD=inv(DTD);
//       deltas=np.dot(invDTD,DTE);
//
//
//       p2=params + deltas;
//
//       sigma = errorEstimate(magScaled,accNorm,avgDot,p2);
//
//       // add some checking here on the behavior of sigma from loop to loop
//       // if satisfied, use the new set of parameters for the next iteration
//
//       params=p2;
//
//       // recalculste gain (magR) and target dot product
//       let magCorrected=applyParams12(magScaled,params);
//       let (mc,mcR)=normalize3(magCorrected);
//       let (avgDot,stdDot)=mgDot(mc,accNorm);
//       magR *= mcR;
//       let magScaled=mag/magR;
//
//       if verbose {
//           println!("iloop: {}, sigma: {}", iloop, sigma);
//       }
//     }
//
//    (params,magR)
//
