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
use na::{Matrix3, Matrix4, RowVector3, RowVector4};

use num_traits::Float;

use crate::attitude::MAG_CAL_DATA_LEN;

use lin_alg2::f32::{Mat3, Mat4, Vec3};

use defmt::println;

/// least squares fit to a 3D-ellipsoid
/// Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz  = 1
///
/// Note that sometimes it is expressed as a solution to
///  Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz  = 1
/// where the last six terms have a factor of 2 in them
/// This is in anticipation of forming a matrix with the polynomial coefficients.
/// Those terms with factors of 2 are all off diagonal elements.  These contribute
/// two terms when multiplied out (symmetric) so would need to be divided by two
fn ls_ellipsoid(
    x: &[f32; MAG_CAL_DATA_LEN],
    y: &[f32; MAG_CAL_DATA_LEN],
    z: &[f32; MAG_CAL_DATA_LEN],
) -> [f32; 10] {
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    // np.hstack performs a loop over all samples and creates
    // a row in J for each x,y,z sample:
    // J[ix,0] = x[ix]*x[ix]
    // J[ix,1] = y[ix]*y[ix]
    // etc.
    // let mut J = [[0.; 9]; MAG_CAL_DATA_LEN];
    // for i in 0..MAG_CAL_DATA_LEN {
    //     let (xi, yi, zi) = (x[i], y[i], z[i]);
    //     J[i] = [xi * xi, yi * yi, zi * zi, xi * yi, xi * zi, yi * zi, xi, yi, zi];
    // }

    // let J_storate = nalgebra::base::ArrayStorage<f32, MAG_CAL_DATA_LEN, 9>

    // let mut J = na::DMatrix::from_diagonal_elements(MAG_CAL_DATA_LEN, 9, 1.);
    let j_data = [0.0_f32; MAG_CAL_DATA_LEN * 9];
    // let mut J = na::MatrixView::from_slice_with_strides_generic(
    //     &j_data,
    //     9,
    //     MAG_CAL_DATA_LEN,
    //     9 * 8,
    //     MAG_CAL_DATA_LEN * 8,
    // );

    let j_storage = na::ArrayStorage::new(j_data);
    let mut J: na::SMatrix<f32, MAG_CAL_DATA_LEN, 9> = na::SMatrix::from_array_storage(j_storage);

    for i in 0..MAG_CAL_DATA_LEN {
        let (xi, yi, zi) = (x[i], y[i], z[i]);
        J[i] = [
            xi * xi,
            yi * yi,
            zi * zi,
            xi * yi,
            xi * zi,
            yi * zi,
            xi,
            yi,
            zi,
        ];
    }

    let K = [1.; MAG_CAL_DATA_LEN]; // column of ones

    let JT = J.transpose();
    let JTJ = JT.dot(J);
    let InvJTJ = JTJ.inverse();
    let ABC = InvJTJ * (JT.dot(K));

    // Rearrange, move the 1 to the other side
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    //    or
    //  Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    //  where J = -1
    let mut result = [0.; 10];
    result[0..9].copy_from_slice(ABC[0..9]);
    result[9] = -1.; // J term

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
fn poly_to_params_3d(coeffs: &[f32; 10]) -> (Vec3, Vec3, Mat3) {
    // #[rustfmt::skip]
    // let amat = Mat4::new([
    //     vec[0], vec[3]/2., vec[4]/2., vec[6]/2.,
    //     vec[3]/2., vec[1], vec[5]/2., vec[7] / 2.,
    //     vec[4]/2., vec[5]/2., vec[2], vec[8] / 2.,
    //     vec[6]/2., vec[7]/2., vec[8]/2., vec[9],
    // ]);

    let c = &coeffs;
    let Amat = na::Matrix4::from_rows(&[
        RowVector4::new(c[0], c[3] / 2.0, c[4] / 2.0, c[6] / 2.0),
        RowVector4::new(c[3] / 2.0, c[1], c[5] / 2.0, c[7] / 2.0),
        RowVector4::new(c[4] / 2.0, c[5] / 2.0, c[2], c[8] / 2.0),
        RowVector4::new(c[6] / 2.0, c[7] / 2.0, c[8] / 2.0, c[9]),
    ]);

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

    let A3 = Matrix3::from_rows(&[
        RowVector3::new(c[0], c[3] / 2.0, c[4] / 2.0),
        RowVector3::new(c[3] / 2.0, c[1], c[5] / 2.0),
        RowVector3::new(c[4] / 2.0, c[5] / 2.0, c[2]),
    ]);

    let A3_inv = A3.inverse();

    let ofs = Vec3::new(c[6], c[7], c[8]) / 2.;
    let center = -(A3_inv.dot(ofs));

    println!("Center: {:?}", center);

    // Center the ellipsoid at the origin
    // let mut tofs = Mat4::new_identity();

    let Tofs = Matrix4::identity();
    // Tofs.data[12] = center.x;
    // Tofs.data[13] = center.y;
    // Tofs.data[14] = center.z;
    Tofs[3] = [center.x, center.y, center.z];

    let R = tofs.dot(Amat.dot(Tofs.transpose()));
    // if printMe: print '\nAlgebraic form translated to center\n',R,'\n'

    // let rd = &R.data;
    // #[rustfmt::skip]
    // let R3 = Mat3::new([
    //     rd[0],rd[1], rd[2],
    //     rd[4], rd[5], rd[6],
    //     rd[8], rd[9], rd[10]
    // ]);

    let R3 = Matrix3::from_rows(&[
        RowVector3::new(R.row[0][0], R.row[0][1], R.row[0][2]),
        RowVector3::new(R.row[1][0], R.row[1][1], R.row[1][2]),
        RowVector3::new(R.row[2][0], R.row[2][1], R.row[2][2]),
    ]);

    let R3_test = R3 / R3.row[0][0];
    // print 'normed \n',R3test

    let s1 = -R.row[3][3];
    let R3S = R3 / s1;

    let eigen = R3S.symmetric_eigen();
    let (eigen_vals, eigen_vecs) = (eigen.eigenvalues, eigen.eigenvectors);

    let recip = 1.0_f32 / eigen_vals.abs();
    let axes = recip.sqrt();
    // if printMe: print '\nAxes are\n',axes  ,'\n'

    let inve = eigen_vecs.inverse(); // inverse is actually the transpose here
                                     // if printMe: print '\nRotation matrix\n',inve
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
