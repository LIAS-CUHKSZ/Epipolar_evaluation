#ifndef MANIFOLD_GN_HPP
#define MANIFOLD_GN_HPP

#include "Eigen/src/Core/Matrix.h"
#include "PoseSovlerWrapper.hpp"
#include "eigen3/Eigen/Eigen"
#include "eigen3/unsupported/Eigen/KroneckerProduct"
#include "eigen3/unsupported/Eigen/MatrixFunctions"
#include <opencv2/calib3d.hpp>
#include <opencv2/core/eigen.hpp>
#include <utils.hpp>
#include <vector>

using namespace std;
using namespace cv;
using namespace Eigen;

#define kron kroneckerProduct

class ManifoldGN : public PoseSolver {
private:
  Matrix3d intrinsic;
  int m;

public:
  ManifoldGN() = default;

  string getName() const override { return "ManifoldGN"; }

  Trans_type GetPose(Matrix3d R_init, Vector3d t_init,
                     Matrix3d &E_out) override {
    Matrix3d E_init;
    Mat intri_cv, es_cv, R_cv, t_cv;
    if (R_init.isZero()) {
      // Initialize using OpenCV if no initial guess provided
      eigen2cv(intrinsic, intri_cv);
      es_cv =
          findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, RANSAC, 0.999, 1.0);
      recoverPose(es_cv, y_cv_pix, z_cv_pix, intri_cv, R_cv, t_cv);
      cv2eigen(R_cv, R_init);
      cv2eigen(t_cv, t_init);
      cv2eigen(es_cv, E_init);
    } else {
      E_init = skew(t_init) * R_init;
    }

    // SVD decomposition
    Matrix3d U, E0, V;
    auto svd_res = E_init.jacobiSvd(ComputeFullU | ComputeFullV);
    U = svd_res.matrixU();
    V = svd_res.matrixV();
    E0 = svd_res.singularValues().asDiagonal();

    // Ensure proper rotation matrices
    if (U.determinant() == -1)
      U.col(2) = -U.col(2);
    if (V.determinant() == -1)
      V.col(2) = -V.col(2);

    // Construct measurement matrix
    MatrixXd bar_M(m, 9);
    for (int i = 0; i < m; i++) {
      bar_M.row(i) = (z_n[i] * y_n[i].transpose()).reshaped(1, 9);
    }
    Matrix<double, 9, 9> m_M = bar_M.transpose() * bar_M / m;

    // Construct generators of so(3)
    Vector<double, 9> Qx, Qy, Qz;
    Qx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
    Qy << 0, 0, -1, 0, 0, 0, 1, 0, 0;
    Qz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
    Qx /= sqrt(2);
    Qy /= sqrt(2);
    Qz /= 2;

    // Construct Q matrices
    Matrix<double, 9, 5> Q1, Q2;
    Matrix<double, 9, 2> zero_stand = Matrix<double, 9, 2>::Zero();
    Q1 << Qx, Qy, Qz, zero_stand;
    Q2 << zero_stand, -Qz, Qx, Qy;

    // Gauss-Newton step
    MatrixXd J = kron(V, U) * (kron(E0, Matrix3d::Identity()) * Q1 -
                               kron(Matrix3d::Identity(), E0) * Q2);
    MatrixXd first_deri = J.transpose() * m_M * E_init.reshaped(9, 1);
    MatrixXd second_deri = J.transpose() * m_M * J;
    MatrixXd x = -second_deri.inverse() * first_deri;

    // Convert to Lie algebra and apply exponential map
    Matrix3d omega1, omega2;
    omega1 << 0, -x(2) / sqrt(2), x(1), x(2) / sqrt(2), 0, -x(0), -x(1), x(0),
        0;
    omega2 << 0, x(2) / sqrt(2), x(4), -x(2) / sqrt(2), 0, -x(3), -x(4), x(3),
        0;
    omega1 /= sqrt(2);
    omega2 /= sqrt(2);

    // Update essential matrix
    Matrix3d U_GN = U * omega1.exp();
    Matrix3d V_GN = V * omega2.exp();
    Matrix3d E_GN = U_GN * E0 * V_GN.transpose();

    // Recover final pose
    Mat E_GN_cv, R_GN_cv, t_GN_cv;
    eigen2cv(E_GN, E_GN_cv);
    recoverPose(E_GN_cv, y_cv_pix, z_cv_pix, intri_cv, R_GN_cv, t_GN_cv);

    Trans_type result;
    cv2eigen(R_GN_cv, result.Rotation);
    cv2eigen(t_GN_cv, result.Translation);
    return result;
  }
};

#endif // MANIFOLD_GN_HPP