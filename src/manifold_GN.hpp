#ifndef MANIFOLD_GN_HPP
#define MANIFOLD_GN_HPP

#include "eigen3/Eigen/Eigen"
#include "eigen3/unsupported/Eigen/KroneckerProduct"
#include "eigen3/unsupported/Eigen/MatrixFunctions"
#include <vector>
#include <iostream>
#include "opencv2/calib3d.hpp"
#include "opencv2/core.hpp"
#include "opencv2/opencv.hpp"
#include "opencv2/core/eigen.hpp"
#include "opencv2/opencv_modules.hpp"
#include <cmath>

using namespace std;
using namespace cv;
using namespace Eigen;

class ManifoldGN
{

#define kron kroneckerProduct

private:
    Matrix3d intrinsic;
    Matrix3d R_GN;
    Vector3d t_GN;
    Matrix3d E_init, E_GN;
    int m;

public:
    ManifoldGN(Matrix3d intrinsic_,Matrix3d Essential_=Matrix3d::Identity()) : intrinsic(intrinsic_),E_init(Essential_){}

    void GetPose(Matrix3d &R, Vector3d &t, vector<Vector3d> &y_h, vector<Vector3d> &z_h)
    {
        /* calc the init E using OpenCV epipolar method, RANSAC */
        


        // convert vector3d to point3d type for OpenCV
        m = y_h.size();
        vector<Point3d> y_h_cv, z_h_cv;
        Mat intri_cv;
        for (int i = 0; i < y_h.size(); i++)
        {
            y_h_cv.push_back(Point3d(y_h[i](0), y_h[i](1), y_h[i](2)));
            z_h_cv.push_back(Point3d(z_h[i](0), z_h[i](1), z_h[i](2)));
        }
        eigen2cv(intrinsic, intri_cv);
        Mat es_cv = findEssentialMat(y_h_cv, z_h_cv, intri_cv, RANSAC, 0.999, 1.0, noArray());

        /* use SVD to get U and V (on SO3)*/
        Matrix3d essential, U, E0, V; // for SVD of E manifold
        cv2eigen(es_cv, essential);
        auto svd_res = essential.jacobiSvd();
        U = svd_res.matrixU();
        V = svd_res.matrixV();
        E0 = svd_res.singularValues().asDiagonal();
        if (U.determinant() < 0)
            U.col(2) = -U.col(2);
        if (V.determinant() < 0)
            V.col(2) = -V.col(2);

        MatrixXd bar_M(m, 9);
        for (int i = 0; i < m; i++)
            bar_M.row(i) = (z_h[i] * y_h[i].transpose()).reshaped(1, 9); // todo check
        Matrix<double, 9, 9> m_M = bar_M.transpose() * bar_M / m;

        Matrix3d Qx, Qy, Qz; // generator of so3
        Qx << 0, 0, 0, 0, 0, -1, 0, 1, 0;
        Qy << 0, 0, 1, 0, 0, 0, -1, 0, 0;
        Qz << 0, -1, 0, 1, 0, 0, 0, 0, 0;
        Qx /= sqrt(2);
        Qy /= sqrt(2);
        Qz /= 2;
        Matrix<double, 9, 5> Q1, Q2; // todo check reshape
        Q1 << Qx.reshaped(1, 9), Qy.reshaped(1, 9), Qz.reshaped(1, 9), Matrix<double, 9, 2>::Zero();
        Q2 << Matrix<double, 9, 2>::Zero(), Qz.reshaped(1, 9), -Qx.reshaped(1, 9), Qy.reshaped(1, 9);

        // assemble the jacobian matrix then calc the GN step
        MatrixXd J = kron(V, U).eval() *
                     (kron(E0, Matrix3d::Identity()) - kron(Matrix3d::Identity(), E0)) * Q2;
        MatrixXd first_deri = J.transpose() * m_M * essential.reshaped(1, 9);
        MatrixXd second_deri = J.transpose() * m_M * J;
        MatrixXd x = -second_deri.inverse() * first_deri; // delta X opt

        // lie algebra to lie group
        Matrix3d omega1, omega2;
        omega1 << 0, -x(2) / sqrt(2), x(1), x(2) / sqrt(2), 0, -x(0), -x(1), x(0), 0;
        omega2 << 0, x(2) / sqrt(2), x(4), -x(2) / sqrt(2), 0, -x(3), -x(4), x(3), 0;
        omega1 /= sqrt(2);
        omega2 /= sqrt(2);

        // calc exp map
        Matrix3d U_GN = U * omega1.exp().eval();
        Matrix3d V_GN = V * omega2.exp().eval();
        E_GN = U_GN * E0 * V_GN.transpose();

        // recover pose from essential mat
        Mat E_GN_cv, t_GN_cv, R_GN_cv;
        eigen2cv(E_GN, E_GN_cv);
        recoverPose(E_GN_cv, y_h_cv, z_h_cv, intri_cv, R_GN_cv, t_GN_cv);
        cv2eigen(R_GN_cv, R_GN);
        cv2eigen(t_GN_cv, t_GN);
        R = R_GN;
        t = t_GN;
    };
};

#endif // !