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

public:
    Matrix3d intrinsic;
    Matrix3d R_GN, R_init;
    Vector3d t_GN, t_init;
    Matrix3d E_init, E_GN;
    int m;

public:
    ManifoldGN(Matrix3d intrinsic_) : intrinsic(intrinsic_){};

    void GetPose(Matrix3d &R, Vector3d &t, vector<Point2d> &y_h_cv, vector<Point2d> &z_h_cv,vector<Vector3d> ynn,vector<Vector3d>znn,Matrix3d E_out = Matrix3d::Zero())
    {
        E_init = E_out;
        // convert vector3d to point3d type for OpenCV
        m = y_h_cv.size();
        Mat intri_cv;
        eigen2cv(intrinsic, intri_cv);
        
        Mat es_cv, R_cv, t_cv;
        if(E_init == Matrix3d::Zero()) // given a initial E or not
            es_cv = findEssentialMat(y_h_cv, z_h_cv, intri_cv, RANSAC, 0.999, 1.0);
        else
            eigen2cv(E_init,es_cv);

        recoverPose(es_cv, y_h_cv, z_h_cv, intri_cv, R_cv, t_cv);
        cv2eigen(R_cv, R_init);
        cv2eigen(t_cv, t_init);

        /* use SVD to get U and V (on SO3)*/
        Matrix3d U, E0, V; // for SVD of E manifold
        cv2eigen(es_cv, E_init);
        auto svd_res = E_init.jacobiSvd(ComputeFullU | ComputeFullV);
        U = svd_res.matrixU();
        V = svd_res.matrixV();
        E0 = svd_res.singularValues().asDiagonal();
        if (U.determinant() == -1)
            U.col(2) = -U.col(2);
        if (V.determinant() == -1)
            V.col(2) = -V.col(2);

        MatrixXd bar_M(m, 9);
        for (int i = 0; i < m; i++)
            bar_M.row(i) = (znn[i] * ynn[i].transpose()).reshaped(1, 9);
        Matrix<double, 9, 9> m_M = bar_M.transpose() * bar_M / m;

        Vector<double, 9> Qx, Qy, Qz; // generator of so3
        Qx << 0, 0, 0, 0, 0, 1, 0, -1, 0;
        Qy << 0, 0, -1, 0, 0, 0, 1, 0, 0;
        Qz << 0, 1, 0, -1, 0, 0, 0, 0, 0;
        Qx /= sqrt(2);
        Qy /= sqrt(2);
        Qz /= 2;
        Matrix<double, 9, 5> Q1, Q2; // todo check reshape
        Matrix<double, 9, 2> zero_stand = Matrix<double, 9, 2>::Zero();
        Q1 << Qx, Qy, Qz, zero_stand;
        Q2 << zero_stand, -Qz, Qx, Qy;

        // assemble the jacobian matrix then calc the GN step
        MatrixXd J = kron(V, U).eval() *
                     (kron(E0, Matrix3d::Identity()).eval() * Q1 - kron(Matrix3d::Identity(), E0).eval() * Q2);
        MatrixXd first_deri = J.transpose() * m_M * E_init.reshaped(9, 1);
        MatrixXd second_deri = J.transpose() * m_M * J;
        MatrixXd x = -second_deri.inverse() * first_deri; // delta X opt

        // lie algebra to lie group
        Matrix3d omega1, omega2;
        omega1 << 0, -x(2) / sqrt(2), x(1),  x(2) / sqrt(2), 0,  -x(0),   -x(1), x(0), 0;
        omega2 << 0, x(2) / sqrt(2), x(4),   -x(2) / sqrt(2), 0, -x(3),   -x(4), x(3), 0;
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