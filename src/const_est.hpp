#ifndef CONST_EST_HPP
#define CONST_EST_HPP

#include "eigen3/Eigen/Dense"
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

class ConsistentEst
{
#define kron kroneckerProduct
public:
    // to solve
    Matrix3d R_svd, R_GN;
    Vector3d t_svd, t_GN;
    int m; // number of correspondences

    Matrix3d Ess_svd, Ess_GN;

    // init param
    Matrix3d intrinsic;
    Matrix<double, 2, 3> W; // projection matrix
    Vector3d e1, e2, e3;    // basis of R^3

private:
    void BiasElimination(vector<Vector3d> &y_h, vector<Vector3d> &z_h)
    {
        /* --------------------get statistics-------------------- */
        MatrixXd A(m, 9);
        Matrix3d Y_bar;
        for (int i = 0; i < m; i++)
        {
            A.block<1, 9>(i, 0) = kron(y_h[i].transpose(), z_h[i].transpose()).eval();
            Y_bar += y_h[i] * y_h[i].transpose();
        }
        Y_bar /= m;
        Matrix<double, 9, 9> Q_m = A.transpose() * A / m;
        Matrix3d W_h;
        W_h << W.row(0), W.row(1), 0, 0, 0;
        Matrix<double, 9, 9> S_m = kron(Y_bar, W_h).eval();

        EigenSolver<Matrix<double, 9, 9>> es(Q_m.inverse() * S_m);
        double var_est = 1 / es.eigenvalues().real().maxCoeff();

        /* --------------------svd to get essential-------------------- */
        Matrix<double, 9, 9> Q_BE = Q_m - var_est * S_m;
        EigenSolver<Matrix<double, 9, 9>> es_svd(Q_BE);
        int idx; // get the idx of the max eigenvalue
        Matrix<double, 9, 1> eigens = es_svd.eigenvalues().real();
        eigens.minCoeff(&idx);
        // @todo check the order of the reshape mat, trans or not
        Ess_svd = es_svd.eigenvectors().real().col(idx).reshaped(3, 3);

        /* --------------------recover pose using OpenCV-------------------- */
        cv::Mat essential, intrinsic_cv, R_cv, t_cv; // for conversion
        eigen2cv(Ess_svd, essential);
        eigen2cv(intrinsic, intrinsic_cv);
        // transform homo coord to cv::Point2d
        vector<Point2d> y_h_cv, z_h_cv;
        for (auto i = 0; i < m; i++)
        {
            y_h_cv.push_back(Point2d(y_h[i](0), y_h[i](1)));
            z_h_cv.push_back(Point2d(z_h[i](0), z_h[i](1)));
        }
        // @todo may be undistortPoints are needed, or the dataset is already undistorted?
        recoverPose(essential, y_h_cv, z_h_cv, intrinsic_cv, R_cv, t_cv);
        cv2eigen(R_cv, R_svd);
        cv2eigen(t_cv, t_svd); // @todo to determine sign of t
    };

    void ManifoldGN(vector<Vector3d> &y_h, vector<Vector3d> &z_h)
    {
        /* -------------------- calc param over s-2 Manifold --------------------*/
        double alpha0 = asin(t_svd(2)), beta0;
        if (t_svd(0) > 0)
            beta0 = atan(t_svd(1) / t_svd(0));
        else
            beta0 = atan(t_svd(1) / t_svd(0)) + CV_PI;
        Vector3d par_t_beta(-cos(alpha0) * sin(beta0), cos(alpha0) * cos(beta0), 0);
        Vector3d par_t_alpha(-sin(alpha0) * cos(beta0), -sin(alpha0) * sin(beta0), cos(alpha0));

        // --------------------jacobian and pred using R_svd and t_svd--------------------
        MatrixXd J(2 * m, 5), z_pred(2 * m, 1);
        Matrix<double, 9, 3> phi;
        phi << 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0, 0;
        Matrix<double, 9, 3> t_tri = kron(Matrix3d::Identity(), t_svd).eval(); // 9 by 3

        // --------------------assemble the jacobi matrix--------------------
        for (int i = 0; i < m; i++)
        {
            Vector3d Ryih = R_svd * y_h[i];
            Matrix<double, 3, 9> C1, C2;
            C1 << -e3.transpose(), 0, 0, 0, z_h[i](0) * e3.transpose(),
                0, 0, 0, -e3.transpose(), z_h[i](1) * e3.transpose(),
                e1.transpose(), e2.transpose(), -(z_h[i](0) * e1.transpose() + z_h[i](1) * e2.transpose());
            C2 << e3.transpose(), 0, 0, 0, z_h[i](0) * e3.transpose() - e1.transpose(),
                0, 0, 0, e3.transpose(), z_h[i](1) * e3.transpose() - e2.transpose(),
                -z_h[i](0) * e3.transpose(), -z_h[i](1) * e3.transpose(), 0, 0, 0;

            Matrix3d pxi = kron(y_h[i].transpose(), R_svd).eval() * phi; // 3 by 3
            Matrix<double, 9, 3> Ryih_I = kron(Ryih, Matrix3d::Identity()).eval();
            Vector<double, 9> Ryih_t = kron(Ryih, t_svd).eval();

            // 1x1
            double den = (t_svd.transpose() * C2 * t_tri * Ryih).value();
            double num = (Ryih.transpose() * C1 * t_tri * Ryih).value();

            Vector3d par_ki_s = (den * (Ryih.transpose() * (C1 * t_tri + t_tri.transpose() * C1.transpose()) * pxi) - num * (t_svd.transpose() * C2 * t_tri * pxi)) / pow(den, 2);                              // 1 by3
            double par_ki_alpha = ((den * (Ryih.transpose() * C1 * Ryih_I * par_t_alpha) - num * (t_svd.transpose() * C2 * Ryih_I + Ryih_t.transpose() * C2.transpose()) * par_t_alpha) / pow(den, 2)).value(); // scalar
            double par_ki_beta = ((den * (Ryih.transpose() * C1 * Ryih_I * par_t_beta) - num * (t_svd.transpose() * C2 * Ryih_I + Ryih_t.transpose() * C2.transpose()) * par_t_beta) / pow(den, 2)).value();    // scalar

            double scale = num / den;
            Vector2d g = W * (R_svd * y_h[i] + scale * t_svd);
            double h = (e3.transpose() * (R_svd * y_h[i] + scale * t_svd)).value();
            Vector2d z_pred_i = g / h;
            z_pred(2 * i, 0) = z_pred_i(0); // assemble the predicted point in the second img
            z_pred(2 * i + 1, 0) = z_pred_i(1);

            J.block<2, 3>(2 * i, 0) = (h * W * (pxi + t_svd * par_ki_s.transpose()) - g * e3.transpose() * (pxi + t_svd * par_ki_s.transpose())) / (h * h);
            J.block<2, 1>(2 * i, 3) = (h * W * (scale * par_t_beta + par_ki_beta * t_svd) - g * e3.transpose() * (scale * par_t_beta + par_ki_beta * t_svd)) / (h * h); // @todo dimension wrong
            J.block<2, 1>(2 * i, 4) = (h * W * (scale * par_t_alpha + par_ki_alpha * t_svd) - g * e3.transpose() * (scale * par_t_alpha + par_ki_alpha * t_svd)) / (h * h);
        }
        Vector<double, 5> init(0, 0, 0, beta0, alpha0);
        // assemble the measurement
        VectorXd z(2 * m, 1);
        for (int i = 0; i < m; i++)
            z.block<2, 1>(2 * i, 0) << z_h[i](0), z_h[i](1);
        Vector<double, 5> res = init + (J.transpose() * J).inverse() * J.transpose() * (z - z_pred);
        Matrix3d lie;
        lie << 0, -res(2), res(1), res(2), 0, -res(0), -res(1), res(0), 0;
        R_GN = R_svd * lie.exp();
        t_GN = {cos(res(4)) * cos(res(3)), cos(res(4)) * sin(res(3)), sin(res(4))};
        t_GN.normalize();

        // for debug
        Matrix3d skew_t;
        skew_t << 0, -t_GN(2), t_GN(1), t_GN(2), 0, -t_GN(0), -t_GN(1), t_GN(0), 0;
        Ess_GN = skew_t * R_GN;
    };

public:
    ConsistentEst(Matrix3d intrinsic_) : intrinsic(intrinsic_)
    {
        W << 1, 0, 0, 0, 1, 0;
        e1 << 1, 0, 0;
        e2 << 0, 1, 0;
        e3 << 0, 0, 1;
    };

    /**
     * @brief public method for pose estimation
     *
     * @param R
     * @param t
     * @param y_h
     * @param z_h
     */
    void GetPose(Matrix3d &R, Vector3d &t, vector<Vector3d> &y_h, vector<Vector3d> &z_h)
    {
        m = y_h.size();

        BiasElimination(y_h, z_h);

        ManifoldGN(y_h, z_h);

        R = R_GN;
        t = t_GN;
    };
};

#endif // !CONSISTENT_EST_H