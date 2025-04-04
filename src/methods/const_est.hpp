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
    Matrix3d R_svd, R_GN;
    Vector3d t_svd, t_GN;
    int m;
    double var_est;

    Matrix3d Ess_svd, Ess_GN;
    Matrix3d intrinsic;
    Matrix<double, 2, 3> W;
    Vector3d e1, e2, e3;

private:
    void BiasElimination(vector<Vector3d> &y_h, vector<Vector3d> &z_h, vector<Point2d> y_p_cv, vector<Point2d> z_p_cv)
    {
        MatrixXd A(m, 9);
        Matrix3d Y_bar = Matrix3d::Zero();
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
        var_est = 1 / es.eigenvalues().real().maxCoeff();

        Matrix<double, 9, 9> Q_BE = Q_m - var_est * S_m;
        EigenSolver<Matrix<double, 9, 9>> es_svd(Q_BE);
        int idx;
        Matrix<double, 9, 1> eigens = es_svd.eigenvalues().real();
        eigens.minCoeff(&idx);
        Ess_svd = es_svd.eigenvectors().real().col(idx).reshaped(3, 3);

        cv::Mat essential, intrinsic_cv, R_cv, t_cv;
        eigen2cv(Ess_svd, essential);
        eigen2cv(intrinsic, intrinsic_cv);
        recoverPose(essential, y_p_cv, z_p_cv, intrinsic_cv, R_cv, t_cv);
        cv2eigen(R_cv, R_svd);
        cv2eigen(t_cv, t_svd);
    }

    void ManifoldGN(vector<Vector3d> &y_h, vector<Vector3d> &z_h, int max_iterations)
    {
        Matrix3d R_curr = R_svd;
        Vector3d t_curr = t_svd;
        double alpha_curr = asin(t_curr(2));
        double beta_curr = (t_curr(0) > 0) ? atan(t_curr(1) / t_curr(0)) 
                                          : atan(t_curr(1) / t_curr(0)) + CV_PI;

        for (int iter = 0; iter < max_iterations; iter++)
        {
            Vector3d par_t_beta(-cos(alpha_curr) * sin(beta_curr),
                              cos(alpha_curr) * cos(beta_curr), 0);
            Vector3d par_t_alpha(-sin(alpha_curr) * cos(beta_curr),
                               -sin(alpha_curr) * sin(beta_curr),
                               cos(alpha_curr));

            MatrixXd J(2 * m, 5), z_pred(2 * m, 1);
            Matrix<double, 9, 3> phi;
            phi << 0, 0, 0,
                   0, 0, 1,
                   0, -1, 0,
                   0, 0, -1,
                   0, 0, 0,
                   1, 0, 0,
                   0, 1, 0,
                   -1, 0, 0,
                   0, 0, 0;

            Matrix<double, 9, 3> t_tri = kron(Matrix3d::Identity(), t_curr).eval();

            for (int i = 0; i < m; i++)
            {
                Vector3d Ryih = R_curr * y_h[i];
                Matrix<double, 3, 9> C1, C2;
                C1 << -e3.transpose(), 0, 0, 0, z_h[i](0) * e3.transpose(),
                      0, 0, 0, -e3.transpose(), z_h[i](1) * e3.transpose(),
                      e1.transpose(), e2.transpose(), -(z_h[i](0) * e1.transpose() + z_h[i](1) * e2.transpose());
                C2 << e3.transpose(), 0, 0, 0, z_h[i](0) * e3.transpose() - e1.transpose(),
                      0, 0, 0, e3.transpose(), z_h[i](1) * e3.transpose() - e2.transpose(),
                      -z_h[i](0) * e3.transpose(), -z_h[i](1) * e3.transpose(), 0, 0, 0;

                Matrix3d pxi = kron(y_h[i].transpose(), R_curr).eval() * phi;
                Matrix<double, 9, 3> Ryih_I = kron(Ryih, Matrix3d::Identity()).eval();
                Vector<double, 9> Ryih_t = kron(Ryih, t_curr).eval();

                double den = (t_curr.transpose() * C2 * t_tri * Ryih).value();
                double num = (Ryih.transpose() * C1 * t_tri * Ryih).value();

                Vector3d par_ki_s = (den * (Ryih.transpose() * (C1 * t_tri + t_tri.transpose() * C1.transpose()) * pxi) -
                                   num * (t_curr.transpose() * C2 * t_tri * pxi)) / pow(den, 2);
                double par_ki_alpha = ((den * (Ryih.transpose() * C1 * Ryih_I * par_t_alpha) -
                                    num * (t_curr.transpose() * C2 * Ryih_I + Ryih_t.transpose() * C2.transpose()) * par_t_alpha) /
                                    pow(den, 2)).value();
                double par_ki_beta = ((den * (Ryih.transpose() * C1 * Ryih_I * par_t_beta) -
                                   num * (t_curr.transpose() * C2 * Ryih_I + Ryih_t.transpose() * C2.transpose()) * par_t_beta) /
                                   pow(den, 2)).value();

                double scale = num / den;
                Vector2d g = W * (R_curr * y_h[i] + scale * t_curr);
                double h = (e3.transpose() * (R_curr * y_h[i] + scale * t_curr)).value();
                Vector2d z_pred_i = g / h;
                z_pred(2 * i, 0) = z_pred_i(0);
                z_pred(2 * i + 1, 0) = z_pred_i(1);

                J.block<2, 3>(2 * i, 0) = (h * W * (pxi + t_curr * par_ki_s.transpose()) -
                                         g * e3.transpose() * (pxi + t_curr * par_ki_s.transpose())) / (h * h);
                J.block<2, 1>(2 * i, 3) = (h * W * (scale * par_t_beta + par_ki_beta * t_curr) -
                                         g * e3.transpose() * (scale * par_t_beta + par_ki_beta * t_curr)) / (h * h);
                J.block<2, 1>(2 * i, 4) = (h * W * (scale * par_t_alpha + par_ki_alpha * t_curr) -
                                         g * e3.transpose() * (scale * par_t_alpha + par_ki_alpha * t_curr)) / (h * h);
            }

            Vector<double, 5> init = Vector<double, 5>::Zero();
            VectorXd z(2 * m, 1);
            for (int i = 0; i < m; i++)
                z.block<2, 1>(2 * i, 0) << z_h[i](0), z_h[i](1);

            Vector<double, 5> delta = (J.transpose() * J).inverse() * J.transpose() * (z - z_pred);

            Matrix3d lie;
            lie << 0, -delta(2), delta(1),
                   delta(2), 0, -delta(0),
                   -delta(1), delta(0), 0;

            R_curr = R_curr * lie.exp().eval();
            alpha_curr += delta(4);
            beta_curr += delta(3);
            t_curr = Vector3d(cos(alpha_curr) * cos(beta_curr),
                            cos(alpha_curr) * sin(beta_curr),
                            sin(alpha_curr));
            t_curr.normalize();
        }

        R_GN = R_curr;
        t_GN = t_curr;

        Matrix3d skew_t;
        skew_t << 0, -t_GN(2), t_GN(1),
                  t_GN(2), 0, -t_GN(0),
                  -t_GN(1), t_GN(0), 0;
        Ess_GN = skew_t * R_GN;
    }

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
    void GetPose(Matrix3d &R, Vector3d &t, vector<Vector3d> &y_h, vector<Vector3d> &z_h, vector<Point2d> &ypix, vector<Point2d> &zpix, int iter=1)
    {
        m = y_h.size();

        BiasElimination(y_h, z_h, ypix, zpix);

        ManifoldGN(y_h, z_h,iter);

        R = R_GN;
        t = t_GN;
    };
};

#endif // !CONSISTENT_EST_H