// Copyright 2017 Thomas Schöps
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

#include <cstring>
#include <iostream>
#include <cmath>
#include <string>
#include <chrono>
#include <cstdlib>
#include <ctime>

#include "cameras.h"
#include "images.h"
#include "pt3d.hpp"
#include "const_est.hpp"
#include "manifold_GN.hpp"
#include "npt_pose.h"
#include "utils.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        std::cout << "Usage: " << argv[0] << "  <datasetname>  <img_windows_size> <pts>" << std::endl;
        std::cout << "you should put your dataset in the <dataset> dir of this proj";
        return EXIT_FAILURE; // <path_to_pt3ds_txt>
    }
    string dtset(argv[1]), windows_size_str(argv[2]), pts_str(argv[3]);
    int wds = stoi(windows_size_str);
    int pts = stoi(pts_str);
    string cameras_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/cameras.txt";
    string images_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/images.txt";

    // Load cameras (indexed by: camera_id).
    ColmapCameraPtrMap cameras;
    bool success = ReadColmapCameras(cameras_txt_path, &cameras);
    if (success)
        std::cout << "Successfully loaded " << cameras.size() << " camera(s)." << std::endl;
    else
    {
        std::cout << "Error: could not load cameras.Exit" << std::endl;
        return EXIT_FAILURE;
    }

    // Load images (indexed by: image_id).
    ColmapImagePtrMap images;
    success = ReadColmapImages(images_txt_path, true, &images, cameras);
    if (success)
        std::cout << "Successfully loaded " << images.size() << " image info(s)." << std::endl;
    else
    {
        std::cout << "Error: could not load image info.Exit" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "[" << dtset << "]eval on going, please wait..." << endl;

    // accumulate the error

    eval c_est("c_est"), e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"), pt7ransac("pt7ransac");
    int valid_round = 0; // the number of valid pairs, only when the number of covisible points is larger than 20, the pair is valid
    int progress = images.size() * wds - (wds * wds + wds / 2), tmp_finish = 0;

    vector<Vector3d> yptmp, zptmp, yntmp, zntmp;
    vector<Vector3d> y_pix, z_pix, y_n, z_n;
    vector<Point2d> y_cv_pix, z_cv_pix;
    vector<Point2d> ycvtmp, zcvtmp;
    yptmp.reserve(4000);
    zptmp.reserve(4000);
    yntmp.reserve(4000);
    zntmp.reserve(4000);
    y_pix.reserve(4000);
    z_pix.reserve(4000);
    y_n.reserve(4000);
    z_n.reserve(4000);
    y_cv_pix.reserve(4000);
    z_cv_pix.reserve(4000);
    ycvtmp.reserve(4000);
    zcvtmp.reserve(4000);

    vector<int> num_pts = {10, 20, 40, 80, 160, 320, 640, 1280, 2560, 3000, -1};
    

    for (int img1 = 1; img1 < progress; ++img1)
    {
        for (int img2 = img1 + 1; img2 <= img1 + wds - 1 && img2 <= images.size(); ++img2, ++tmp_finish)
        {
            printProg(tmp_finish, progress);

            // clear all vector
            yptmp.clear();
            zptmp.clear();
            yntmp.clear();
            zntmp.clear();
            y_pix.clear();
            z_pix.clear();
            y_n.clear();
            z_n.clear();
            y_cv_pix.clear();
            z_cv_pix.clear();
            ycvtmp.clear();
            zcvtmp.clear();

            if (img1 == img2 || images[img1]->camera_id != images[img2]->camera_id) // only consider imgs with the same camera for convenience
                continue;

            for (auto yi : images[img1]->observations) // traverse obs in first img
            {
                if (images[img2]->observations.count(yi.first)) // if the point is visible in the another image
                {
                    point_pair p1 = yi.second, p2 = images[img2]->observations[yi.first];
                    yptmp.push_back(p1.pixel);
                    zptmp.push_back(p2.pixel);
                    ycvtmp.push_back(Point2d(p1.pixel.x(), p1.pixel.y()));
                    zcvtmp.push_back(Point2d(p2.pixel.x(), p2.pixel.y()));
                    yntmp.push_back(p1.normalized);
                    zntmp.push_back(p2.normalized);
                }
            }
            int covisible_num = yptmp.size();
            if (covisible_num < 3000)
                continue; // for progressive evaluation
            ++valid_round;

            // compute the ground truth relative pose
            Matrix3d R1 = images[img1]->global_T_image.linear().transpose().cast<double>();
            Matrix3d R_gt = ((images[img1]->global_T_image.linear().transpose() * images[img1 + 1]->global_T_image.linear()).cast<double>()).transpose();
            Vector3d t_gt = -R1 * (images[img1 + 1]->global_T_image.translation() - images[img1]->global_T_image.translation()).cast<double>();
            Vector3d t_with_scale = t_gt;
            t_gt.normalize(); // only the bearing vector is needed
            Matrix3d t_skew;
            t_skew << 0, -t_gt(2), t_gt(1),
                t_gt(2), 0, -t_gt(0),
                -t_gt(1), t_gt(0), 0;
            Matrix3d E_ground = t_skew * R_gt;
            Matrix3d R_log = R_gt.log().eval();
            Vector3d R_lie(R_log(2, 1), R_log(0, 2), R_log(1, 0));

            // temp variables to store res
            Matrix3d R_estimated;
            Vector3d t_estimated;

            double r_err_this_round[5] = {0}, t_err_this_round[5] = {0};
            double time_elapse[5] = {0};
            double est_vars = 0;

            for (auto now_pts : num_pts)
            {
                for (auto iter = 1; iter < 11; iter++)
                {
                    // randomly choose num points,sample from the covisible points
                    if (pts == -1)
                    {
                        y_pix = yptmp;
                        z_pix = zptmp;
                        y_n = yntmp;
                        z_n = zntmp;
                        y_cv_pix = ycvtmp;
                        z_cv_pix = zcvtmp;
                        iter += 10; // only run once
                    }
                    else
                    {
                        // change rand seed
                        std::srand(std::time(nullptr));
                        y_pix.clear();
                        z_pix.clear();
                        y_n.clear();
                        z_n.clear();
                        y_cv_pix.clear();
                        z_cv_pix.clear();
                        for (int i = 0; i < pts; ++i)
                        {
                            int idx = std::rand() % covisible_num;
                            y_pix.push_back(yptmp[idx]);
                            z_pix.push_back(zptmp[idx]);
                            y_n.push_back(yntmp[idx]);
                            z_n.push_back(zntmp[idx]);
                            y_cv_pix.push_back(ycvtmp[idx]);
                            z_cv_pix.push_back(zcvtmp[idx]);
                        }
                    }
                    /* ↓------------------consistent estimator------------------↓ */
                    ConsistentEst est(cameras[images[img1]->camera_id]->intrinsic);
                    time_elapse[0] += TIME_IT(est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix););
                    r_err_this_round[0] += (R_estimated - R_gt).norm();
                    t_err_this_round[0] += abs(t_estimated.dot(t_gt));
                    est_vars += est.var_est;

                    /* ↑------------------consistent estimator------------------↑ */

                    /* ↓------------------E Manifold GN------------------↓ */
                    ManifoldGN MGN(cameras[images[img1]->camera_id]->intrinsic);
                    time_elapse[1] += TIME_IT(MGN.GetPose(R_estimated, t_estimated, y_cv_pix, z_cv_pix, y_n, z_n););
                    r_err_this_round[1] += (R_estimated - R_gt).norm();
                    t_err_this_round[1] += abs(t_estimated.dot(t_gt));

                    /* ↑------------------E Manifold GN------------------↑ */

                    Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars

                    /* ↓------------------RANSAC method------------------↓ */
                    time_elapse[2] += TIME_IT(eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intrinsic_cv);
                                              E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
                                              recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                              cv2eigen(R_cv, R_estimated);
                                              cv2eigen(t_cv, t_estimated););
                    r_err_this_round[2] += (R_estimated - R_gt).norm();
                    t_err_this_round[2] += abs(t_estimated.dot(t_gt));

                    time_elapse[3] += TIME_IT(eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intrinsic_cv);
                                              E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::LMEDS, 0.999, 1.0);
                                              recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                              cv2eigen(R_cv, R_estimated);
                                              cv2eigen(t_cv, t_estimated););
                    r_err_this_round[3] += (R_estimated - R_gt).norm();
                    t_err_this_round[3] += abs(t_estimated.dot(t_gt));

                    /* ↑------------------RANSAC method------------------↑ */

                    /* ↓------------------SDP on Essential Mat------------------↓ */
                    // prepare data
                    double *npt_p1 = new double[3 * y_n.size()], *npt_p2 = new double[3 * y_n.size()];
                    for (size_t i = 0; i < y_n.size(); ++i)
                    {
                        npt_p1[3 * i] = y_n[i](0);
                        npt_p1[3 * i + 1] = y_n[i](1);
                        npt_p1[3 * i + 2] = y_n[i](2);
                        npt_p2[3 * i] = z_n[i](0);
                        npt_p2[3 * i + 1] = z_n[i](1);
                        npt_p2[3 * i + 2] = z_n[i](2);
                    }
                    // solve
                    double *C = new double[81];
                    Eigen::Matrix<double, 12, 12> X_sol;
                    Matrix3d sdp_E;
                    time_elapse[4] += TIME_IT(npt_pose(npt_p1, npt_p2, C, y_n.size(), X_sol, sdp_E, true);
                                              eigen2cv(sdp_E, E_cv);
                                              cv::recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                              cv2eigen(R_cv.t(), R_estimated);
                                              cv2eigen(t_cv, t_estimated););
                    r_err_this_round[4] += (R_estimated - R_gt).norm();
                    t_err_this_round[4] += abs(t_estimated.dot(t_gt));

                    /* ↑------------------SDP on Essential Mat------------------↑ */
                }
            }

            if (pts != -1)
            {
                for (int i = 0; i < 5; ++i)
                {
                    r_err_this_round[i] /= 10;
                    t_err_this_round[i] /= 10;
                    time_elapse[i] /= 10;
                }
                est_vars /= 10;
            }
        }
    }

    return EXIT_SUCCESS;
}
