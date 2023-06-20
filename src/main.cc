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

#include "cameras.h"
#include "images.h"
#include "pt3d.hpp"
#include "const_est.hpp"
#include "manifold_GN.hpp"
#include "utils.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		std::cout << "Usage: " << argv[0] << "datasetname" << std::endl;
		std::cout << "you should put your dataset in the <dataset> dir of this proj";
		return EXIT_FAILURE; // <path_to_pt3ds_txt>
	}
	string dtset(argv[1]);
	string cameras_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/cameras.txt";
	string images_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/images.txt";
	// string points3D_txt_path = " ../" + dtset + "/dslr_calibration_undistorted/points3D.txt";

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
	eval c_est("c_est"),
		e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"), pt7ransac("pt7ransac");
	int valid_round = 0; // the number of valid pairs, only when the number of covisible points is larger than 20, the pair is valid
	int progress = images.size() * 7 - 28, tmp_finish = 0;

	for (int img1 = 1; img1 < progress; ++img1)
	{
		// find the covisible points
		// surfix _pix means homogeneous coord in pixel plane, _n means normalized coordinates
		vector<Vector3d> y_pix, z_pix, y_n, z_n;
		vector<Point2d> y_cv_pix, z_cv_pix;
		for (int img2 = img1 + 1; img2 <= img1 + 7 && img2 <= images.size(); ++img2, ++tmp_finish)
		{
			printProg(tmp_finish, progress);
			if (img1 == img2 || images[img1]->camera_id != images[img2]->camera_id) // only consider imgs with the same camera for convenience
				continue;

			/**
			 * @brief imgs in this dataset have covisible points only when they are close enough
			 * 	      so we just check closest 6 imgs
			 */
			for (auto yi : images[img1]->observations) // traverse obs in first img
			{
				if (images[img2]->observations.count(yi.first)) // if the point is visible in the another image
				{
					point_pair p1 = yi.second, p2 = images[img2]->observations[yi.first];
					y_pix.push_back(p1.pixel);
					z_pix.push_back(p2.pixel);
					y_cv_pix.push_back(Point2d(p1.pixel.x(), p1.pixel.y()));
					z_cv_pix.push_back(Point2d(p2.pixel.x(), p2.pixel.y()));
					y_n.push_back(p1.normalized);
					z_n.push_back(p2.normalized);
				}
			}
			if (y_pix.size() < 20)
			{
				std::cout << "[point too few,only " << y_pix.size() << " available,skip]" << endl;
				continue;
			}
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

			// temp variables to store res
			Matrix3d R_estimated;
			Vector3d t_estimated;
			double r_err_this_round, t_err_this_round;

			/* ---------------------------remove outlier if necessary--------------------------------------*/
			vector<Vector3d> y_pix_in, z_pix_in, y_n_in, z_n_in; // inliers
			cv::Mat inlier_mask, intri_cv;
			cv::eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intri_cv);
			cv::findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, cv::RANSAC, 0.999, 3.0, inlier_mask);
			for (int i = 0; i < inlier_mask.rows; ++i) // remove outlier using inlier_mask
			{
				if (inlier_mask.at<uchar>(i, 0))
				{
					y_pix_in.push_back(y_pix[i]);
					z_pix_in.push_back(z_pix[i]);
					y_n_in.push_back(y_n[i]);
					z_n_in.push_back(z_n[i]);
				}
			}
			int inlier_num = y_pix_in.size();

			//------------------------------ solve the epipolar problem here-----------------------------------

			/* ↓ ↓ ↓ ↓ ↓ ↓---consistent estimator--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
			ConsistentEst est(cameras[images[img1]->camera_id]->intrinsic);
			est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix);
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = abs(t_estimated.dot(t_gt));
			calcEval(c_est, img1, img2, t_err_this_round, r_err_this_round, inlier_num, est.var_est);
			// std::cout << valid_round << "--------" << endl;
			// std::cout << "[R_est]\n"
			// 		  << R_estimated << endl;
			// std::cout << "[GT_R]\n"
			// 		  << R_gt << endl;
			// std::cout << " [E_est]\n"
			// 		  << est.Ess_svd.normalized() << endl;
			// std::cout << "[E_GT]\n"
			// 		  << E_ground.normalized() << endl;

			// std::cout << "[T_est] " << t_estimated.transpose() << endl;
			// std::cout << "[GT_T] " << t_gt.transpose() << endl;
			// std::cout << "[scale T]" << t_with_scale.transpose() << endl;
			// std::cout << " [norm of T]" << t_with_scale.norm() << endl;
			/* ↑↑↑↑↑↑↑↑↑↑↑↑--consistent estimator--↑↑↑↑↑↑↑↑↑↑↑↑ */

			/* ↓ ↓ ↓ ↓ ↓ ↓---E Manifold GN--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
			ManifoldGN MGN(cameras[images[img1]->camera_id]->intrinsic);
			MGN.GetPose(R_estimated, t_estimated, y_pix, z_pix, y_n, z_n);
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = abs(t_estimated.dot(t_gt));
			calcEval(e_m_gn, img1, img2, t_err_this_round, r_err_this_round, inlier_num);

			/* ↑↑↑↑↑↑↑↑↑↑↑↑--E Manifold GN--↑↑↑↑↑↑↑↑↑↑↑↑ */

			/* ↓ ↓ ↓ ↓ ↓ ↓--- RANSAC method--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
			Mat E_cv, intrinsic_cv, R_cv, t_cv;
			eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intrinsic_cv);
			E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
			recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
			cv2eigen(R_cv, R_estimated);
			cv2eigen(t_cv, t_estimated);
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = abs(t_estimated.dot(t_gt));
			calcEval(pt5ransac, img1, img2, t_err_this_round, r_err_this_round, inlier_num);

			E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::LMEDS, 0.999, 1.0);
			recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
			cv2eigen(R_cv, R_estimated);
			cv2eigen(t_cv, t_estimated);
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = abs(t_estimated.dot(t_gt));
			calcEval(pt7ransac, img1, img2, t_err_this_round, r_err_this_round, inlier_num);
			/* ↑↑↑↑↑↑↑↑↑↑↑↑-- RANSAC method--↑↑↑↑↑↑↑↑↑↑↑↑ */
			// std::cout << "round: " << valid_round << endl;
		}
	}

	// ------------------------------------ save results here----------------------------------------
	std::cout << "[c_est] R:" << c_est.total_R_Fn << " t: " << c_est.total_t_cos << endl;
	// std::cout << "[sdp] R:" << c_est.total_R_Fn << " t: " << c_est.total_t_cos << endl;
	std::cout << "[e_m_gn] R:" << e_m_gn.total_R_Fn << " t: " << e_m_gn.total_t_cos << endl;
	std::cout << "[pt5ransac] R:" << pt5ransac.total_R_Fn << " t: " << pt5ransac.total_t_cos << endl;
	std::cout << "[pt7ransac] R:" << pt7ransac.total_R_Fn << " t: " << pt7ransac.total_t_cos << endl;
	std::cout << "--------------------------------------" << endl;

	saveRes(c_est, dtset);
	// saveRes(sdp, dtset);
	saveRes(e_m_gn, dtset);
	saveRes(pt5ransac, dtset);
	saveRes(pt7ransac, dtset);

	// Load 3d points ,index by point3d_id
	// ColmapPoint3DPtrMap point3ds;
	// success = ReadColmapPt3d(points3d_txt_path, point3ds);
	// if (success)
	// 	std::cout << "Successfully loaded " << point3ds.size() << " point3d(s)." << std::endl;
	// else
	// 	std::cout << "Error: could not load point3ds." << std::endl;

	return EXIT_SUCCESS;
}
