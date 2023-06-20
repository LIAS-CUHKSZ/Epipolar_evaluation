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
	if (argc != 3)
	{
		std::cout << "Usage: " << argv[0] << " <path_to_cameras_txt> <path_to_images_txt> " << std::endl;
		return EXIT_FAILURE; // <path_to_pt3ds_txt>
	}

	const char *cameras_txt_path = argv[1];
	const char *images_txt_path = argv[2];
	// const char *points3d_txt_path = argv[3];

	// Load cameras (indexed by: camera_id).
	ColmapCameraPtrMap cameras;
	bool success = ReadColmapCameras(cameras_txt_path, &cameras);
	if (success)
		std::cout << "Successfully loaded " << cameras.size() << " camera(s)." << std::endl;
	else
		std::cout << "Error: could not load cameras." << std::endl;

	// Load images (indexed by: image_id).
	ColmapImagePtrMap images;
	success = ReadColmapImages(images_txt_path, /* read_observations */ true, &images, cameras);
	if (success)
		std::cout << "Successfully loaded " << images.size() << " image info(s)." << std::endl;
	else
		std::cout << "Error: could not load image infos." << std::endl;

	// accumulate the error
	eval c_est("c_est"), e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"), pt7ransac("pt7ransac");
	int valid_round = 0;

	// for each piar of imgs in images, find covisible points and compute the relative pose
	for (int i = 1; i < images.size(); ++i)
	{
		// find the covisible points
		// surfix _h means homogeneous coordinates, _n means normalized coordinates
		vector<Vector3d> y_pix, z_pix, y_n, z_n;
		vector<Point2d> y_cv_pix, z_cv_pix;

		vector<Vector3d> y_pix_in, z_pix_in, y_n_in, z_n_in; // inliers

		for (auto yi : images[i]->observations)
		{
			if (images[i + 1]->observations.count(yi.first)) // if the point is visible in the next image
			{
				point_pair p1 = yi.second, p2 = images[i + 1]->observations[yi.first];
				y_pix.push_back(Eigen::Vector3d(p1.pixel.x(), p1.pixel.y(), 1.0));
				z_pix.push_back(Eigen::Vector3d(p2.pixel.x(), p2.pixel.y(), 1.0));
				y_n.push_back(Eigen::Vector3d(p1.normalized.x(), p1.normalized.y(), 1.0));
				z_n.push_back(Eigen::Vector3d(p2.normalized.x(), p2.normalized.y(), 1.0));
				y_cv_pix.push_back(Point2d(p1.pixel.x(), p1.pixel.y()));
				z_cv_pix.push_back(Point2d(p2.pixel.x(), p2.pixel.y()));
			}
		}
		if (y_pix.size() < 20)
		{
			std::cout << "[************+++++++point too few,only" << y_pix.size() << endl;
			continue;
		}
		++valid_round;
		std::cout << "total pt numbers:" << y_pix.size() << endl;

		// compute the ground truth relative pose
		Matrix3d R1 = images[i]->global_T_image.linear().transpose().cast<double>();
		Matrix3d R_gt = ((images[i]->global_T_image.linear().transpose() * images[i + 1]->global_T_image.linear()).cast<double>()).transpose();
		Vector3d t_gt = -R1 * (images[i + 1]->global_T_image.translation() - images[i]->global_T_image.translation()).cast<double>();
		Vector3d t_with_scale = t_gt;
		t_gt.normalize(); // only the bearing vector is needed

		// temp variables
		Matrix3d R_estimated;
		Vector3d t_estimated;
		double r_this, t_this;

		/* remove outlier if necessary */
		cv::Mat inlier_mask, intri_cv;
		cv::eigen2cv(cameras[images[i]->camera_id]->intrinsic, intri_cv);
		cv::findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, cv::RANSAC, 0.999, 3.0, inlier_mask);
		// remove outlier using inlier_mask
		for (int i = 0; i < inlier_mask.rows; ++i)
		{
			if (inlier_mask.at<uchar>(i, 0))
			{
				y_pix_in.push_back(y_pix[i]);
				z_pix_in.push_back(z_pix[i]);
				y_n_in.push_back(y_n[i]);
				z_n_in.push_back(z_n[i]);
			}
		}
		std::cout << "inlier numbers:" << y_pix_in.size() << endl;
		std::cout << "----------" << endl;

		//------------------------------ solve the epipolar problem here-----------------------------------

		/* ↓ ↓ ↓ ↓ ↓ ↓---consistent estimator--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
		ConsistentEst est(cameras[images[i]->camera_id]->intrinsic);
		est.GetPose(R_estimated, t_estimated, y_n_in, z_n_in);
		r_this = (R_estimated - R_gt).norm();
		t_this = abs(t_estimated.dot(t_gt));
		calcEval(c_est, r_this, t_this, est.var_est);
		std::cout << "[R_est]\n"
				  << R_estimated << endl;
		std::cout << "[GT_R]\n"
				  << R_gt << endl;
		std::cout << "[T_est] " << t_estimated.transpose() << endl;
		std::cout << "[GT_T] " << t_gt.transpose() << endl;
		std::cout << "[scale T]" << t_with_scale.transpose() << endl;
		std::cout << "[norm of T]" << t_with_scale.norm() << endl;
		/* ↑↑↑↑↑↑↑↑↑↑↑↑--consistent estimator--↑↑↑↑↑↑↑↑↑↑↑↑ */

		/* ↓ ↓ ↓ ↓ ↓ ↓---E Manifold GN--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
		ManifoldGN MGN(cameras[images[i]->camera_id]->intrinsic);
		MGN.GetPose(R_estimated, t_estimated, y_pix, z_pix, y_n, z_n);
		// R_estimated = MGN.R_init;
		// t_estimated = MGN.t_init;
		r_this = (R_estimated - R_gt).norm();
		t_this = abs(t_estimated.dot(t_gt));
		calcEval(e_m_gn, r_this, t_this);
		/* ↑↑↑↑↑↑↑↑↑↑↑↑--E Manifold GN--↑↑↑↑↑↑↑↑↑↑↑↑ */

		/* ↓ ↓ ↓ ↓ ↓ ↓--- RANSAC method--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
		Mat E_cv, intrinsic_cv, R_cv, t_cv;
		eigen2cv(cameras[images[i]->camera_id]->intrinsic, intrinsic_cv);
		E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
		recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
		cv2eigen(R_cv, R_estimated);
		cv2eigen(t_cv, t_estimated);
		r_this = (R_estimated - R_gt).norm();
		t_this = abs(t_estimated.dot(t_gt));
		calcEval(pt5ransac, r_this, t_this);

		E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::LMEDS, 0.999, 1.0);
		recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
		cv2eigen(R_cv, R_estimated);
		cv2eigen(t_cv, t_estimated);
		r_this = (R_estimated - R_gt).norm();
		t_this = abs(t_estimated.dot(t_gt));
		calcEval(pt7ransac, r_this, t_this);
		/* ↑↑↑↑↑↑↑↑↑↑↑↑-- RANSAC method--↑↑↑↑↑↑↑↑↑↑↑↑ */

		// ------------------------------ save results here-----------------------------------

		std::cout << "[R_est]\n"
				  << R_estimated << endl;
		std::cout << "[GT_R]\n"
				  << R_gt << endl;
		std::cout << "[T_est] " << t_estimated.transpose() << endl;
		std::cout << "[GT_T] " << t_gt.transpose() << endl;
		std::cout << "[scale T]" << t_with_scale.transpose() << endl;
		std::cout << "[norm of T]" << t_with_scale.norm() << endl;
	}
	std::cout << "[c_est] R:" << c_est.total_R_Fn << " t: " << c_est.total_t_cos << endl;
	std::cout << "[e_m_gn] R:" << e_m_gn.total_R_Fn << " t: " << e_m_gn.total_t_cos << endl;
	std::cout << "[pt5ransac] R:" << pt5ransac.total_R_Fn << " t: " << pt5ransac.total_t_cos << endl;
	std::cout << "[pt7ransac] R:" << pt7ransac.total_R_Fn << " t: " << pt7ransac.total_t_cos << endl;

	for (auto i : c_est.noise)
		std::cout << i << " " << endl;

	// Load 3d points ,index by point3d_id
	// ColmapPoint3DPtrMap point3ds;
	// success = ReadColmapPt3d(points3d_txt_path, point3ds);
	// if (success)
	// 	std::cout << "Successfully loaded " << point3ds.size() << " point3d(s)." << std::endl;
	// else
	// 	std::cout << "Error: could not load point3ds." << std::endl;

	return EXIT_SUCCESS;
}
