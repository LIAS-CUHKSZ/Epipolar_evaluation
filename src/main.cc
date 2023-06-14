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

using namespace cv;
using namespace std;
using namespace Eigen;

struct eval
{
	string method_name;
	vector<double> t_err_per_round;
	vector<double> R_err_per_round;
	double total_R_Fn;
	double total_t_cos;
};

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
	double frobenius_error = 0, total_cosine = 0;
	unordered_map<string, eval> methods_eval;

	// for each piar of imgs in images, find covisible points and compute the relative pose
	for (int i = 1; i < images.size(); ++i)
	{
		// find the covisible points
		vector<Vector3d> y_h, z_h;
		for (auto yi : images[i]->observations)
		{
			if (images[i + 1]->observations.count(yi.first)) // if the point is visible in the next image
			{
				y_h.push_back(Eigen::Vector3d(yi.second.x(), yi.second.y(), 1.0));
				z_h.push_back(Eigen::Vector3d(images[i + 1]->observations[yi.first].x(), images[i + 1]->observations[yi.first].y(), 1.0));
			}
		}
		if (y_h.size() < 20)
		{
			std::cout << "[************+++++++point too few,only" << y_h.size() << endl;
			continue;
		}
		std::cout << "total pt numbers:" << y_h.size() << endl;

		// compute the relative pose
		Matrix3d R1 = images[i]->global_T_image.linear().transpose().cast<double>();
		Matrix3d R_gt = ((images[i]->global_T_image.linear().transpose() * images[i + 1]->global_T_image.linear()).cast<double>()).transpose();
		Vector3d t_gt = -R1 * (images[i + 1]->global_T_image.translation() - images[i]->global_T_image.translation()).cast<double>();
		Vector3d t_with_scale = t_gt;
		t_gt.normalize(); // only the bearing vector is needed

		Matrix3d R_estimated;
		Vector3d t_estimated;
		double r_this, t_this;

		/* solve the epipolar problem here */
		// ConsistentEst est(cameras[images[i]->camera_id]->intrinsic);
		// est.GetPose(R_estimated, t_estimated, y_h, z_h);
		// R_estimated = est.R_svd;
		// t_estimated = est.t_svd;

		/* ↓ ↓ ↓ ↓ ↓ ↓---five point RANSAC method--- ↓ ↓ ↓ ↓ ↓ ↓ ↓ */
		vector<cv::Point2d> ycv, zcv;
		for (int k = 0; k < y_h.size(); ++k)
		{
			ycv.push_back(cv::Point2d(y_h[k].x(), y_h[k].y()));
			zcv.push_back(cv::Point2d(z_h[k].x(), z_h[k].y()));
		}
		Mat E_cv, intrinsic_cv, R_cv, t_cv;
		eigen2cv(cameras[images[i]->camera_id]->intrinsic, intrinsic_cv);
		E_cv = findEssentialMat(ycv, zcv, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
		recoverPose(E_cv, ycv, zcv, intrinsic_cv, R_cv, t_cv);
		cv2eigen(R_cv, R_estimated);
		cv2eigen(t_cv, t_estimated);
		/* ↑↑↑↑↑↑↑↑↑↑↑↑--five point RANSAC method--↑↑↑↑↑↑↑↑↑↑↑↑ */

		std::cout << "[R_est]\n"
				  << R_estimated << endl;
		std::cout << "[GT_R]\n"
				  << R_gt << endl;
		std::cout << "[T_est] " << t_estimated.transpose() << endl;
		std::cout << "[GT_T] " << t_gt.transpose() << endl;
		std::cout << "[scale T]" << t_with_scale.transpose() << endl;
		std::cout << "[norm of T]" << t_with_scale.norm() << endl;

		r_this = (R_estimated - R_gt).norm();
		t_this = abs(t_estimated.dot(t_gt));
		frobenius_error += r_this;
		total_cosine += 1/t_this;

		// todo : export those value for plotting
	}

	std::cout << "frobenius_error(L2):" << frobenius_error << endl;
	std::cout << "total_cosine(smaller better):" << total_cosine << endl;
	// Load 3d points ,index by point3d_id
	// ColmapPoint3DPtrMap point3ds;
	// success = ReadColmapPt3d(points3d_txt_path, point3ds);
	// if (success)
	// 	std::cout << "Successfully loaded " << point3ds.size() << " point3d(s)." << std::endl;
	// else
	// 	std::cout << "Error: could not load point3ds." << std::endl;

	return EXIT_SUCCESS;
}
