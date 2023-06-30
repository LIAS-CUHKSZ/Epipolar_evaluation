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

#include "draw_corre.hpp"
#include "cameras.h"
#include "images.h"
#include "pt3d.hpp"
#include "const_est.hpp"
#include "manifold_GN.hpp"
#include "npt_pose.h"
#include "eigensolver_wrapper.hpp"
#include "utils.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
	if (argc < 3)
	{
		std::cout << "Usage: " << argv[0] << "  <datasetname>  <img_windows_size>" << std::endl;
		std::cout << "you should put your dataset in the <dataset> dir of this proj";
		return EXIT_FAILURE; // <path_to_pt3ds_txt>
	}
	bool debug_img = false;
	if (argc == 4)
	{
		debug_img = true;
	}
	string dtset(argv[1]), windows_size_str(argv[2]);
	int wds = stoi(windows_size_str) - 1;
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
	double time_elapse;
	eval c_est("c_est"), e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"), pt7ransac("pt7ransac"), egsolver("eigenSolver");
	int valid_round = 0; // the number of valid pairs, only when the number of covisible points is larger than 20, the pair is valid
	// int progress = images.size() * wds - (wds * wds + wds) / 2, tmp_finish = 0;
	int progress = images.size() * (images.size() - 1) / 2, tmp_finish = 0;
	string dir_name = getTimeDir(dtset);

	vector<Vector3d> y_pix, z_pix, y_n, z_n;
	vector<Point2d> y_cv_pix, z_cv_pix;

	for (int img1 = 1; img1 < images.size(); ++img1)
	{
		// find the covisible points
		// surfix _pix means homogeneous coord in pixel plane, _n means normalized coordinates
		// img2 <= img1 + wds
		for (size_t img2 = img1 + 1; img2 <= images.size(); ++img2, ++tmp_finish)
		{
			printProg(tmp_finish, progress);
			if (images[img1]->camera_id != images[img2]->camera_id) // only consider imgs with the same camera for convenience
				continue;
			y_pix.clear();
			z_pix.clear();
			y_n.clear();
			z_n.clear();
			y_cv_pix.clear();
			z_cv_pix.clear();

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
			int total_covisible = y_pix.size();
			if (total_covisible < 200)
			{
				// std::cout << "[point too few,only " << y_pix.size() << " available,skip]" << endl;
				continue;
			}

			// compute the ground truth relative pose
			Matrix3d R1 = images[img1]->global_T_image.linear().transpose().cast<double>();
			Matrix3d R_gt = ((images[img1]->global_T_image.linear().transpose() * images[img2]->global_T_image.linear()).cast<double>()).transpose();
			Vector3d t_gt = -R_gt * R1 * (images[img2]->global_T_image.translation() - images[img1]->global_T_image.translation()).cast<double>();
			Vector3d t_with_scale = t_gt;
			if (t_with_scale.norm() < 0.075)
			{
				continue;
			}
			++valid_round;
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
			double r_err_this_round, t_err_this_round;

			/* ---------------------------remove outlier if necessary--------------------------------------*/
			// vector<Vector3d> y_pix_in, z_pix_in, y_n_in, z_n_in; // inliers
			// vector<Point2d> y_cv_pix_in, z_cv_pix_in;
			// cv::Mat inlier_mask, intri_cv;
			// cv::eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intri_cv);
			// cv::findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, cv::RANSAC, 0.999, 3.0, inlier_mask);
			// for (int i = 0; i < inlier_mask.rows; ++i) // remove outlier using inlier_mask
			// {
			// 	if (inlier_mask.at<uchar>(i, 0))
			// 	{
			// 		y_pix_in.push_back(y_pix[i]);
			// 		z_pix_in.push_back(z_pix[i]);
			// 		y_n_in.push_back(y_n[i]);
			// 		z_n_in.push_back(z_n[i]);
			// 		y_cv_pix_in.push_back(y_cv_pix[i]);
			// 		z_cv_pix_in.push_back(z_cv_pix[i]);
			// 	}
			// }
			// int inlier_num = y_pix_in.size();

			// get img pair name
			string img1path = GetFileName(images[img1]->file_path), img2path = GetFileName(images[img2]->file_path);
			string img1showpath = "../dataset/" + dtset + "/images/" + images[img1]->file_path;
			string img2showpath = "../dataset/" + dtset + "/images/" + images[img2]->file_path;
			//-------------------------------------- solve the epipolar problem here------------------------------------------

			/* ↓------------------consistent estimator------------------↓ */
			ConsistentEst est(cameras[images[img1]->camera_id]->intrinsic);
			time_elapse = TIME_IT(est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix););
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = (t_estimated - t_gt).norm();
			calcEval(R_lie, t_with_scale, c_est, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse, est.var_est);
			if ((r_err_this_round > 0.02 || t_err_this_round > 0.01) && debug_img)
				DrawCorreImg(dir_name, img1showpath, img2showpath, y_cv_pix, z_cv_pix, valid_round, true);

			/* ↑------------------consistent estimator------------------↑ */

			Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars
			eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intrinsic_cv);
			/* ↓------------------RANSAC method------------------↓ */
			time_elapse = TIME_IT(E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
								  recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
								  cv2eigen(R_cv, R_estimated);
								  cv2eigen(t_cv, t_estimated););
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = (t_estimated - t_gt).norm();
			calcEval(R_lie, t_with_scale, pt5ransac, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse);

			/* ↓------------------eigensolver estimator------------------↓ */
			//  @todo
			//  there is a problem with eigensolver, sometimes the translation recover from E is reversed
			//  after modifying the implementation of translation direction test, it still exists.
			// now just calc the Essnetial then use cv::recoverPose() to distinguish the direction
			eigenSolverWrapper gv_esv(y_n, z_n);
			time_elapse += TIME_IT(gv_esv.GetPose(R_estimated, t_estimated, R_estimated.transpose()););
			Matrix3d t_est_skew;
			t_est_skew << 0, -t_estimated(2), t_estimated(1),
				t_estimated(2), 0, -t_estimated(0),
				-t_estimated(1), t_estimated(0), 0;
			Matrix3d E_recover = t_est_skew * R_estimated;
			r_err_this_round = (R_estimated.transpose() - R_gt).norm();
			t_err_this_round = (-R_estimated.transpose() * t_estimated - t_gt).norm();
			if (t_err_this_round > 1.9)
			{
				t_err_this_round = (R_estimated.transpose() * t_estimated - t_gt).norm();
			}
			calcEval(R_lie, t_with_scale, egsolver, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse, est.var_est);
			/* ↑------------------eigensolver estimator------------------↑ */

			/* ↓------------------E Manifold GN------------------↓ */
			Matrix3d E_MGN_Init;
			cv2eigen(E_cv, E_MGN_Init);
			ManifoldGN MGN(cameras[images[img1]->camera_id]->intrinsic);
			time_elapse = TIME_IT(MGN.GetPose(R_estimated, t_estimated, y_cv_pix, z_cv_pix, y_n, z_n, E_MGN_Init););
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = (t_estimated - t_gt).norm();
			calcEval(R_lie, t_with_scale, e_m_gn, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse);
			/* ↑------------------E Manifold GN------------------↑ */

			// LMEDS
			time_elapse = TIME_IT(E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::LMEDS, 0.999, 1.0);
								  recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
								  cv2eigen(R_cv, R_estimated);
								  cv2eigen(t_cv, t_estimated););
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = (t_estimated - t_gt).norm();
			calcEval(R_lie, t_with_scale, pt7ransac, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse);
			/* ↑------------------RANSAC method------------------↑ */

			/* ↓------------------SDP on Essential Mat------------------↓ */
			// prepare data
			double *npt_p1 = new double[3 * y_n.size()], *npt_p2 = new double[3 * y_n.size()];
			for (size_t i = 0; i < y_n.size(); ++i)
			{
				y_n[i].normalize();
				z_n[i].normalize();
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
			time_elapse = TIME_IT(npt_pose(npt_p2, npt_p1, C, y_n.size(), X_sol, sdp_E, true);
								  eigen2cv(sdp_E, E_cv);
								  cv::recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
								  cv2eigen(R_cv, R_estimated);
								  cv2eigen(t_cv, t_estimated););
			r_err_this_round = (R_estimated - R_gt).norm();
			t_err_this_round = (t_estimated - t_gt).norm();
			calcEval(R_lie, t_with_scale, sdp, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse);
			/* ↑------------------SDP on Essential Mat------------------↑ */
		}
	}

	// ------------------------------------ save results here----------------------------------------
	std::cout << "----------------" << endl;
	std::cout << std::setw(15) << "[c_est] R:" << std::setw(12) << c_est.total_R_Fn << std::setw(12) << " t: " << c_est.total_t_cos << endl;
	std::cout << std::setw(15) << "[sdp] R:" << std::setw(12) << sdp.total_R_Fn << std::setw(12) << " t: " << sdp.total_t_cos << endl;
	std::cout << std::setw(15) << "[egsolver] R:" << std::setw(12) << egsolver.total_R_Fn << std::setw(12) << " t: " << egsolver.total_t_cos << endl;
	std::cout << std::setw(15) << "[e_m_gn] R:" << std::setw(12) << e_m_gn.total_R_Fn << std::setw(12) << " t: " << e_m_gn.total_t_cos << endl;
	std::cout << std::setw(15) << "[pt5ransac] R:" << std::setw(12) << pt5ransac.total_R_Fn << std::setw(12) << " t: " << pt5ransac.total_t_cos << endl;
	std::cout << std::setw(15) << "[pt7ransac] R:" << std::setw(12) << pt7ransac.total_R_Fn << std::setw(12) << " t: " << pt7ransac.total_t_cos << endl;
	std::cout << "--------------------------------------" << endl;

	saveRes(c_est, dir_name);
	saveRes(sdp, dir_name);
	saveRes(egsolver, dir_name);
	saveRes(e_m_gn, dir_name);
	saveRes(pt5ransac, dir_name);
	saveRes(pt7ransac, dir_name);

	// save method error;
	std::ofstream file(dir_name + "/errors.txt");
	file << "method,      avr_time,      R_err ,      t_err" << endl;
	file << std::setw(12) << "c_est: " << std::setw(10) << c_est.average_time << std::setw(10) << c_est.total_R_Fn << std::setw(10) << c_est.total_t_cos << endl;
	file << std::setw(12) << "sdp: " << std::setw(10) << sdp.average_time << std::setw(10) << sdp.total_R_Fn << std::setw(10) << sdp.total_t_cos << endl;
	file << std::setw(12) << "egsolver: " << std::setw(10) << egsolver.average_time << std::setw(10) << egsolver.total_R_Fn << std::setw(10) << egsolver.total_t_cos << endl;
	file << std::setw(12) << "e_m_gn: " << std::setw(10) << e_m_gn.average_time << std::setw(10) << e_m_gn.total_R_Fn << std::setw(10) << e_m_gn.total_t_cos << endl;
	file << std::setw(12) << "pt5ransac: " << std::setw(10) << pt5ransac.average_time << std::setw(10) << pt5ransac.total_R_Fn << std::setw(10) << pt5ransac.total_t_cos << endl;
	file << std::setw(12) << "pt7ransac: " << std::setw(10) << pt7ransac.average_time << std::setw(10) << pt7ransac.total_R_Fn << std::setw(10) << pt7ransac.total_t_cos << endl;
	// if the error of c_est is the smallest of the five methods, print it
	int flag = 0;
	if (c_est.total_R_Fn < sdp.total_R_Fn && c_est.total_R_Fn < e_m_gn.total_R_Fn && c_est.total_R_Fn < pt5ransac.total_R_Fn && c_est.total_R_Fn < pt7ransac.total_R_Fn && c_est.total_R_Fn < egsolver.total_R_Fn)
	{
		file << "c_est sota R" << endl;
		flag |= 1;
	}
	if (c_est.total_t_cos < sdp.total_t_cos && c_est.total_t_cos < e_m_gn.total_t_cos && c_est.total_t_cos < pt5ransac.total_t_cos && c_est.total_t_cos < egsolver.total_t_cos && c_est.total_t_cos < egsolver.total_t_cos)
	{
		file << "c_est sota t" << endl;
		flag |= 2;
	}
	if (flag == 3)
	{
		file << "c_est best!!!" << endl;
	}
	file.close();

	std::string filename = "best_estimation.txt";
	std::ofstream eval_file;
	// save the best method
	std::ifstream infile(filename);
	if (infile.good())
		eval_file.open(filename, std::ios::app);
	else
		eval_file.open(filename);
	eval_file << flag << "  " << dtset << endl;
	eval_file.close();

	return EXIT_SUCCESS;
}

// Load 3d points ,index by point3d_id
// ColmapPoint3DPtrMap point3ds;
// success = ReadColmapPt3d(points3d_txt_path, point3ds);
// if (success)
// 	std::cout << "Successfully loaded " << point3ds.size() << " point3d(s)." << std::endl;
// else
// 	std::cout << "Error: could not load point3ds." << std::endl;
