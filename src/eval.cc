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

#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

#include "cameras.h"
#include "const_est.hpp"
#include "draw_corre.hpp"
#include "eigensolver_wrapper.hpp"
#include "images.h"
#include "lm_method_wrapper.hpp"
#include "manifold_GN.hpp"
#include "npt_pose.h"
#include "pt3d.hpp"
#include "utils.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << "  <datasetname>  <img_windows_size>"
              << std::endl;
    std::cout << "you should put your dataset in the <dataset> dir of this proj"
              << std::endl;
    std::cout << "if <img_window_size> is set to -1, then each pair of img in "
                 "the dataset will be taken into account"
              << endl;
    return EXIT_FAILURE;
  }

  bool debug_img = false;
  if (argc == 4) // input the 4th param(anything is ok) to enable debug image
    debug_img = true;
  string dtset(argv[1]), windows_size_str(argv[2]);
  int wds = stoi(windows_size_str) - 1;
  if (wds == -2)
    std::cout << "evaluate all img pairs" << std::endl;
  string cameras_txt_path =
      "../dataset/" + dtset + "/dslr_calibration_undistorted/cameras.txt";
  string images_txt_path =
      "../dataset/" + dtset + "/dslr_calibration_undistorted/images.txt";

  // Load cameras (indexed by: camera_id).
  ColmapCameraPtrMap cameras;
  bool success = ReadColmapCameras(cameras_txt_path, &cameras);
  if (success)
    std::cout << "Successfully loaded " << cameras.size() << " camera(s)."
              << std::endl;
  else {
    std::cout << "Error: could not load cameras.Exit" << std::endl;
    return EXIT_FAILURE;
  }

  // Load images (indexed by: image_id).
  ColmapImagePtrMap images;
  success = ReadColmapImages(images_txt_path, true, &images, cameras);
  if (success)
    std::cout << "Successfully loaded " << images.size() << " image info(s)."
              << std::endl;
  else {
    std::cout << "Error: could not load image info. Exit." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "[" << dtset << "]eval on going..." << endl;

  int valid_round = 0;

  // under "build" named time(m-d-h-m)-dataset_name
  string dir_name = getTimeDir(dtset);

  // store point pairs
  vector<Vector3d> y_n, z_n;
  vector<Point2d> y_cv_pix, z_cv_pix;

  ResultAggregator results;
    vector<unique_ptr<PoseSolver>> solvers;
    solvers.emplace_back(make_unique<ManifoldGN>());

  for (size_t img1 = 1; img1 < images.size(); ++img1) {
    // surfix _pix means homogeneous coord in pixel plane
    // _n means normalized coordinates(z=1)
    for (size_t img2 = img1 + 1;
         img2 <= images.size() && (img2 <= img1 + wds || wds == -2); ++img2) {

      ImagePair eval_pair;
      eval_pair.img1_id = img1;
      eval_pair.img2_id = img2;
      // compute the ground truth relative pose
      Matrix3d R1 =
          images[img1]->global_T_image.linear().transpose().cast<double>();
      Matrix3d R_gt = ((images[img1]->global_T_image.linear().transpose() *
                        images[img2]->global_T_image.linear())
                           .cast<double>())
                          .transpose();
      Vector3d t_gt = -R_gt * R1 *
                      (images[img2]->global_T_image.translation() -
                       images[img1]->global_T_image.translation())
                          .cast<double>();
      if (t_gt.norm() < 0.075)
        continue;
      eval_pair.t_gt = t_gt.normalized();
      eval_pair.R_gt = R_gt;

      // traverse obs in first img
      for (auto yi : images[img1]->observations)
        if (images[img2]->observations.count(
                yi.first)) // if the point is visible in the another image
        {
          point_pair p1 = yi.second, p2 = images[img2]->observations[yi.first];
          eval_pair.y_cv_pix.emplace_back(Point2d(p1.pixel.x(), p1.pixel.y()));
          eval_pair.z_cv_pix.emplace_back(Point2d(p2.pixel.x(), p2.pixel.y()));
          eval_pair.y_n.emplace_back(p1.normalized);
          eval_pair.z_n.emplace_back(p2.normalized);
        }
      if (eval_pair.y_cv_pix.size() < 200) // we focus on large number case
        continue;
      ++valid_round;

      // get img pair name
      eval_pair.img1_path =
          to_string(img1) + "-" + GetFileName(images[img1]->file_path);
      eval_pair.img2_path =
          to_string(img2) + "-" + GetFileName(images[img2]->file_path);

    //   string img1showpath =
    //       "../dataset/" + dtset + "/images/" + images[img1]->file_path;
    //   string img2showpath =
    //       "../dataset/" + dtset + "/images/" + images[img2]->file_path;

      // remove outlier if necessary
      // @note the low res datasets have much more outlier than high-res you can
      // add definition REMOVE_OUTLIER whenu sing low-res
#ifdef REMOVE_OUTLIER
      vector<Vector3d> y_n_in, z_n_in; // inliers
      vector<Point2d> y_cv_pix_in, z_cv_pix_in;
      cv::Mat inlier_mask, intri_cv;
      cv::eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intri_cv);
      cv::findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, cv::LMEDS, 0.999, 20,
                           inlier_mask);
      for (size_t i = 0; i < inlier_mask.rows;
           ++i) // remove outlier using inlier_mask
      {
        if (inlier_mask.at<uchar>(i, 0)) {
          y_n_in.emplace_back(y_n[i]);
          z_n_in.emplace_back(z_n[i]);
          y_cv_pix_in.emplace_back(y_cv_pix[i]);
          z_cv_pix_in.emplace_back(z_cv_pix[i]);
        }
      }
      y_n = std::move(y_n_in);
      z_n = std::move(z_n_in);
      y_cv_pix = std::move(y_cv_pix_in);
      z_cv_pix = std::move(z_cv_pix_in);
      total_covisible = z_cv_pix.size();
#endif // REMOVE_OUTLIER
    }


  }

  return EXIT_SUCCESS;
}