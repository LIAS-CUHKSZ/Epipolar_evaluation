

#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <string>
#include <yaml-cpp/yaml.h>

#include "cameras.h"
#include "const_est.hpp"
#include "draw_corre.hpp"
#include "eigensolver_wrapper.hpp"
#include "images.h"
#include "lm_wrapper.hpp"
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

  // Load configuration using yaml-cpp
  YAML::Node config = YAML::LoadFile("../config.yaml");
  if (!config) {
    std::cerr << "Failed to open configuration file: ../config.yaml"
              << std::endl;
    return EXIT_FAILURE;
  }

  double t_min = config["config"]["eth3d"]["t_min"].as<double>();
  double r_min = config["config"]["eth3d"]["r_min"].as<double>();
  int pts_min = config["config"]["eth3d"]["pts_min"].as<int>();
  bool use_lie = config["config"]["eth3d"]["use_lie"].as<bool>();

  std::cout << "Configuration for eth3d loaded:" << std::endl;
  std::cout << "t_min: " << t_min << ", r_min: " << r_min
            << ", pts_min: " << pts_min << ", use_lie: " << use_lie
            << std::endl;

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
    std::cout << "Error: could not load image info.Exit" << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "[" << dtset << "]eval on going, please wait..." << endl;

  // accumulate the error
  eval c_est("c_est0"), c1("c_est1"), c2("c_est2"), c3("c_est3"), c4("c_est4"),
      c5("c_est5");
  // the number of valid pairs, only when the number of covisible points is
  // larger than 20 and the norm of translation is bigger than 0.075m will this
  // pair to be evaluate
  int valid_round = 0;
  int progress, tmp_finish = 0; // for printing progress bar
  if (wds == -2)
    progress = images.size() * (images.size() - 1) / 2;
  else
    progress = images.size() * wds - (wds * wds + wds) / 2;

  // create a folder under "build" named time(m-d-h-m)-dataset_name
  string dir_name = getTimeDir(dtset);

  // store point pairs
  vector<Vector3d> y_n, z_n;
  vector<Point2d> y_cv_pix, z_cv_pix;

  for (size_t img1 = 1; img1 < images.size(); ++img1) {
    // find the covisible points
    // surfix _pix means homogeneous coord in pixel plane, _n means normalized
    // coordinates(z=1)
    for (size_t img2 = img1 + 1;
         img2 <= images.size() && (img2 <= img1 + wds || wds == -2);
         ++img2, ++tmp_finish) {
      printProg(tmp_finish, progress);

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
      Vector3d t_with_scale = t_gt;
      Vector3d R_lie = unskew(R_gt);
      if (t_with_scale.norm() < t_min ||
          R_lie.norm() < r_min) // two-view geometry cannot evaluate the
                                // translation when it is too small
        continue;
      t_gt.normalize(); // only the bearing vector is needed
      // Matrix3d E_ground = skew(t_gt) * R_gt;

      y_n.clear();
      z_n.clear();
      y_cv_pix.clear();
      z_cv_pix.clear();

      // traverse obs in first img
      for (auto yi : images[img1]->observations) {
        if (images[img2]->observations.count(
                yi.first)) // if the point is visible in the another image
        {
          point_pair p1 = yi.second, p2 = images[img2]->observations[yi.first];

          y_cv_pix.emplace_back(Point2d(p1.pixel.x(), p1.pixel.y()));
          z_cv_pix.emplace_back(Point2d(p2.pixel.x(), p2.pixel.y()));
          y_n.emplace_back(p1.normalized);
          z_n.emplace_back(p2.normalized);
        }
      }
      int total_covisible = z_n.size();
      if (total_covisible < pts_min) // we focus on large number case, you can
                                     // modify it to a smaller num
        continue;
      ++valid_round;

      // temp variables to store res
      Matrix3d R_estimated;
      Vector3d t_estimated;
      double r_err_this_round, t_err_this_round, time_elapse = 0;

      // get img pair name
      string img1path =
                 to_string(img1) + "-" + GetFileName(images[img1]->file_path),
             img2path = to_string(img2) + "-" +
                        GetFileName(images[img2]->file_path); // for res saving
      // for debug img saving(name)
      string img1showpath =
          "../dataset/" + dtset + "/images/" + images[img1]->file_path;
      string img2showpath =
          "../dataset/" + dtset + "/images/" + images[img2]->file_path;

      //-------------------------------------- solve the epipolar problem
      // here------------------------------------------

      // temp vars
      Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars
      Matrix3d R_r5pt;
      Vector3d t_r5pt;
      eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intrinsic_cv);

      ConsistentEst est(cameras[images[img1]->camera_id]->intrinsic);
      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 0,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c_est, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);

      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 1,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c1, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);
      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 2,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c2, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);
      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 3,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c3, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);
      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 4,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c4, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);
      est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix, 5,
                  -1.0, true);
      calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
              t_estimated, use_lie);
      calcEval(R_lie, t_with_scale, c5, img1path, img2path, t_err_this_round,
               r_err_this_round, total_covisible, time_elapse, est.var_est);
    }
  }

  // ------------------------------------ save results
  // here----------------------------------------
  saveRes(c_est, dir_name);
  saveRes(c2, dir_name);
  saveRes(c5, dir_name);
  saveRes(c1, dir_name);
  saveRes(c4, dir_name);
  saveRes(c3, dir_name);

  std::cout << "\n---------------- Evaluation Summary ----------------"
            << std::endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << std::setw(15) << "Method" << std::setw(20) << "Avg Time (s)"
            << std::setw(15) << "R_err" << std::setw(15) << "t_err"
            << std::endl;

  std::cout << std::setw(15) << "c_est" << std::setw(20) << c_est.average_time
            << std::setw(15) << c_est.total_R_Fn << std::setw(15)
            << c_est.total_t_cos << std::endl;

  std::cout << std::setw(15) << "c1" << std::setw(20) << c1.average_time
            << std::setw(15) << c1.total_R_Fn << std::setw(15) << c1.total_t_cos
            << std::endl;

  std::cout << std::setw(15) << "c2" << std::setw(20) << c2.average_time
            << std::setw(15) << c2.total_R_Fn << std::setw(15) << c2.total_t_cos
            << std::endl;
  std::cout << std::setw(15) << "c3" << std::setw(20) << c3.average_time
            << std::setw(15) << c3.total_R_Fn << std::setw(15) << c3.total_t_cos
            << std::endl;
  std::cout << std::setw(15) << "c4" << std::setw(20) << c4.average_time
            << std::setw(15) << c4.total_R_Fn << std::setw(15) << c4.total_t_cos
            << std::endl;
  std::cout << std::setw(15) << "c5" << std::setw(20) << c5.average_time
            << std::setw(15) << c5.total_R_Fn << std::setw(15) << c5.total_t_cos
            << std::endl;

  std::cout << "----------------------------------------------------"
            << std::endl;

  // save method error;
  std::ofstream file(dir_name + "/errors.txt");
  file << std::setw(15) << "Method" << std::setw(12) << "Average Time (s)"
       << std::setw(12) << "R_err" << std::setw(12) << "t_err" << endl;
  file << std::setw(15) << "c_est: " << std::setw(15) << c_est.average_time
       << std::setw(15) << c_est.total_R_Fn << std::setw(15)
       << c_est.total_t_cos << endl;
  file << std::setw(15) << "c1: " << std::setw(15) << c1.average_time
       << std::setw(15) << c1.total_R_Fn << std::setw(15) << c1.total_t_cos
       << endl;
  file << std::setw(15) << "c2: " << std::setw(15) << c2.average_time
       << std::setw(15) << c2.total_R_Fn << std::setw(15) << c2.total_t_cos
       << endl;
  file << std::setw(15) << "c3: " << std::setw(15) << c3.average_time
       << std::setw(15) << c3.total_R_Fn << std::setw(15) << c3.total_t_cos
       << endl;
  file << std::setw(15) << "c4: " << std::setw(15) << c4.average_time
       << std::setw(15) << c4.total_R_Fn << std::setw(15) << c4.total_t_cos
       << endl;

  file << std::setw(15) << "c5: " << std::setw(15) << c5.average_time
       << std::setw(15) << c5.total_R_Fn << std::setw(15) << c5.total_t_cos
       << endl;

  return EXIT_SUCCESS;
}