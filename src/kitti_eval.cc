#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "const_est.hpp"
#include "draw_corre.hpp"
#include "eigensolver_wrapper.hpp"
#include "lm_wrapper.hpp"
#include "manifold_GN.hpp"
#include "npt_pose.h"
#include "pt3d.hpp"
#include "utils.hpp"

using namespace cv;
using namespace std;
using namespace Eigen;

// Structure to hold point correspondences
struct PointCorrespondence {
  int img1_idx;
  int img2_idx;
  vector<Point2d> points1;
  vector<Point2d> points2;
};

// Function to load point correspondences from matches.txt
vector<PointCorrespondence>
loadPointCorrespondences(const string &matchesFile) {
  vector<PointCorrespondence> correspondences;
  ifstream file(matchesFile);

  if (!file.is_open()) {
    cerr << "Failed to open matches file: " << matchesFile << endl;
    return correspondences;
  }

  string line;
  while (getline(file, line)) {
    PointCorrespondence corr;
    stringstream ss(line);
    string token;

    // Parse image indices
    getline(ss, token, ',');
    stringstream idxStream(token);
    idxStream >> corr.img1_idx >> corr.img2_idx;

    // Parse points from first image
    getline(ss, token, ',');
    stringstream pts1Stream(token);
    double x, y;
    while (pts1Stream >> x >> y) {
      corr.points1.emplace_back(x, y);
    }

    // Parse points from second image
    getline(ss, token);
    stringstream pts2Stream(token);
    while (pts2Stream >> x >> y) {
      corr.points2.emplace_back(x, y);
    }

    if (corr.points1.size() == corr.points2.size() && !corr.points1.empty()) {
      correspondences.push_back(corr);
    }
  }

  return correspondences;
}

// Function to load ground truth poses from file
vector<Matrix<double, 3, 4>> loadGroundTruthPoses(const string &posesFile) {
  vector<Matrix<double, 3, 4>> poses;
  ifstream file(posesFile);

  if (!file.is_open()) {
    cerr << "Failed to open poses file: " << posesFile << endl;
    return poses;
  }

  string line;
  while (getline(file, line)) {
    stringstream ss(line);
    Matrix<double, 3, 4> pose;

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 4; ++j) {
        ss >> pose(i, j);
      }
    }

    poses.push_back(pose);
  }

  return poses;
}

// Convert pixel coordinates to normalized coordinates
void pixelToNormalized(const vector<Point2d> &pixelPoints,
                       vector<Vector3d> &normalizedPoints, const Matrix3d &K) {

  for (const auto &p : pixelPoints) {
    double x = (p.x - K(0, 2)) / K(0, 0);
    double y = (p.y - K(1, 2)) / K(1, 1);
    Vector3d normalized(x, y, 1.0);
    normalizedPoints.emplace_back(normalized);
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <sequence_number>  <eval_dir> "
              << std::endl;
    std::cout << "Example: " << argv[0] << " 00  eval" << std::endl;
    return EXIT_FAILURE;
  }

  // Load configuration using yaml-cpp
  YAML::Node config = YAML::LoadFile("../config.yaml");
  if (!config) {
    std::cerr << "Failed to open configuration file: ../config.yaml"
              << std::endl;
    return EXIT_FAILURE;
  }

  double t_min = config["config"]["kitti"]["t_min"].as<double>();
  double r_min = config["config"]["kitti"]["r_min"].as<double>();
  int pts_min = config["config"]["kitti"]["pts_min"].as<int>();
  bool use_lie = config["config"]["kitti"]["use_lie"].as<bool>();

  std::cout << "Configuration for kitti loaded:" << std::endl;
  std::cout << "t_min: " << t_min << ", r_min: " << r_min
            << ", pts_min: " << pts_min << ", use_lie: " << use_lie
            << std::endl;

  string sequenceNum(argv[1]);
  string eval_dir(argv[2]);
  bool debug_img = false;

  if (argc == 4) // input the 4th param(anything is ok) to enable debug image
    debug_img = true;

  cout << "Evaluating KITTI sequence " << sequenceNum << "..." << endl;

  // Set paths
  string seqDir = "../dataset/" + sequenceNum;
  string matchesFile = seqDir + "/" + eval_dir + "/matches.txt";
  string posesFile = "../dataset/poses/" + sequenceNum + "_rel.txt";

  // use p0
  Matrix3d K = Matrix3d::Identity();
  K(0, 0) = 718.856; // fx
  K(1, 1) = 718.856; // fy
  K(0, 2) = 607.193; // cx
  K(1, 2) = 185.216; // cy
  cout << "Camera intrinsics:" << endl << K << endl;

  vector<PointCorrespondence> correspondences =
      loadPointCorrespondences(matchesFile);
  cout << "Loaded " << correspondences.size() << " point correspondences"
       << endl;
  vector<Matrix<double, 3, 4>> gtPoses = loadGroundTruthPoses(posesFile);
  cout << "Loaded " << gtPoses.size() << " ground truth poses" << endl;

  if (correspondences.size() != gtPoses.size()) {
    cout << "Warning: Number of point correspondences ("
         << correspondences.size() << ") doesn't match number of poses ("
         << gtPoses.size() << ")" << endl;
  }
  // gtPoses.resize(50);
  // correspondences.resize(50);

  // create a folder under "build" named time(m-d-h-m)-dataset_name
  string dir_name = getTimeDir("kitti-" + sequenceNum + "-" + eval_dir);

  // accumulate the error
  eval c_est("c_est"), e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"),
      egsolver("eigenSolver"), lm("lmsolver");

  // the number of valid pairs
  int valid_round = 0;
  int progress = correspondences.size();
  int tmp_finish = 0; // for printing progress bar

  long long total_covisible_average = 0;

  for (size_t i = 0; i < correspondences.size() && i < gtPoses.size(); ++i) {
    printProg(tmp_finish++, progress);

    const PointCorrespondence &corr = correspondences[i];

    // Skip if too few points
    if (corr.points1.size() < pts_min) {
      continue;
    }

    // Get ground truth pose
    Matrix<double, 3, 4> gt_pose = gtPoses[i];
    Matrix3d R_gt = gt_pose.block<3, 3>(0, 0).transpose();
    Vector3d R_lie = unskew(R_gt);
    Vector3d t_gt = -R_gt * gt_pose.block<3, 1>(0, 3);
    Vector3d t_with_scale = t_gt;

    // Skip if translation or rotation is too small
    if (t_with_scale.norm() < t_min || R_lie.norm() < r_min) {
      continue;
    }

    t_gt.normalize(); // only the bearing vector is needed

    // Convert pixel coordinates to normalized coordinates
    vector<Vector3d> y_n, z_n;
    vector<Point2d> y_cv_pix = corr.points1;
    vector<Point2d> z_cv_pix = corr.points2;

    pixelToNormalized(y_cv_pix, y_n, K);
    pixelToNormalized(z_cv_pix, z_n, K);

    int total_covisible = z_n.size();
    total_covisible_average += total_covisible;
    ++valid_round;

    // temp variables to store res
    Matrix3d R_estimated;
    Vector3d t_estimated;
    double r_err_this_round, t_err_this_round, time_elapse;

    // get img pair name
    string img1path = to_string(corr.img1_idx),
           img2path = to_string(corr.img2_idx);

    // for debug img saving(name)
    string img1showpath =
        seqDir + "/image_0/" + to_string(corr.img1_idx) + ".png";
    string img2showpath =
        seqDir + "/image_0/" + to_string(corr.img2_idx) + ".png";

    // temp vars
    Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars
    Matrix3d R_r5pt;
    Vector3d t_r5pt;
    eigen2cv(K, intrinsic_cv);

    /* ↓------------------RANSAC-5pt method------------------↓ */
    double ransac_time =
        TIME_IT(E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv,
                                        cv::RANSAC, .999, 1.);
                recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                cv2eigen(R_cv, R_r5pt); cv2eigen(t_cv, t_r5pt););
    time_elapse = ransac_time;
    R_estimated = R_r5pt;
    t_estimated = t_r5pt;
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_r5pt, t_r5pt,
            use_lie); // Use calcErr
    calcEval(R_lie, t_with_scale, pt5ransac, img1path, img2path,
             t_err_this_round, r_err_this_round, total_covisible, time_elapse);

    /* ↓------------------consistent estimator------------------↓ */
    ConsistentEst est(K);
    time_elapse = TIME_IT(est.GetPose(R_estimated, t_estimated, y_n, z_n,
                                      y_cv_pix, z_cv_pix, 1, 0.008, true););
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
            t_estimated, use_lie); // Use calcErr
    calcEval(R_lie, t_with_scale, c_est, img1path, img2path, t_err_this_round,
             r_err_this_round, total_covisible, time_elapse, est.var_est);
    if ((r_err_this_round > 0.02 || t_err_this_round > 0.01) && debug_img)
      DrawCorreImg(dir_name, img1showpath, img2showpath, y_cv_pix, z_cv_pix,
                   valid_round, false);
    /* ↑------------------consistent estimator------------------↑ */

    EigenWrapper gv_esv(y_n, z_n);
    time_elapse = TIME_IT(gv_esv.GetPose(R_estimated, t_estimated, R_r5pt);) +
                  ransac_time;
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
            t_estimated, use_lie); // Use calcErr
    calcEval(R_lie, t_with_scale, egsolver, img1path, img2path,
             t_err_this_round, r_err_this_round, total_covisible, time_elapse);

    /* ↑------------------eigensolver estimator------------------↑ */

    // lm solver, using 5pt as init value
    lmSolverWrapper lm_solver(y_n, z_n);
    time_elapse =
        TIME_IT(lm_solver.GetPose(R_estimated, t_estimated, R_r5pt, t_r5pt);) +
        ransac_time;
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
            t_estimated, use_lie);
    calcEval(R_lie, t_with_scale, lm, img1path, img2path, t_err_this_round,
             r_err_this_round, total_covisible, time_elapse);

    /* ↓------------------E Manifold GN------------------↓ */
    Matrix3d E_MGN_Init;
    cv2eigen(E_cv, E_MGN_Init); // use RANSAC-5pt result as initial value
    ManifoldGN MGN(K);
    time_elapse = TIME_IT(MGN.GetPose(R_estimated, t_estimated, y_cv_pix,
                                      z_cv_pix, y_n, z_n, E_MGN_Init);) +
                  ransac_time;
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
            t_estimated, use_lie);
    calcEval(R_lie, t_with_scale, e_m_gn, img1path, img2path, t_err_this_round,
             r_err_this_round, total_covisible, time_elapse);
    /* ↑------------------E Manifold GN------------------↑ */

    /* ↓------------------SDP on Essential Mat------------------↓ */
    // prepare data
    double *npt_p1 = new double[3 * y_n.size()],
           *npt_p2 = new double[3 * y_n.size()];
    for (size_t i = 0; i < y_n.size(); ++i) {
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
    time_elapse = TIME_IT(
        npt_pose(npt_p2, npt_p1, C, y_n.size(), X_sol, sdp_E, true);
        eigen2cv(sdp_E, E_cv);
        cv::recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
        cv2eigen(R_cv, R_estimated); cv2eigen(t_cv, t_estimated););
    calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated,
            t_estimated, use_lie); // Use calcErr
    calcEval(R_lie, t_with_scale, sdp, img1path, img2path, t_err_this_round,
             r_err_this_round, total_covisible, time_elapse);
    /* ↑------------------SDP on Essential Mat------------------↑ */

    // Clean up
    delete[] npt_p1;
    delete[] npt_p2;
    delete[] C;
  }

  // ------------------------------------ save results
  // here----------------------------------------
  saveRes(c_est, dir_name);
  saveRes(sdp, dir_name);
  saveRes(egsolver, dir_name);
  saveRes(e_m_gn, dir_name);
  saveRes(lm, dir_name);
  saveRes(pt5ransac, dir_name);

  std::cout << endl
            << "---------------- Time and Error Summary ----------------"
            << endl;
  std::cout << std::fixed << std::setprecision(6);
  std::cout << std::setw(15) << "Method" << std::setw(20) << "Average Time (s)"
            << std::setw(15) << "R_err" << std::setw(15) << "t_err" << endl;

  std::cout << std::setw(15) << "c_est" << std::setw(20) << c_est.average_time
            << std::setw(15) << c_est.total_R_Fn << std::setw(15)
            << c_est.total_t_cos << endl;

  std::cout << std::setw(15) << "sdp" << std::setw(20) << sdp.average_time
            << std::setw(15) << sdp.total_R_Fn << std::setw(15)
            << sdp.total_t_cos << endl;

  std::cout << std::setw(15) << "egsolver" << std::setw(20)
            << egsolver.average_time << std::setw(15) << egsolver.total_R_Fn
            << std::setw(15) << egsolver.total_t_cos << endl;

  std::cout << std::setw(15) << "e_m_gn" << std::setw(20) << e_m_gn.average_time
            << std::setw(15) << e_m_gn.total_R_Fn << std::setw(15)
            << e_m_gn.total_t_cos << endl;

  std::cout << std::setw(15) << "lm" << std::setw(20) << lm.average_time
            << std::setw(15) << lm.total_R_Fn << std::setw(15) << lm.total_t_cos
            << endl;

  std::cout << std::setw(15) << "pt5ransac" << std::setw(20)
            << pt5ransac.average_time << std::setw(15) << pt5ransac.total_R_Fn
            << std::setw(15) << pt5ransac.total_t_cos << endl;
  std::cout << "--------------------------------------------------------"
            << endl;

  // save method error;
  std::ofstream file(dir_name + "/errors.txt");
  file << std::setw(15) << "Method" << std::setw(12) << "Average Time (s)"
       << std::setw(12) << "R_err" << std::setw(12) << "t_err" << endl;
  file << std::fixed << std::setw(15) << "c_est: " << std::setw(12)
       << c_est.average_time << std::setw(12) << c_est.total_R_Fn
       << std::setw(12) << c_est.total_t_cos << endl;
  file << std::fixed << std::setw(15) << "sdp: " << std::setw(12)
       << sdp.average_time << std::setw(12) << sdp.total_R_Fn << std::setw(12)
       << sdp.total_t_cos << endl;
  file << std::fixed << std::setw(15) << "egsolver: " << std::setw(12)
       << egsolver.average_time << std::setw(12) << egsolver.total_R_Fn
       << std::setw(12) << egsolver.total_t_cos << endl;
  file << std::fixed << std::setw(15) << "e_m_gn: " << std::setw(12)
       << e_m_gn.average_time << std::setw(12) << e_m_gn.total_R_Fn
       << std::setw(12) << e_m_gn.total_t_cos << endl;
  file << std::fixed << std::setw(15) << "lm: " << std::setw(12)
       << lm.average_time << std::setw(12) << lm.total_R_Fn << std::setw(12)
       << lm.total_t_cos << endl;
  file << std::fixed << std::setw(15) << "pt5ransac: " << std::setw(12)
       << pt5ransac.average_time << std::setw(12) << pt5ransac.total_R_Fn
       << std::setw(12) << pt5ransac.total_t_cos << endl;
  file << "average_covisible: " << (double)total_covisible_average / valid_round
       << endl;

  // if the error of c_est is the smallest of the five methods, print it
  int flag = 0;
  if (c_est.total_R_Fn < sdp.total_R_Fn && c_est.total_R_Fn < lm.total_R_Fn &&
      c_est.total_R_Fn < e_m_gn.total_R_Fn &&
      c_est.total_R_Fn < egsolver.total_R_Fn) {
    flag |= 1;
  }
  if (c_est.total_t_cos < sdp.total_t_cos &&
      c_est.total_t_cos < lm.total_t_cos &&
      c_est.total_t_cos < e_m_gn.total_t_cos &&
      c_est.total_t_cos < egsolver.total_t_cos) {
    flag |= 2;
  }

  if (flag == 3) {
    file << "c_est best!!!" << endl;
    cout << "c_est best!!!" << endl;
  } else if (flag == 1) {
    file << "c_est sota R" << endl;
    cout << "c_est sota R" << endl;
  } else if (flag == 2) {
    file << "c_est sota t" << endl;
    cout << "c_est sota t" << endl;
  } else {
    file << "no best" << endl;
    cout << "no best" << endl;
  }
  file.close();

  std::string filename = "kitti_best_estimation.txt";
  std::ofstream eval_file;
  // save the best method
  std::ifstream infile(filename);
  if (infile.good())
    eval_file.open(filename, std::ios::app);
  else
    eval_file.open(filename);
  eval_file << flag << "  "
            << sequenceNum + " " + eval_dir + (use_lie ? " true" : " false")
            << endl;
  eval_file.close();

  return EXIT_SUCCESS;
}