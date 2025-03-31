#pragma once
#include "methods/PoseSovlerWrapper.hpp"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>

#include <Eigen/Dense>

namespace fs = std::filesystem;
using namespace std;
using namespace Eigen;
using namespace cv;

#define TIME_IT(code)                                                          \
  [&]() -> long long {                                                         \
    auto start_time = std::chrono::high_resolution_clock::now();               \
    code;                                                                      \
    auto end_time = std::chrono::high_resolution_clock::now();                 \
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(     \
        end_time - start_time);                                                \
    return duration.count();                                                   \
  }()

Eigen::Matrix3d skew(Eigen::Vector3d a) {
  Eigen::Matrix3d skew_mat;
  skew_mat << 0, -a(2), a(1), a(2), 0, -a(0), -a(1), a(0), 0;
  return skew_mat;
}

Eigen::Vector3d unskew(Eigen::Matrix3d M) {
  Eigen::Matrix3d R_log = M.log().eval();
  return Eigen::Vector3d(R_log(2, 1), R_log(0, 2), R_log(1, 0));
}

struct PoseResult {
  Matrix3d R_estimated;
  Vector3d t_estimated;
  double time_taken;
  double r_error;
  double t_error;
};

struct ImagePair {
  size_t img1_id;
  size_t img2_id;
  vector<Eigen::Vector3d> y_n, z_n;
  vector<cv::Point2d> y_cv_pix, z_cv_pix;
  Matrix3d R_gt;
  Vector3d t_gt;
  string img1_path;
  string img2_path;
};

struct eval {
  std::string method_name;
  std::vector<std::pair<std::string, std::string>> img_pair;
  std::vector<double> t_err_per_round;
  std::vector<double> R_err_per_round;
  std::vector<double> t_gt_norm;
  std::vector<Eigen::Vector3d> R_gt;
  double total_R_Fn = 0;
  double total_t_cos = 0;
  double average_time = 0;
  std::vector<int> num_pts;
  std::vector<double> time;
  std::vector<double> noise;
  eval(std::string name) : method_name(name){};
};

void calcEval(Eigen::Vector3d rgt, Eigen::Vector3d tgt, eval &Method,
              string img_idx1, string img_idx2, double t_err, double r_err,
              int num, double time_cost, double noise = -1) {
  Method.img_pair.push_back(std::make_pair(img_idx1, img_idx2));
  Method.t_err_per_round.push_back(t_err);
  Method.R_err_per_round.push_back(r_err);
  Method.total_R_Fn += r_err;
  Method.total_t_cos += t_err;
  if (noise != -1)
    Method.noise.push_back(noise);
  Method.num_pts.push_back(num);
  Method.time.push_back(time_cost);
  Method.average_time += time_cost;
  Method.R_gt.push_back(rgt);
  Method.t_gt_norm.push_back(tgt.norm());
}

template <typename T>
int removeElements(std::vector<T> &vec, const std::vector<int> &mask) {
  int num = 0;
  auto it = vec.begin();
  for (auto m : mask) {
    if (m)
      ++it;
    else
      it = vec.erase(it);
    ++num;
  }
  return num;
}

std::string getTimeDir(std::string dataset_name, std::string opt = "") {
  auto now = std::chrono::system_clock::now();
  std::time_t now_c = std::chrono::system_clock::to_time_t(now);

  // Convert current time to string
  std::stringstream ss;
  ss << std::put_time(std::localtime(&now_c), "%m%d-%H%M");
  std::string time_str = ss.str();

  std::string dir_name = time_str + "_" + dataset_name + opt;
  fs::create_directory(dir_name);
  return dir_name;
}

void saveRes(eval &evals, std::string time_dir) {
  // Construct filename with current time
  std::string csvname = time_dir + "/" + evals.method_name + ".csv";
  std::ofstream file(csvname);
  file << "idx,  img1_idx,  img2_idx,  num_pts,  time,  t_err,  R_err, "
          "T_GTnorm, rgt_lie\n";
  evals.average_time /= evals.time.size();
  for (size_t i = 0; i < evals.R_err_per_round.size(); ++i) {
    file << i + 1 << "," << evals.img_pair[i].first << ","
         << evals.img_pair[i].second << "," << evals.num_pts[i] << ","
         << evals.time[i] << "," << evals.t_err_per_round[i] << ","
         << evals.R_err_per_round[i] << "," << evals.t_gt_norm[i] << ","
         << evals.R_gt[i].transpose();
    // if (evals.noise[i] != 0)
    //     file << "," << evals.noise[i];
    file << std::endl;
  }
  file.close();
}

string GetFileName(string path) {
  size_t pos = path.find_last_of('/');
  if (pos != std::string::npos)
    path.erase(0, pos + 1);
  pos = path.find_last_of('/');
  if (pos != std::string::npos)
    path.erase(0, pos + 1);
  return path;
}

class ResultAggregator {
public:
  void addResult(const string &solver_name, const ImagePair &pair,
                 const PoseResult &result) {
    eval &e = getOrCreateEval(solver_name);
    calcEval(unskew(pair.R_gt), pair.t_gt, e, to_string(pair.img1_id),
             to_string(pair.img2_id), result.t_error, result.r_error,
             pair.y_n.size(), result.time_taken);
  }

  void saveAllResults(const string &dir_name) {
    for (auto &[name, e] : evaluations) {
      saveRes(e, dir_name);
    }
    saveErrorSummary(dir_name);
  }

private:
  map<string, eval> evaluations;

  eval &getOrCreateEval(const string &name) {
    if (evaluations.find(name) == evaluations.end()) {
      evaluations.emplace(name, eval(name));
    }
    return evaluations[name];
  }

  void saveErrorSummary(const string &dir_name) {
    ofstream file(dir_name + "/errors.txt");
    file << "method,      avr_time,      R_err ,      t_err" << endl;

    for (const auto &[name, e] : evaluations) {
      file << setw(12) << name << ": " << setw(10) << e.average_time << setw(10)
           << e.total_R_Fn << setw(10) << e.total_t_cos << endl;
    }
    file.close();
  }
};

void evaluateImagePair(const ImagePair &pair,
                       vector<unique_ptr<PoseSolver>> &solvers,
                       ResultAggregator &results) {
  for (auto &solver : solvers) {
    PoseResult result;
    result.time_taken = TIME_IT(solver->GetPose(
        result.R_estimated, result.t_estimated, pair.y_n, pair.z_n,
        pair.y_cv_pix, pair.z_cv_pix, cameras[pair.img1_id]->intrinsic));

    result.r_error = (result.R_estimated - pair.R_gt).norm();
    result.t_error = (result.t_estimated - pair.t_gt).norm();

    results.addResult(solver->getName(), pair, result);
  }
}
