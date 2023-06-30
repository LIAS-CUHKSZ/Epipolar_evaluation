#pragma once
#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>

#include <opencv2/core.hpp>
#include <opencv2/calib3d.hpp>

#include <Eigen/Dense>

#define TIME_IT(code)                                                                                 \
    [&]() -> long long                                                                                \
    {                                                                                                 \
        auto start_time = std::chrono::high_resolution_clock::now();                                  \
        code;                                                                                         \
        auto end_time = std::chrono::high_resolution_clock::now();                                    \
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); \
        return duration.count();                                                                      \
    }()

namespace fs = std::filesystem;

struct eval
{
    std::string method_name;
    std::vector<double> t_err_per_round;
    std::vector<double> R_err_per_round;
    std::vector<double> t_gt_norm;
    std::vector<Eigen::Vector3d> R_gt;
    std::vector<int> num_pts;
    double total_R_Fn;
    double total_t_cos;
    std::vector<double> noise;
    std::vector<std::pair<std::string, std::string>> img_pair;
    std::vector<double> time;
    double average_time;
    eval(std::string name) : method_name(name), total_R_Fn(0), total_t_cos(0) {}
};

void calcEval(Eigen::Vector3d rgt, Eigen::Vector3d tgt, eval &Method, string img_idx1, string img_idx2, double t_err, double r_err, int num, double time_cost, double noise = -1)
{
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
int removeElements(std::vector<T> &vec, const std::vector<int> &mask)
{
    int num = 0;
    auto it = vec.begin();
    for (auto m : mask)
    {
        if (m)
            ++it;
        else
            it = vec.erase(it);
        ++num;
    }
    return num;
}

string getTimeDir(std::string dataset_name, std::string opt = "")
{
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

void saveRes(eval &evals, std::string time_dir)
{
    // Construct filename with current time
    std::string csvname = time_dir + "/" + evals.method_name + ".csv";
    std::ofstream file(csvname);
    file << "idx,  img1_idx,  img2_idx,  num_pts,  time,  t_err,  R_err, T_GTnorm, rgt_lie\n";
    evals.average_time /= evals.time.size();
    for (int i = 0; i < evals.R_err_per_round.size(); ++i)
    {
        file << i + 1 << "," << evals.img_pair[i].first << "," << evals.img_pair[i].second << "," << evals.num_pts[i] << "," << evals.time[i] << "," << evals.t_err_per_round[i] << "," << evals.R_err_per_round[i] << "," << evals.t_gt_norm[i] << "," << evals.R_gt[i].transpose() << std::endl;
    }
    file.close();
}

void printProg(int now, int total)
{
    int progress = now * 40 / total;

    if (now == total)
    {
        std::cout << "\r[";
        for (int j = 0; j < 40; ++j)
        {
            std::cout << "=";
        }
        std::cout << "] 100%" << std::endl;
        std::cout << "\n"
                  << std::flush << std::endl;
    }
    else
    {
        std::cout << "\r[";
        for (int j = 0; j < progress; ++j)
            std::cout << "=";
        for (int j = progress; j < 40; ++j)
            std::cout << " ";
        std::cout << "] " << now * 100.0 / total << "%" << std::flush;
    }
}

int DecomposeEssential(cv::Mat &_rotation, cv::Mat &_translation, cv::Mat &_essential, std::vector<cv::Point2d> pts1, std::vector<cv::Point2d> pts2)
{

    return 0;
}

string GetFileName(string path)
{
    size_t pos = path.find_last_of('/');
    if (pos != std::string::npos)
        path.erase(0, pos + 1);
    pos = path.find_last_of('/');
    if (pos != std::string::npos)
        path.erase(0, pos + 1);
    return path;
}

// std::cout << valid_round << "--------" << endl;
// std::cout << "[R_est]\n"
// 		  << R_estimated << endl;
// std::cout << "[GT_R]\n"
// 		  << R_gt << endl;
// std::cout << " [E_est]\n"
// 		  << sdp_E << endl;
// std::cout << "[E_GT]\n"
// 		  << E_ground << endl;

// std::cout << "[T_est] " << t_estimated.transpose() << endl;
// std::cout << "[GT_T] " << t_gt.transpose() << endl;
// std::cout << "[scale T]" << t_with_scale.transpose() << endl;
// std::cout << " [norm of T]" << t_with_scale.norm() << endl;