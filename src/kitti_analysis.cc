#include <cstring>
#include <iostream>
#include <cmath>
#include <opencv2/calib3d.hpp>
#include <opencv2/core.hpp>
#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include <yaml-cpp/yaml.h>
#include <iomanip> // For std::setw
#include <filesystem> // For creating directories

#include "draw_corre.hpp"
#include "pt3d.hpp"
#include "const_est.hpp"
#include "manifold_GN.hpp"
#include "npt_pose.h"
#include "eigensolver_wrapper.hpp"
#include "lm_wrapper.hpp"
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
vector<PointCorrespondence> loadPointCorrespondences(const string& matchesFile) {
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
vector<Matrix<double, 3, 4>> loadGroundTruthPoses(const string& posesFile) {
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
void pixelToNormalized(const vector<Point2d>& pixelPoints, vector<Vector3d>& normalizedPoints, const Matrix3d& K) {

    for (const auto& p : pixelPoints) {
        double x = (p.x-K(0, 2)) / K(0, 0);
        double y = (p.y-K(1, 2)) / K(1, 1);
        Vector3d normalized(x, y, 1.0);
        normalizedPoints.emplace_back(normalized);
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <sequence_number>  <eval_dir> " << std::endl;
        std::cout << "Example: " << argv[0] << " 00  eval" << std::endl;
        return EXIT_FAILURE;
    }

    // Load configuration using yaml-cpp
    YAML::Node config = YAML::LoadFile("../config.yaml");
    if (!config)
    {
        std::cerr << "Failed to open configuration file: ../config.yaml" << std::endl;
        return EXIT_FAILURE;
    }

    double t_min = config["config"]["kitti"]["t_min"].as<double>();
    double r_min = config["config"]["kitti"]["r_min"].as<double>();
    int pts_min = config["config"]["kitti"]["pts_min"].as<int>();
    bool use_lie = config["config"]["kitti"]["use_lie"].as<bool>();

    std::cout << "Configuration for kitti loaded:" << std::endl;
    std::cout << "t_min: " << t_min << ", r_min: " << r_min << ", pts_min: " << pts_min << ", use_lie: " << use_lie << std::endl;

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
    
    vector<PointCorrespondence> correspondences = loadPointCorrespondences(matchesFile);
    cout << "Loaded " << correspondences.size() << " point correspondences" << endl;
    vector<Matrix<double, 3, 4>> gtPoses = loadGroundTruthPoses(posesFile);
    cout << "Loaded " << gtPoses.size() << " ground truth poses" << endl;
    
    if (correspondences.size() != gtPoses.size()) {
        cout << "Warning: Number of point correspondences (" << correspondences.size() 
             << ") doesn't match number of poses (" << gtPoses.size() << ")" << endl;

    }
    // gtPoses.resize(50);
    // correspondences.resize(50);
    
    // create a folder under "build" named time(m-d-h-m)-dataset_name
    string dir_name = getTimeDir("kitti" + sequenceNum + eval_dir);
    
    // accumulate the error
    eval c_est("c_est"), e_m_gn("e_m_gn"), sdp("sdp"), pt5ransac("pt5ransac"), egsolver("eigenSolver"), lm("lmsolver");
    
    // the number of valid pairs
    int valid_round = 0;
    int progress = correspondences.size();
    int tmp_finish = 0; // for printing progress bar
    
    std::string result_dir = dir_name + "/result";
    std::filesystem::create_directory(result_dir);

    for (size_t i = 0; i < correspondences.size() && i < gtPoses.size(); ++i) {
        printProg(tmp_finish++, progress);
        
        const PointCorrespondence& corr = correspondences[i];
        
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
        ++valid_round;
        
        // temp variables to store res
        Matrix3d R_estimated;
        Vector3d t_estimated;
        double r_err_this_round, t_err_this_round, time_elapse;
        
        // get img pair name
        string img1path = to_string(corr.img1_idx),
               img2path = to_string(corr.img2_idx);
               
        // for debug img saving(name)
        std::ostringstream img1name, img2name;
        img1name << std::setw(6) << std::setfill('0') << corr.img1_idx << ".png";
        img2name << std::setw(6) << std::setfill('0') << corr.img2_idx << ".png";
        string img1showpath = seqDir + "/image_0/" + img1name.str();
        string img2showpath = seqDir + "/image_0/" + img2name.str();

                // temp vars
                Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars
                Matrix3d R_r5pt;
                Vector3d t_r5pt;
                eigen2cv(K, intrinsic_cv);

        /* ↓------------------RANSAC-5pt method------------------↓ */
        double ransac_time = TIME_IT(E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
                                     recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                     cv2eigen(R_cv, R_r5pt);
                                     cv2eigen(t_cv, t_r5pt););
        time_elapse = ransac_time;
        R_estimated = R_r5pt;
        t_estimated = t_r5pt;
        calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_r5pt, t_r5pt, use_lie); // Use calcErr
        calcEval(R_lie, t_with_scale, pt5ransac, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse);


        /* ↓------------------consistent estimator------------------↓ */
        ConsistentEst est(K);
        time_elapse = TIME_IT(est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix,1););
        calcErr(t_err_this_round, r_err_this_round, R_gt, t_gt, R_estimated, t_estimated, use_lie); // Use calcErr
        calcEval(R_lie, t_with_scale, c_est, img1path, img2path, t_err_this_round, r_err_this_round, total_covisible, time_elapse, est.var_est);
        if ((r_err_this_round > 0.02 || t_err_this_round > 0.01) && debug_img)
            DrawCorreImg(dir_name, img1showpath, img2showpath, y_cv_pix, z_cv_pix, valid_round, false);

        if (t_err_this_round > 0.3) {
            // Create folder for this pair
            std::ostringstream folder_name;
            folder_name << result_dir << "/" << std::setw(6) << std::setfill('0') << valid_round;
            std::filesystem::create_directory(folder_name.str());

            // Save image pair
            std::string img1_save_path = folder_name.str() +"/"+img1name.str();
            std::string img2_save_path = folder_name.str() +"/"+img2name.str();
            std::filesystem::copy(img1showpath, img1_save_path, std::filesystem::copy_options::overwrite_existing);
            std::filesystem::copy(img2showpath, img2_save_path, std::filesystem::copy_options::overwrite_existing);

            // Save ground truth and estimated poses
            std::ofstream pose_file(folder_name.str() + "/poses.txt");
            pose_file << std::fixed << std::setprecision(6);
            pose_file << "GT_R:\n" << R_gt << "\nGT_t:\n" << t_gt.transpose() << "\n";
            pose_file << "Estimated_R:\n" << R_estimated << "\nEstimated_t:\n" << t_estimated.transpose() << "\n";
            pose_file.close();

            // Save point pairs
            std::ofstream points_file(folder_name.str() + "/points.txt");
            points_file << std::fixed << std::setprecision(6);
            for (const auto& pt : y_cv_pix) {
                points_file << pt.x << " " << pt.y << " ";
            }
            points_file << "\n";
            for (const auto& pt : z_cv_pix) {
                points_file << pt.x << " " << pt.y << " ";
            }
            points_file.close();
        }
        /* ↑------------------consistent estimator------------------↑ */
    }
      

    // ------------------------------------ save results here----------------------------------------
    saveRes(c_est, dir_name);

    // save method error;
    std::ofstream file(dir_name + "/errors.txt");
    file << "method,      avr_time,      R_err ,      t_err" << endl;
    file << std::setw(12) << "c_est: " << std::setw(10) << c_est.average_time << std::setw(10) << c_est.total_R_Fn << std::setw(10) << c_est.total_t_cos << endl;


    std::cout << "----------------" << endl;
    std::cout << std::setw(15) << "[c_est] R:" << std::setw(12) << c_est.total_R_Fn << std::setw(12) << " t: " << c_est.total_t_cos << endl;

    return EXIT_SUCCESS;
}