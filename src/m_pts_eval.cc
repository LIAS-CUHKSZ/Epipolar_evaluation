
#include <cstring>
#include <iostream>
#include <cmath>
#include <string>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>

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

namespace fs = std::filesystem;
void appendRes(const std::string &res_path, const std::string &res_name, double t, double r, int pt_num, double time, double est_vars = 0);

int main(int argc, char **argv)
{
    if (argc != 6)
    {
        std::cout << "Usage: " << argv[0] << "  <datasetname>  <num_pts>  <pair_num1>  <pair_num2>  <sample_num>" << std::endl;
        return EXIT_FAILURE;
    }
    string dtset(argv[1]), pts_str(argv[2]), imgidx1(argv[3]), imgidx2(argv[4]), spn(argv[5]);
    if (imgidx1 == imgidx2)
    {
        std::cout << "Error: image indices must be different.Exit" << std::endl;
        return EXIT_FAILURE;
    }
    int pnum = stoi(pts_str), idx1 = stoi(imgidx1), idx2 = stoi(imgidx2), sample_num = stoi(spn);
    string cameras_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/cameras.txt";
    string images_txt_path = "../dataset/" + dtset + "/dslr_calibration_undistorted/images.txt";

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

    std::string res_path = "m-consistent-" + imgidx1 + "-" + imgidx2 + "-" + dtset;
    if (!fs::exists(res_path))
        fs::create_directory(res_path);

    vector<Vector3d> yptmp, zptmp, yntmp, zntmp;
    vector<Vector3d> y_pix, z_pix, y_n, z_n;
    vector<Point2d> y_cv_pix, z_cv_pix;
    vector<Point2d> ycvtmp, zcvtmp;
    yptmp.reserve(6000);
    zptmp.reserve(6000);
    yntmp.reserve(6000);
    zntmp.reserve(6000);
    y_pix.reserve(6000);
    z_pix.reserve(6000);
    y_n.reserve(6000);
    z_n.reserve(6000);
    y_cv_pix.reserve(6000);
    z_cv_pix.reserve(6000);
    ycvtmp.reserve(6000);
    zcvtmp.reserve(6000);

    for (auto yi : images[idx1]->observations) // traverse obs in first img
    {
        if (images[idx2]->observations.count(yi.first)) // if the point is visible in the another image
        {
            point_pair p1 = yi.second, p2 = images[idx2]->observations[yi.first];
            yptmp.emplace_back(p1.pixel);
            zptmp.emplace_back(p2.pixel);
            ycvtmp.emplace_back(Point2d(p1.pixel.x(), p1.pixel.y()));
            zcvtmp.emplace_back(Point2d(p2.pixel.x(), p2.pixel.y()));
            yntmp.emplace_back(p1.normalized);
            zntmp.emplace_back(p2.normalized);
        }
    }
    int covisible_num = yptmp.size();
    if (covisible_num < 200)
    {
        std::cout << "Error: too few covisible points.Exit" << std::endl;
        return EXIT_FAILURE;
    }
    if (covisible_num < pnum)
    {
        std::cout << "Warning: covisible points are less than the given target" << std::endl;
        std::cout << "All %d available pts would be used for eval" << covisible_num << std::endl;
        return EXIT_FAILURE;
    }

    /* ---------------------------remove outlier if necessary--------------------------------------*/
    // @note the low res datasets have much more outlier(mis-match points pair) than high-res one
    // you can add definition REMOVE_OUTLIER when using low-res one
    // This option is recommended when you want to evaluate the m-consistency
#ifdef REMOVE_OUTLIER
    vector<Vector3d> y_pix_in, z_pix_in, y_n_in, z_n_in; // inliers
    vector<Point2d> y_cv_pix_in, z_cv_pix_in;
    cv::Mat inlier_mask, intri_cv;
    cv::eigen2cv(cameras[images[img1]->camera_id]->intrinsic, intri_cv);
    cv::findEssentialMat(y_cv_pix, z_cv_pix, intri_cv, cv::LMEDS, 0.999, 20, inlier_mask);
    for (size_t i = 0; i < inlier_mask.rows; ++i) // remove outlier using inlier_mask
    {
        if (inlier_mask.at<uchar>(i, 0))
        {
            y_pix_in.emplace_back(yptmp[i]);
            z_pix_in.emplace_back(zptmp[i]);
            y_n_in.emplace_back(yntmp[i]);
            z_n_in.emplace_back(zntmp[i]);
            y_cv_pix_in.emplace_back(ycvtmp[i]);
            z_cv_pix_in.emplace_back(zcvtmp[i]);
        }
    }
    yptmp = std::move(y_pix_in);
    zptmp = std::move(z_pix_in);
    yntmp = std::move(y_n_in);
    zntmp = std::move(z_n_in);
    ycvtmp = std::move(y_cv_pix_in);
    zcvtmp = std::move(z_cv_pix_in);
    covisible_num = y_pix.size();
#endif // REMOVE_OUTLIER

    // compute the ground truth relative pose
    Matrix3d R1 = images[idx1]->global_T_image.linear().transpose().cast<double>();
    Matrix3d R_gt = ((images[idx1]->global_T_image.linear().transpose() * images[idx2]->global_T_image.linear()).cast<double>()).transpose();
    Vector3d t_gt = -R_gt * R1 * (images[idx2]->global_T_image.translation() - images[idx1]->global_T_image.translation()).cast<double>();
    if (t_gt.norm() < 0.075) // two-view geometry cannot evaluate the translation when it is too small
    {
        std::cout << "Error: translation too small.Exit" << std::endl;
        return EXIT_FAILURE;
    }
    t_gt.normalize(); // only the bearing vector is needed

    // temp variables to store res
    Matrix3d R_estimated;
    Vector3d t_estimated;
    double r_err_this_round[6] = {0}, t_err_this_round[6] = {0}, time_cons_this_round[6] = {0};
    string method_name[6] = {
        "c_est",
        "E_GN",
        "RANSAC-5pt",
        "LMEDS-5pt",
        "E_SDP",
        "E_EIGEN"};
    double est_vars = 0;

    for (auto iter = 1; iter <= sample_num; iter++)
    {
        // randomly choose num points,sample from the covisible points
        if (pnum == -1)
        {
            y_pix = std::move(yptmp);
            z_pix = std::move(zptmp);
            y_n = std::move(yntmp);
            z_n = std::move(zntmp);
            y_cv_pix = std::move(ycvtmp);
            z_cv_pix = std::move(zcvtmp);
            iter += sample_num; // only run once
        }
        else
        {
            // change rand seed
            std::srand(std::time(nullptr));
            y_pix.clear();
            z_pix.clear();
            y_n.clear();
            z_n.clear();
            y_cv_pix.clear();
            z_cv_pix.clear();
            for (auto i = 0; i < pnum; ++i)
            {
                int idx = std::rand() % covisible_num;
                y_pix.emplace_back(yptmp[idx]);
                z_pix.emplace_back(zptmp[idx]);
                y_n.emplace_back(yntmp[idx]);
                z_n.emplace_back(zntmp[idx]);
                y_cv_pix.emplace_back(ycvtmp[idx]);
                z_cv_pix.emplace_back(zcvtmp[idx]);
            }
        }

        /* ↓------------------consistent estimator------------------↓ */
        ConsistentEst est(cameras[images[idx1]->camera_id]->intrinsic);
        time_cons_this_round[0] += TIME_IT(
            est.GetPose(R_estimated, t_estimated, y_n, z_n, y_cv_pix, z_cv_pix););
        r_err_this_round[0] += (R_estimated - R_gt).norm();
        t_err_this_round[0] += (t_estimated - t_gt).norm();
        est_vars += est.var_est;

        /* ↑------------------consistent estimator------------------↑ */

        /* ↓------------------E Manifold GN------------------↓ */
        ManifoldGN MGN(cameras[images[idx1]->camera_id]->intrinsic);
        time_cons_this_round[1] += TIME_IT(MGN.GetPose(R_estimated, t_estimated, y_cv_pix, z_cv_pix, y_n, z_n););
        r_err_this_round[1] += (R_estimated - R_gt).norm();
        t_err_this_round[1] += (t_estimated - t_gt).norm();
        /* ↑------------------E Manifold GN------------------↑ */

        Mat E_cv, intrinsic_cv, R_cv, t_cv; // tmp vars

        /* ↓------------------RANSAC method------------------↓ */
        time_cons_this_round[2] += TIME_IT(
            eigen2cv(cameras[images[idx1]->camera_id]->intrinsic, intrinsic_cv);
            E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::RANSAC, 0.999, 1.0);
            recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
            cv2eigen(R_cv, R_estimated);
            cv2eigen(t_cv, t_estimated););
        r_err_this_round[2] += (R_estimated - R_gt).norm();
        t_err_this_round[2] += (t_estimated - t_gt).norm();

        double init_time = TIME_IT(eigen2cv(cameras[images[idx1]->camera_id]->intrinsic, intrinsic_cv);
                                   E_cv = findEssentialMat(y_cv_pix, z_cv_pix, intrinsic_cv, cv::LMEDS, 0.999, 1.0);
                                   recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                   cv2eigen(R_cv, R_estimated);
                                   cv2eigen(t_cv, t_estimated););
        time_cons_this_round[3] += init_time;
        r_err_this_round[3] += (R_estimated - R_gt).norm();
        t_err_this_round[3] += (t_estimated - t_gt).norm();

        /* ↑------------------RANSAC method------------------↑ */

        /* ↓------------------eigensolver estimator------------------↓ */
        // with the initial value given by RANSAC LMEDS
        eigenSolverWrapper gv_esv(y_n, z_n);
        time_cons_this_round[5] += TIME_IT(gv_esv.GetPose(R_estimated, t_estimated, R_estimated.transpose());) + init_time;
        r_err_this_round[5] += (R_estimated.transpose() - R_gt).norm();
        double tmp = (-R_estimated.transpose() * t_estimated - t_gt).norm();
        t_err_this_round[5] += tmp;
        /* ↑------------------eigensolver estimator------------------↑ */

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
        time_cons_this_round[4] += TIME_IT(npt_pose(npt_p2, npt_p1, C, y_n.size(), X_sol, sdp_E, true);
                                           eigen2cv(sdp_E, E_cv);
                                           cv::recoverPose(E_cv, y_cv_pix, z_cv_pix, intrinsic_cv, R_cv, t_cv);
                                           cv2eigen(R_cv, R_estimated);
                                           cv2eigen(t_cv, t_estimated););
        r_err_this_round[4] += (R_estimated - R_gt).norm();
        t_err_this_round[4] += (t_estimated - t_gt).norm();
        /* ↑------------------SDP on Essential Mat------------------↑ */
    }

    // save res to res_path dir
    if (pnum != -1)
    {
        for (size_t i = 0; i < 6; ++i)
        {
            r_err_this_round[i] /= sample_num;
            t_err_this_round[i] /= sample_num;
            time_cons_this_round[i] /= sample_num;
        }
        est_vars /= sample_num;
    }
    else
        pnum = covisible_num;
    for (size_t i = 1; i < 6; ++i)
    {
        appendRes(res_path, method_name[i], t_err_this_round[i], r_err_this_round[i], pnum, time_cons_this_round[i]);
    }
    appendRes(res_path, method_name[0], t_err_this_round[0], r_err_this_round[0], pnum, time_cons_this_round[0], est_vars);

    return EXIT_SUCCESS;
}

// 将结果写入csv文件
void appendRes(const std::string &res_path, const std::string &res_name, double t, double r, int pt_num, double time, double est_vars)
{
    std::string file_path = res_path + "/" + res_name + ".csv";
    bool file_exists = fs::exists(file_path);

    std::ofstream fout(file_path, std::ios::app);
    if (est_vars != 0)
    {
        if (!file_exists)
            fout << "pt_num,t,r,cost_time,est_vars" << std::endl;
        fout << pt_num << "," << t << "," << r << "," << time << "," << est_vars << std::endl;
    }
    else
    {
        if (!file_exists)
            fout << "pt_num,t,r,time_cost" << std::endl;
        fout << pt_num << "," << t << "," << r << "," << time << std::endl;
    }
    fout.close();
}
