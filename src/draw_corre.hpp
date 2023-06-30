#pragma once
#include "opencv2/highgui.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/imgproc.hpp"

#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <ctime>

using namespace std;
using namespace cv;

/**
 * @brief 从path1和path2中读取图片,显示在同一个窗口中,并将pix1和pix2的每一组对应点用不同颜色进行连线
 *
 * @param path1
 * @param path2
 * @param pix1
 * @param pix2
 */
void DrawCorreImg(std::string dts, std::string path1, std::string path2, std::vector<cv::Point2d> &pix1, std::vector<cv::Point2d> &pix2, int test_idx, bool show_now = false)
{
    static int cnt;
    cnt++;
    cv::Mat img1 = imread(path1);
    cv::Mat img2 = imread(path2);
    cv::Mat img;
    hconcat(img1, img2, img);
    for (int i = 0; i < pix1.size(); i++)
    {
        // 生成随机种子
        std::srand(std::time(nullptr));
        cv::Scalar color((std::rand() + i) % 256, (std::rand() + 100 + i) % 256, (std::rand() + 200 + 2 * i) % 256);
        cv::line(img, pix1[i], pix2[i] + cv::Point2d(img1.cols, 0), color, 2);
    }
    // if (show_now)
    // {
    //     cv::imshow("Correspondences", img);
    //     cv::waitKey(0);
    // }
    // 保存图片到路径下

    cv::imwrite(dts + "/" + std::to_string(test_idx) + ".jpg", img);
}