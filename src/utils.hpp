#include <vector>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <filesystem>

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
    std::vector<int> num_pts;
    double total_R_Fn;
    double total_t_cos;
    std::vector<double> noise;
    std::vector<std::pair<int, int>> img_pair;
    std::vector<double> time;
    double average_time;
    eval(std::string name) : method_name(name), total_R_Fn(0), total_t_cos(0) {}
};

void calcEval(eval &Method, int img_idx1, int img_idx2, double t_err, double r_err, int num, double time_cost, double noise = -1)
{
    Method.img_pair.push_back(std::make_pair(img_idx1, img_idx2));
    Method.t_err_per_round.push_back(1 - t_err);
    Method.R_err_per_round.push_back(r_err);
    Method.total_R_Fn += r_err;
    Method.total_t_cos += 1 - t_err;
    if (noise != -1)
        Method.noise.push_back(noise);
    Method.num_pts.push_back(num);
    Method.time.push_back(time_cost);
    Method.average_time += time_cost;
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

string getTimeDir(std::string dataset_name)
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    // Convert current time to string
    std::stringstream ss;
    ss << std::put_time(std::localtime(&now_c), "%m%d-%H%M");
    std::string time_str = ss.str();

    std::string dir_name = time_str + "_" + dataset_name;
    fs::create_directory(dir_name);
    return dir_name;
}

void saveRes(eval &evals, std::string time_dir)
{
    // Construct filename with current time
    std::string csvname = time_dir + "/" + evals.method_name + ".csv";
    std::ofstream file(csvname);
    file << "idx,img1_idx,img2_idx,num_pts,time,t_err,R_err\n";
    evals.average_time /= evals.time.size();
    for (int i = 0; i < evals.R_err_per_round.size(); ++i)
    {
        file << i << "," << evals.img_pair[i].first << "," << evals.img_pair[i].second << "," << evals.num_pts[i] << "," << evals.time[i] << "," << evals.t_err_per_round[i] << "," << evals.R_err_per_round[i] << endl;
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