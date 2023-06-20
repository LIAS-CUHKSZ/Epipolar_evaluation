#include <vector>
#include <algorithm>

struct eval
{
	string method_name;
	vector<double> t_err_per_round;
	vector<double> R_err_per_round;
	double total_R_Fn;
	double total_t_cos;
	vector<double> noise;
	eval(string name) : method_name(name), total_R_Fn(0), total_t_cos(0) {}
};

void calcEval(eval &Method, double t_err, double r_err, double noise = -1)
{
	Method.t_err_per_round.push_back(t_err);
	Method.R_err_per_round.push_back(r_err);
	Method.total_R_Fn += r_err;
	Method.total_t_cos += 1 - t_err;
	if (noise != -1)
		Method.noise.push_back(noise);
}

template <typename T>
int removeElements(std::vector<T> &vec, const std::vector<int> &mask)
{
    int num=0;
    auto it = vec.begin();
    for (auto m : mask)
    {
        if (m)
        {
            ++it;
        }
        else
        {
            it = vec.erase(it);
            ++num;
        }
    }
    return num;
}

