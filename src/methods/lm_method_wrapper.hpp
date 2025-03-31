#pragma once
#include "opengv/relative_pose/methods.hpp"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"
#include <opengv/sac/Ransac.hpp>
#include <vector>
#include "Eigen/Dense"
#include "opengv/types.hpp"

using namespace opengv;
using namespace Eigen;
using namespace std;

class mlSolverWrapper
{
public:
    mlSolverWrapper(vector<Vector3d> &y1, vector<Vector3d> &z1)
    {
        y.resize(y1.size());
        z.resize(z1.size());
        for(size_t i = 0; i < y1.size(); i++)
        {
            y[i] = y1[i].normalized();
            z[i] = z1[i].normalized();
        }
    }

    void GetPose(Matrix3d &R_est, Vector3d &t_est, Matrix3d R_init, Vector3d T_init)
    {

        relative_pose::CentralRelativeAdapter adapter(y, z, T_init, R_init);
        transformation_t output;
        output = relative_pose::optimize_nonlinear(adapter);
    
        R_est = output.block(0, 0, 3, 3);
        t_est = output.block(0, 3, 3, 1).normalized();
    }

private:
    vector<Vector3d, Eigen::aligned_allocator<Vector3d>> y, z;
};