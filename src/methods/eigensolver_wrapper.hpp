#include "opengv/relative_pose/methods.hpp"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"
#include <opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp>
#include <opengv/sac/Ransac.hpp>
#include <memory>
#include <vector>
#include "Eigen/Dense"

using namespace opengv;
using namespace Eigen;
using namespace std;
class eigenSolverWrapper
{
public:
    eigenSolverWrapper(vector<Vector3d> &y1, vector<Vector3d> &z1)
    {
        y.resize(y1.size());
        z.resize(z1.size());
        for(size_t i = 0; i < y1.size(); i++)
        {
            y[i] = y1[i].normalized();
            z[i] = z1[i].normalized();
        }
    }

    void GetPose(Matrix3d &R_est, Vector3d &t_est, Matrix3d R_init)
    {

        relative_pose::CentralRelativeAdapter adapter(y, z, R_init);
        eigensolverOutput_t output;
        output.rotation = R_init;
        relative_pose::eigensolver(adapter, output);
        R_est = output.rotation;
        t_est = output.translation.normalized();
    }

private:
    vector<Vector3d, Eigen::aligned_allocator<Vector3d>> y, z;
};