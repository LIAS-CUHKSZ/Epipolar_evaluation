#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include "opengv/relative_pose/CentralRelativeAdapter.hpp"
#include "opengv/relative_pose/methods.hpp"
#include <memory>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp>
#include <vector>

using namespace opengv;
using namespace Eigen;
using namespace std;
class lmSolverWrapper {
public:
  lmSolverWrapper(vector<Vector3d> &y1, vector<Vector3d> &z1){
    y.resize(y1.size());
    z.resize(z1.size());
    for (size_t i = 0; i < y1.size(); i++) {
      y[i] = y1[i].normalized();
      z[i] = z1[i].normalized();
    }
  }

  void GetPose(Matrix3d &R_est, Vector3d &t_est, Matrix3d R_init,
               Vector3d t_init) {

    relative_pose::CentralRelativeAdapter adapter(y, z);
    adapter.sett12(-R_init.transpose()*t_init);
    adapter.setR12(R_init.transpose());
    auto res = relative_pose::optimize_nonlinear(adapter);
    R_est = res.block(0, 0, 3, 3).transpose();
    t_est = -R_est*res.block(0, 3, 3, 1).normalized();
    // auto res = relative_pose::sevenpt(adapter);
    // R_est = res.block(0, 0, 3, 3);
    // t_est = res.block(0, 3, 3, 1).normalized();
  }

private:
  vector<Vector3d, Eigen::aligned_allocator<Vector3d>> y;
  vector<Vector3d, Eigen::aligned_allocator<Vector3d>> z;
};