#pragma once
#include "Eigen/src/Core/Matrix.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <opencv2/core.hpp>
#include <vector>

using namespace Eigen;
using namespace cv;
using namespace std;

struct Trans_type {
  Matrix3d Rotation;
  Vector3d Translation;
};

class PoseSolver {
public:
  virtual ~PoseSolver() = default;

  virtual Trans_type GetPose(Matrix3d R_init, Vector3d t_init,
                             Matrix3d &E_out) = 0;
  virtual string getName() const = 0;
};
