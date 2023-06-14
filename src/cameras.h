// Copyright 2017 Thomas Schöps
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
#pragma once
#include <cstring>
#include <memory>
#include <unordered_map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

// Holds data of a COLMAP camera.
struct ColmapCamera {
  // Unique camera id.
  int camera_id;
  
  // Name of the distortion model. Determines the number of parameters.
  std::string model_name;
  
  // Image width in pixels.
  int width;
  
  // Image height in pixels.
  int height;
  
  // Distortion parameters. Their number and interpretation depends on the
  // distortion model.
  std::vector<double> parameters;

  Eigen::Matrix3d intrinsic;
};

typedef std::shared_ptr<ColmapCamera> ColmapCameraPtr;
typedef std::shared_ptr<const ColmapCamera> ColmapCameraConstPtr;

typedef std::vector<ColmapCameraPtr> ColmapCameraPtrVector;
typedef std::unordered_map<int, ColmapCameraPtr> ColmapCameraPtrMap;


// Loads ColmapCameraPtr from a COLMAP cameras.txt file and appends
// them to the cameras map (indexed by camera_id). Returns true if successful.
bool ReadColmapCameras(const std::string& cameras_txt_path,
                       ColmapCameraPtrMap* cameras);
