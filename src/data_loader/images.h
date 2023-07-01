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

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/StdVector>

#include "cameras.h"

// Type for transformations in SE3. This example code sticks to Eigen, however
// the SE3f type from Sophus may be better suited, for example (since it uses
// quaternions for storing the rotation in contrast to the class chosen here,
// which uses a matrix for storage and allows all affine transformations for
// it in principle).
typedef Eigen::Transform<float, 3, Eigen::Affine> SE3f;

struct point_pair
{
  Eigen::Vector3d pixel;
  Eigen::Vector3d normalized;
};

// Holds data of a COLMAP image.
struct ColmapImage
{
  // Unique image id.
  int image_id;

  // Id of the camera model for this image.
  int camera_id;

  // Path to the image file, may be a relative path.
  std::string file_path;

  // Global-to-image transformation.
  SE3f image_T_global;

  // Image-to-global transformation.
  SE3f global_T_image;

  // 2D feature observations in this image.
  std::unordered_map<int, point_pair> observations;
};

typedef std::shared_ptr<ColmapImage> ColmapImagePtr;
typedef std::shared_ptr<const ColmapImage> ColmapImageConstPtr;

typedef std::vector<ColmapImagePtr> ColmapImagePtrVector;
typedef std::unordered_map<int, ColmapImagePtr> ColmapImagePtrMap;

// Loads ColmapImagePtr from a COLMAP images.txt file and appends them
// to the images map (indexed by image_id). Returns true if successful.
bool ReadColmapImages(const std::string &images_txt_path,
                      bool read_observations,
                      ColmapImagePtrMap *images,
                      ColmapCameraPtrMap& cameras);
