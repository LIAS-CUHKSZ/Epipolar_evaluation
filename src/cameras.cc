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

#include "cameras.h"

#include <fstream>

bool ReadColmapCameras(const std::string& cameras_txt_path,
                       ColmapCameraPtrMap* cameras) {
  std::ifstream cameras_file_stream(cameras_txt_path, std::ios::in);
  if (!cameras_file_stream) {
    return false;
  }
  
  while (!cameras_file_stream.eof() && !cameras_file_stream.bad()) {
    std::string line;
    std::getline(cameras_file_stream, line);
    if (line.size() == 0 || line[0] == '#') {
      continue;
    }
    
    ColmapCamera* new_camera = new ColmapCamera();
    std::istringstream line_stream(line);
    line_stream >> new_camera->camera_id >> new_camera->model_name
                >> new_camera->width >> new_camera->height;
    while (!line_stream.eof() && !line_stream.bad()) {
      new_camera->parameters.emplace_back();
      line_stream >> new_camera->parameters.back();
    }

    new_camera->intrinsic << new_camera->parameters[0],0,new_camera->parameters[2],
                            0,new_camera->parameters[1], new_camera->parameters[3],
                            0,0,1;
    
    cameras->insert(std::make_pair(new_camera->camera_id,
                                   ColmapCameraPtr(new_camera)));
  }
  return true;
}
