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

#include "images.h"
#include "Eigen/Dense"

#include <fstream>

bool ReadColmapImages(const std::string &images_txt_path,
					  bool read_observations,
					  ColmapImagePtrMap *images,
					  ColmapCameraPtrMap &cameras)
{
	std::ifstream images_file_stream(images_txt_path, std::ios::in);
	if (!images_file_stream)
	{
		return false;
	}
	while (!images_file_stream.eof() && !images_file_stream.bad())
	{
		std::string line;
		std::getline(images_file_stream, line);
		if (line.size() == 0 || line[0] == '#')
		{
			continue;
		}

		// Read image info line.
		ColmapImage *new_image = new ColmapImage();
		Eigen::Quaternionf image_R_global;
		std::istringstream image_stream(line);
		image_stream >>
			new_image->image_id >>
			image_R_global.w() >>
			image_R_global.x() >>
			image_R_global.y() >>
			image_R_global.z() >>
			new_image->image_T_global.translation()[0] >>
			new_image->image_T_global.translation()[1] >>
			new_image->image_T_global.translation()[2] >>
			new_image->camera_id >>
			new_image->file_path;
		new_image->image_T_global.linear() = image_R_global.toRotationMatrix();
		new_image->global_T_image = new_image->image_T_global.inverse();

		float du = cameras[new_image->camera_id]->parameters[2]; // pixel dist to x origin
		float dv = cameras[new_image->camera_id]->parameters[3]; // pixel dist to y origin

		// Read feature observations line.
		std::getline(images_file_stream, line);
		if (read_observations)
		{
			std::istringstream observations_stream(line);
			while (!observations_stream.eof() && !observations_stream.bad())
			{
				int pt3d_idx;
				Eigen::Vector2f pt2d;
				observations_stream >> pt2d.x() >> pt2d.y() >> pt3d_idx;
				if (pt3d_idx != -1) // if there exists a corresponding point
				{
					// transform pixel coordinates to image normalized coordinates
					pt2d.x() = (pt2d.x() - du) / cameras[new_image->camera_id]->parameters[0];
					pt2d.y() = (pt2d.y() - dv) / cameras[new_image->camera_id]->parameters[1];
					new_image->observations[pt3d_idx] = pt2d;
				}
			}
		}

		images->insert(std::make_pair(new_image->image_id, ColmapImagePtr(new_image)));
	}
	return true;
}