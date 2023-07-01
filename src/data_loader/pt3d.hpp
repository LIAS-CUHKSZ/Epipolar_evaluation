#pragma once
#include "Eigen/Dense"
#include <unordered_map>
#include <memory>
#include <vector>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;

struct ColmapPoint3D
{
    int point3D_id;

    Vector3d xyz;

    Vector3d rgb;

    double error; // reprojection error

    // key: image_id | value: feature index in corresponding image
    unordered_map<int, int> observations;
};

typedef shared_ptr<ColmapPoint3D> ColmapPoint3DPtr;
typedef shared_ptr<const ColmapPoint3D> ColmapPoint3DConstPtr;

typedef vector<ColmapPoint3DPtr> ColmapPoint3DPtrVector;
typedef unordered_map<int, ColmapPoint3DPtr> ColmapPoint3DPtrMap;

/**
 * @brief Read COLMAP point3D.txt file, the format is:
 *        point3D_id, X, Y, Z, R, G, B, error, n x [IMAGE_ID, POINT2D_IDX] (end with \n)
 *
 * @param images_txt_path path to point3D.txt
 * @param pt3ds ptr to map of point3D
 * @return true success
 * @return false fail
 */
bool ReadColmapPt3d(const std::string &pt3ds_txt_path, ColmapPoint3DPtrMap &pt3ds)
{
    ifstream pt3ds_txt(pt3ds_txt_path);
    if (!pt3ds_txt.is_open())
    {
        cerr << "Error: cannot open file " << pt3ds_txt_path << endl;
        return false;
    }

    string line;
    while (getline(pt3ds_txt, line))
    {
        if (line[0] == '#')
            continue;

        stringstream ss(line);
        ColmapPoint3DPtr pt3d = make_shared<ColmapPoint3D>();

        ss >> pt3d->point3D_id >> pt3d->xyz[0] >> pt3d->xyz[1] >> pt3d->xyz[2] >>
            pt3d->rgb[0] >> pt3d->rgb[1] >> pt3d->rgb[2] >> pt3d->error;

        int image_id, feature_idx;
        while (ss >> image_id >> feature_idx)
        {
            pt3d->observations[image_id] = feature_idx;
        }

        pt3ds[pt3d->point3D_id] = pt3d;
    }

    return true;
}
