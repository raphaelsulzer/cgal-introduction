#include <cgal_typedefs.h>



// for colmap fused.ply.vis(binary format)
bool LoadVisibilityFile(const char* pszVisFile, std::vector<std::vector<int>>& vis_info)
{
    FILE * fp = fopen(pszVisFile, "rb");
    if (fp == NULL){
        printf("Can not open the file %s\n", pszVisFile);
        return false;
    }

    uint64_t vis_num_points = 0;
    fread(& vis_num_points, sizeof(uint64_t), 1, fp);

    if (vis_num_points < 0){
        printf("non visibilty information\n");
        fclose(fp);

        return true;
    }

    vis_info.reserve(vis_num_points);
    for (int i = 0; i < vis_num_points; ++i){
        uint32_t num_visible_images = 0;
        fread(& num_visible_images, sizeof(uint32_t), 1, fp);

        std::vector<int> image_view;
        image_view.reserve(num_visible_images);
        for (int j = 0; j < num_visible_images; ++j){
            uint32_t image_idx = 0;
            fread(& image_idx, sizeof(uint32_t), 1, fp);
            image_view.push_back(image_idx);
        }
        vis_info.push_back(image_view);
    }

    fclose(fp);
    return true;
}

// because colmap save qvec and tvec
// we need cvec
Eigen::Vector4d NormalizeQuaternion(const Eigen::Vector4d& qvec){
  const double norm = qvec.norm();
  if (norm == 0) {
    // We do not just use (1, 0, 0, 0) because that is a constant and when used
    // for automatic differentiation that would lead to a zero derivative.
    return Eigen::Vector4d(1.0, qvec(1), qvec(2), qvec(3));
  } else {
    return qvec / norm;
  }
}

Eigen::Vector3d ProjectionCenterFromPose(const Eigen::Vector4d& qvec,
                                         const Eigen::Vector3d& tvec){
  // Inverse rotation as conjugate quaternion.
  const Eigen::Vector4d normalized_qvec = NormalizeQuaternion(qvec);
  const Eigen::Quaterniond quat(normalized_qvec(0), -normalized_qvec(1),
                                -normalized_qvec(2), -normalized_qvec(3));
  return quat * -tvec;
}

bool LoadImageFile(const char* pszImageFile, std::map<int, Point>& image_proj_center)
{
    FILE * fp = fopen(pszImageFile, "rb");
    if (fp == NULL){
        printf("Can not open the file %s\n", pszImageFile);
        return false;
    }

    uint64_t num_reg_images = 0;

    fread(& num_reg_images, sizeof(uint64_t), 1, fp);

    if (num_reg_images < 0){
        printf("non image information\n");
        fclose(fp);

        return true;
    }

    for (int i = 0; i < num_reg_images; ++i){
        uint32_t num_image_index = 0;
        fread(& num_image_index, sizeof(uint32_t), 1, fp);

        double qvec_array[4] = { 0 };
        fread(qvec_array, sizeof(double), 4, fp);

        double tvec_array[3] = { 0 };
        fread(tvec_array, sizeof(double), 3, fp);

        Eigen::Vector4d qvec;
        for (int j = 0; j < 4; ++j){
            qvec(j) = qvec_array[j];
        }

        Eigen::Vector3d tvec;
        for (int j = 0; j < 3; ++j){
            tvec(j) = tvec_array[j];
        }

        Eigen::Vector3d cvec = ProjectionCenterFromPose(qvec, tvec);

        Point proj_center(cvec(0), cvec(1), cvec(2));

        uint32_t num_camera_index = 0;
        fread(& num_camera_index, sizeof(uint32_t), 1, fp);

        // read char in binary
        std::string img_name;

        char name_char;
        do {
            fread(& name_char, sizeof(char), 1, fp);
            if (name_char != '\0') {
                img_name += name_char;
            }
        } while (name_char != '\0');

        uint64_t num_points2D = 0;
        fread(& num_points2D, sizeof(uint64_t), 1, fp);

        for (int j = 0; j < num_points2D; ++j){
            double xy[2] = { 0 };
            fread(xy, sizeof(double), 2, fp);

            uint64_t num_point3d_index = 0;
            fread(& num_point3d_index, sizeof(uint64_t), 1, fp);
        }

        image_proj_center.insert(std::map<int, Point>::value_type(num_image_index, proj_center));
    }
    fclose(fp);

    return true;
}
