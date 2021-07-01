#ifndef RAY_H_INCLUDED
#define RAY_H_INCLUDED
#include "Sphere.h"
#include "Triangle.h"
#include <Eigen/Dense>

class Ray {
private:
    Eigen::RowVector3d ray_start_position;
    Eigen::RowVector3d ray_direction;
    mutable double best_t;
    mutable Eigen::Matrix<double, 1, 17> best_sph;
    /* Edit #4 */
    mutable Triangle best_tri;
    mutable Eigen::RowVector3d best_pt;
    /* Edit #3 */
    mutable std::string best_obj;
    mutable double best_beta, best_gamma;

public:
    Ray(Eigen::RowVector3d &L, Eigen::RowVector3d &D,
        Eigen::Matrix<double, 1, 17> &init_matrix, Eigen::RowVector3d &init_rowvec);
    Eigen::RowVector3d get_best_pt() const;
    Eigen::Matrix<double, 1, 17> get_best_sph() const;
    bool sphere_test(const Sphere &sph) const;
    /* Edit #1 */
    bool tri_test(const std::vector<Triangle> &tris) const;
    Eigen::RowVector3d get_L() const;
    Eigen::RowVector3d get_D() const;
    Triangle get_best_tri() const;
    std::string get_best_obj() const;
    double get_best_t() const;
    double get_best_beta() const;
    double get_best_gamma() const;
};

#endif /* RAY_H_INCLUDED */