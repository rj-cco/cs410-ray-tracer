#include "Ray.h"
#include "Sphere.h"
//#include <iostream>

Ray::Ray(Eigen::RowVector3d &L, Eigen::RowVector3d &D, Eigen::Matrix<double, 1, 17> &init_matrix, Eigen::RowVector3d &init_rowvec) :
        ray_start_position(L), ray_direction(D.normalized()),
        best_t(std::numeric_limits<double>::infinity()), best_sph(init_matrix), best_pt(init_rowvec),
        /* Edit #4 */best_tri(Triangle("", std::vector<Eigen::RowVector3d>(), 0, Eigen::RowVector3d(0,0,0), Eigen::RowVector3d (0,0,0), Eigen::RowVector3d (0,0,0), 0, 0)),
        best_beta(0.0), best_gamma(0.0){ }

bool Ray::sphere_test(const Sphere &sph) const {
    Eigen::RowVector3d Tv(sph.center() - this->ray_start_position);
    double v = Tv.dot(this->ray_direction);
    double csq = Tv.dot(Tv);
    double disc = (sph.get_radius() * sph.get_radius()) - (csq - (v*v));
    if (disc > 0) {
        double t_val = v - sqrt(disc);
        if ((t_val < this->best_t) && (t_val > 0.00001)) {
            this->best_t = t_val;
            this->best_sph << sph.center(), sph.get_radius(), sph.ambient_materials(), sph.diffuse_materials(), sph.specular_materials(), sph.reflective_materials(), sph.get_refraction_index();
            this->best_pt << this->ray_start_position + t_val * this->ray_direction;
            /* Edit #3 */
            this->best_obj = "sphere";
        }
        return true;
    }
    else { return false; }
}

bool Ray::tri_test(const std::vector<Triangle> &tris) const {
    bool hit = false;
    for (const auto &t: tris) {
        bool res = t.rayIntersect(this->ray_start_position, this->ray_direction);
        if (res) {
            double t_val = t.get_tval();
            if ((t_val < this->best_t) && (t_val > 0.00001)) {
                this->best_t = t_val;
                this->best_beta = t.get_beta();
                this->best_gamma = t.get_gamma();
                this->best_pt << this->ray_start_position + t_val * this->ray_direction;
                this->best_tri = t;
                this->best_obj = "triangle";
                hit = true;
            }
        }
    }
    return hit;
}

Eigen::Matrix<double, 1, 17> Ray::get_best_sph() const {
    return best_sph;
}

Eigen::RowVector3d Ray::get_best_pt() const {
    return best_pt;
}

Eigen::RowVector3d Ray::get_L() const {
    return ray_start_position;
}

Eigen::RowVector3d Ray::get_D() const {
    return ray_direction;
}

Triangle Ray::get_best_tri() const {
    return best_tri;
}

std::string Ray::get_best_obj() const {
    return best_obj;
}

double Ray::get_best_t() const {
    return best_t;
}

double Ray::get_best_beta() const {
    return best_beta;
}
double Ray::get_best_gamma() const {
    return best_gamma;
}
