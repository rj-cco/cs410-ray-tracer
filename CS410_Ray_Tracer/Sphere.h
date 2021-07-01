#ifndef SPHERE_H
#define SPHERE_H
// INCLUDING RAY.H MAY INTRODUCE CIRCULAR DEPENDENCY< KEEP AN EYE OUT
//#include "Ray.h"
#include <Eigen/Dense>
#include <vector>

class Sphere {
private:
    double x, y, z;
    double r;
    double Ka_r, Ka_g, Ka_b, Kd_r, Kd_g, Kd_b;
    double Ks_r, Ks_g, Ks_b, Kr_r, Kr_g, Kr_b;
    double Ni;
public:
    Sphere(const double &X, const double &Y, const double &Z, const double &R,
           const double &a_red, const double &a_green, const double &a_blue,
           const double &d_red, const double &d_green, const double &d_blue,
           const double &s_red, const double &s_green, const double &s_blue,
           const double &r_red, const double &r_green, const double &r_blue,
           const double &ref_index);
    /* Generate Sphere Center and Materials as Vectors */
    Eigen::RowVector3d center() const;
    double get_radius() const;
    Eigen::RowVector3d ambient_materials() const;
    Eigen::RowVector3d diffuse_materials() const;
    Eigen::RowVector3d specular_materials() const;
    Eigen::RowVector3d reflective_materials() const;
    double get_refraction_index() const;
    Eigen::RowVector3d opacity_materials() const;
    static Eigen::RowVector3d refract_tray(const Eigen::RowVector3d &W, const Eigen::RowVector3d &pt, const Eigen::RowVector3d &N, const double &eta1, const double &eta2) ;
    std::vector<Eigen::RowVector3d> refract_exit(const Eigen::RowVector3d &W, const Eigen::RowVector3d &pt, const double &eta_in, const double &eta_out) const;
    static int get_alpha() ;
};

#endif /* SPHERE_H */