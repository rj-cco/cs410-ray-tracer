#include "Sphere.h"
#include <iostream>

Sphere::Sphere(const double &X, const double &Y, const double &Z, const double &R,
               const double &a_red, const double &a_green, const double &a_blue,
               const double &d_red, const double &d_green, const double &d_blue,
               const double &s_red, const double &s_green, const double &s_blue,
               const double &r_red, const double &r_green, const double &r_blue,
               const double &ref_i) :
        x(X), y(Y), z(Z), r(R),
        Ka_r(a_red), Ka_g(a_green), Ka_b(a_blue),
        Kd_r(d_red), Kd_g(d_green), Kd_b(d_blue),
        Ks_r(s_red), Ks_g(s_green), Ks_b(s_blue),
        Kr_r(r_red), Kr_g(r_green), Kr_b(r_blue),
        Ni(ref_i) { }

Eigen::RowVector3d Sphere::center() const {
    Eigen::RowVector3d center(x, y, z);
    return center;
}

double Sphere::get_radius() const {
    return r;
}

Eigen::RowVector3d Sphere::ambient_materials() const {
    Eigen::RowVector3d ambient_mats(Ka_r, Ka_g, Ka_b);
    return ambient_mats;
}

Eigen::RowVector3d Sphere::diffuse_materials() const {
    Eigen::RowVector3d diffuse_mats(Kd_r, Kd_g, Kd_b);
    return diffuse_mats;
}

Eigen::RowVector3d Sphere::specular_materials() const {
    Eigen::RowVector3d spec_mats(Ks_r, Ks_g, Ks_b);
    return spec_mats;
}

Eigen::RowVector3d Sphere::reflective_materials() const {
    Eigen::RowVector3d reflect_mats(Kr_r, Kr_g, Kr_b);
    return reflect_mats;
}

double Sphere::get_refraction_index() const {
    return Ni;
}

Eigen::RowVector3d Sphere::opacity_materials() const {
    Eigen::RowVector3d opacity_mats(1, 1, 1);
    Eigen::RowVector3d reflect_mats1(Kr_r, Kr_g, Kr_b);
    if (Ni == 0) {
        return opacity_mats;
    }
    else {
        return opacity_mats - reflect_mats1;
    }
}

Eigen::RowVector3d Sphere::refract_tray(const Eigen::RowVector3d &W, const Eigen::RowVector3d &pt, const Eigen::RowVector3d &N, const double &eta1, const double &eta2) {
    Eigen::RowVector3d T;
    double etar = eta1 / eta2;
    //std::cout << "Etar: " << etar << '\n';
    double a = etar * -1;
    //std::cout << "A: " << a << '\n';
    double wn = W.dot(N);
    //std::cout << "WN: " << wn << '\n';
    double radsq = (etar * etar) * ((wn*wn)-1) + 1;
    //std::cout << "Radsq: " << radsq << '\n';
    if (radsq < 0.0) {
        T << 0, 0, 0;
    }
    else {
        double b = (etar * wn) - sqrt(radsq);
        T = a * W + b * N;
    }
    return T;
}

std::vector<Eigen::RowVector3d> Sphere::refract_exit(const Eigen::RowVector3d &W, const Eigen::RowVector3d &pt, const double &eta_in, const double &eta_out) const {
    Eigen::RowVector3d Sn = pt - this->center();
    Sn.normalize();
    Eigen::RowVector3d T1 = this->refract_tray(W, pt, Sn, eta_out, eta_in);
    Eigen::Matrix<double, 1, 17> TEMP_INIT_MATRIX;
    TEMP_INIT_MATRIX << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    Eigen::RowVector3d TEMP_INIT_VEC(0, 0, 0);
    if (T1.sum() == 0.0) {
        Eigen::RowVector3d TEMP_L(0, 0, 0);
        Eigen::RowVector3d TEMP_D(0, 0, 0);
        //Ray refR(TEMP_L, TEMP_D, TEMP_INIT_MATRIX, TEMP_INIT_VEC);
        std::vector<Eigen::RowVector3d> temp_ray_holder;
        temp_ray_holder.push_back(TEMP_L);
        temp_ray_holder.push_back(TEMP_D);

        return temp_ray_holder;
    }
    else {
        Eigen::RowVector3d exit = pt + 2 * ((this->center() - pt).dot(T1)) * T1;
        Eigen::RowVector3d Nin = this->center() - exit;
        Nin.normalize();
        Eigen::RowVector3d T2 = this->refract_tray(-T1, exit, Nin, eta_in, eta_out);
        std::vector<Eigen::RowVector3d> temp_ray_holder;
        temp_ray_holder.push_back(exit);
        temp_ray_holder.push_back(T2);
        return temp_ray_holder;
    }
}

int Sphere::get_alpha() {
    return 16;
}
