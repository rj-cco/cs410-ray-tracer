#include "Light.h"

Light::Light(const double &lx, const double &ly, const double &lz, const double &lw,
             const double &l_e_red, const double &l_e_green, const double &l_e_blue) :
        x(lx), y(ly), z(lz), w(lw),
        l_e_r(l_e_red), l_e_g(l_e_green), l_e_b(l_e_blue) { }

Eigen::RowVector3d Light::light_coordinates() const {
    Eigen::RowVector3d location(x, y, z);
    return location;
}

Eigen::RowVector3d Light::light_color_values() const {
    Eigen::RowVector3d light(l_e_r, l_e_g, l_e_b);
    return light;
}
