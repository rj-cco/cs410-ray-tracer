#ifndef LIGHT_H
#define LIGHT_H
#include <Eigen/Dense>

class Light {
private:
    double x, y, z, w, l_e_r, l_e_g, l_e_b;
public:
    Light(const double &lx, const double &ly, const double &lz, const double &lw,
          const double &l_e_red, const double &l_e_green, const double &l_e_blue);
    /* Generate Light Location (Point), and Color Emission */
    Eigen::RowVector3d light_coordinates() const;
    Eigen::RowVector3d light_color_values() const;
};

#endif /* LIGHT_H */
