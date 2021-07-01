#ifndef CAMERA_H_INCLUDED
#define CAMERA_H_INCLUDED
#include <Eigen/Dense>

class Camera {
private:
    double ex, ey, ez;
    double lx, ly, lz;
    double ux, uy, uz;
    double near;
    double left_b, right_b, bottom_b, top_b;
    double width, height;
public:
    /* Constructor for building the camera object */
    Camera(const double &e_x, const double &e_y, const double &e_z,
           const double &l_x, const double &l_y, const double &l_z,
           const double &u_x, const double &u_y, const double &u_z,
           const double &n,
           const double &l_b, const double &r_b, const double &b_b, const double &t_b,
           const double &w, const double &h);
    /* Creating the Eye, Look at Point, and Upv Vectors */
    Eigen::RowVector3d eye() const;
    Eigen::RowVector3d look() const;
    Eigen::RowVector3d up() const;
    /* Getters for accessing private variables */
    double get_near() const;
    double get_lb() const;
    double get_rb() const;
    double get_bb() const;
    double get_tb() const;
    double get_width() const;
    double get_height() const;
    /* Member Method Declarations */
    Eigen::RowVector3d setupW() const;
    Eigen::RowVector3d setupU() const;
    Eigen::RowVector3d setupV() const;
};

#endif /* CAMERA_H_INCLUDED */
