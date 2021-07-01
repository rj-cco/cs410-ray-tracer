#include "Camera.h"

Camera::Camera(const double &e_x, const double &e_y, const double &e_z,
               const double &l_x, const double &l_y, const double &l_z,
               const double &u_x, const double &u_y, const double &u_z,
               const double &n,
               const double &l_b, const double &r_b, const double &b_b, const double &t_b,
               const double &w, const double &h) :
        ex(e_x), ey(e_y), ez(e_z),
        lx(l_x), ly(l_y), lz(l_z),
        ux(u_x), uy(u_y), uz(u_z),
        near(n),
        left_b(l_b), right_b(r_b), bottom_b(b_b), top_b(t_b),
        width(w), height(h) { }

Eigen::RowVector3d Camera::eye() const {
    Eigen::RowVector3d eye(ex, ey, ez);
    return eye;
}

Eigen::RowVector3d Camera::look() const {
    Eigen::RowVector3d look(lx, ly, lz);
    return look;
}

Eigen::RowVector3d Camera::up() const {
    Eigen::RowVector3d up(ux, uy, uz);
    return up;
}

double Camera::get_near() const {
    return near;
}

double Camera::get_lb() const {
    return left_b;
}

double Camera::get_rb() const {
    return right_b;
}

double Camera::get_bb() const {
    return bottom_b;
}

double Camera::get_tb() const {
    return top_b;
}

double Camera::get_width() const {
    return width;
}

double Camera::get_height() const {
    return height;
}

Eigen::RowVector3d Camera::setupW() const {
    Eigen::RowVector3d Wv(this->eye() - this->look());
    Wv.normalize();
    return Wv;
}

Eigen::RowVector3d Camera::setupU() const {
    Eigen::RowVector3d Uv = this->up().cross(this->setupW());
    Uv.normalize();
    return Uv;
}

Eigen::RowVector3d Camera::setupV() const {
    return this->setupW().cross(this->setupU());
}
