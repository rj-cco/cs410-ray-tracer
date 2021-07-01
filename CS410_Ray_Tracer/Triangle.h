#ifndef TRIANGLE_H
#define TRIANGLE_H
#include <Eigen/Dense>
#include <vector>
#include <string>

class Obj;
class Triangle {
private:
    std::string line;
    std::vector<Eigen::RowVector3d> vertices;
    double Ns, illumination;
    Eigen::RowVector3d Ka, Kd, Ks;
    mutable double beta, gamma, tval;
    double index;
    mutable Eigen::RowVector3d AavgN, BavgN, CavgN;
    mutable Eigen::RowVector3d dir_to_cam;
    mutable std::vector<std::vector<double>> vertex_map;
    mutable double cutoff_ang;

public:
    Triangle(std::string l, std::vector<Eigen::RowVector3d> verts, const double &spec_exp,
             Eigen::RowVector3d ambi, Eigen::RowVector3d diffuse, Eigen::RowVector3d spec,
             const double &illum, const double &i);
    std::vector<std::string> get_parts() const;
    int get_Ai() const;
    int get_Bi() const;
    int get_Ci() const;
    Eigen::RowVector3d get_vert_A() const;
    Eigen::RowVector3d get_vert_B() const;
    Eigen::RowVector3d get_vert_C() const;
    double get_index() const;
    /* Edit: FIX THIS BACK */
    Eigen::RowVector3d triangleNormal(const Eigen::RowVector3d &dir_cam) const;
    bool rayIntersect(const Eigen::RowVector3d &Lv, const Eigen::RowVector3d &Dv) const;
    double get_tval() const;
    double get_beta() const;
    double get_gamma() const;
    void set_cutoff(const double &ca) const;
    double get_cutoff() const;
    Eigen::RowVector3d get_ambient_mats() const;
    Eigen::RowVector3d get_diffuse_mats() const;
    Eigen::RowVector3d get_specular_mats() const;
    double get_spec_exponent() const;
    double get_illumination() const;
    void set_vtmap(std::vector<std::vector<double>> &vmap);
    Eigen::RowVector3d triangleMakeVrtAvergN (const std::string &label, const std::vector<double> &m, const Obj &obj) const;
    void triangleMakeAVrtNormal(const Obj &obj) const;
    void triangleMakeBVrtNormal(const Obj &obj) const;
    void triangleMakeCVrtNormal(const Obj &obj) const;
    Eigen::RowVector3d triangleSmoothedNormal(const double &bestb, const double &bestg) const;
};

#endif /* TRIANGLE_H */
