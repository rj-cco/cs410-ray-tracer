#include "Triangle.h"
#include "Obj.h"
#include <utility>
#include <vector>
#include <sstream>
#include <iostream>

Triangle::Triangle(std::string li, std::vector<Eigen::RowVector3d> v, const double &s_e,
                   Eigen::RowVector3d ambi, Eigen::RowVector3d diffuse, Eigen::RowVector3d spec,
                   const double &illum, const double &ind) :
        line(std::move(li)), vertices(std::move(v)), Ns(s_e),
        Ka(std::move(ambi)), Kd(std::move(diffuse)), Ks(std::move(spec)), illumination(illum),
        beta(0.0), gamma(0.0), tval(0.0), index(ind), cutoff_ang(0.0) { }

std::vector<std::string> Triangle::get_parts() const {
    std::string s = this->line;
    std::vector<std::string> tokens;
    std::stringstream iss(s);
    std::string temp;
    while (iss >> temp) { tokens.push_back(temp); }
    return tokens;
}

int Triangle::get_Ai() const {
    std::vector<std::string> parts = this->get_parts();
    int i = stoi(parts[1]);
    return i - 1;
}
int Triangle::get_Bi() const {
    std::vector<std::string> parts = this->get_parts();
    int i = stoi(parts[2]);
    return i - 1;
}
int Triangle::get_Ci() const {
    std::vector<std::string> parts = this->get_parts();
    int i = stoi(parts[3]);
    return i - 1;
}

Eigen::RowVector3d Triangle::get_vert_A() const {
    Eigen::RowVector3d vertA = vertices[this->get_Ai()];
    return vertA;
}

Eigen::RowVector3d Triangle::get_vert_B() const {
    Eigen::RowVector3d vertB = vertices[this->get_Bi()];
    return vertB;
}

Eigen::RowVector3d Triangle::get_vert_C() const {
    Eigen::RowVector3d vertC = vertices[this->get_Ci()];
    return vertC;
}

double Triangle::get_index() const{
    return index;
}

Eigen::RowVector3d Triangle::triangleNormal(const Eigen::RowVector3d &cam_dir) const {
    this->dir_to_cam = cam_dir;
    Eigen::RowVector3d vertA = this->get_vert_A();
    //std::cout << "VertA: " << vertA << '\n';
    Eigen::RowVector3d vertB = this->get_vert_B();
    //std::cout << "VertB: " << vertB << '\n';
    Eigen::RowVector3d vertC = this->get_vert_C();
    //std::cout << "VertC: " << vertC << '\n';
    Eigen::RowVector3d vectorBA = vertB - vertA;
    Eigen::RowVector3d vectorCB = vertC - vertB;
    Eigen::RowVector3d SurfaceNormal = vectorBA.cross(vectorCB);
    SurfaceNormal.normalize();
    if ((SurfaceNormal.dot(cam_dir)) < 0) {
        SurfaceNormal = -1 * SurfaceNormal;
    }
    return SurfaceNormal;
}

bool Triangle::rayIntersect(const Eigen::RowVector3d &L, const Eigen::RowVector3d &D) const {
    Eigen::RowVector3d Av = this->get_vert_A();
    Eigen::RowVector3d Bv = this->get_vert_B();
    Eigen::RowVector3d Cv = this->get_vert_C();
    const Eigen::RowVector3d& Lv = L;
    const Eigen::RowVector3d& Dv = D;
    Eigen::RowVector3d YV = Av - Lv;
    //std::cout << "YV: " << YV << '\n';
    Eigen::Matrix3d MM;
    MM << Av - Bv,
            Av - Cv,
            Dv;
    //std::cout << "MM: " << MM << '\n';
    Eigen::Matrix3d MMt = MM.transpose();
    //std::cout << "MMt: " << MMt << '\n';
    Eigen::Matrix3d MMi = MMt.inverse();
    //std::cout << "MMi: " << MMi << '\n';
    Eigen::RowVector3d row0 = MMi.row(0);
    Eigen::RowVector3d row1 = MMi.row(1);
    Eigen::RowVector3d row2 = MMi.row(2);
    double XV0 = YV.dot(row0), XV1 = YV.dot(row1), XV2 = YV.dot(row2);
    this->beta = XV0; this->gamma = XV1; this->tval = XV2;
    //std::cout << "Beta: " << this->beta << '\n';
    //std::cout << "Gamma: " << this->gamma << '\n';
    //std::cout << "Tval: " << this->tval << '\n';
    if ((beta > 0.0) and (gamma > 0.0) and (beta + gamma < 1.0) and (tval > 0.0001)) {
        return true;
    }
    else {
        return false;
    }
}

double Triangle::get_tval() const {
    return tval;
}

double Triangle::get_beta() const {
    return beta;
}

double Triangle::get_gamma() const {
    return gamma;
}

double Triangle::get_cutoff() const {
    return cutoff_ang;
}

void Triangle::set_cutoff(const double &cutoffangle) const {
    this->cutoff_ang = cutoffangle;
}

Eigen::RowVector3d Triangle::get_ambient_mats() const {
    return Ka;
}
Eigen::RowVector3d Triangle::get_diffuse_mats() const {
    return Kd;
}
Eigen::RowVector3d Triangle::get_specular_mats() const {
    return Ks;
}
double Triangle::get_spec_exponent() const {
    return Ns;
}
double Triangle::get_illumination() const {
    return illumination;
}

void Triangle::set_vtmap(std::vector<std::vector<double>> &vmap) {
    this->vertex_map = vmap;
}

void Triangle::triangleMakeAVrtNormal(const Obj &obj) const {
    //std::cout << "Vertex A" << '\n';
    Eigen::RowVector3d AvgN = this->triangleMakeVrtAvergN("A", vertex_map[this->get_Ai()], obj);
    AavgN = AvgN;
    //std::cout << "Calculated Avgerage normal A: " << AvgN << '\n';
    //std::cout << this->AavgN << '\n';

}

void Triangle::triangleMakeBVrtNormal(const Obj &obj) const {
    //std::cout << "Vertex B" << '\n';
    Eigen::RowVector3d BvgN = this->triangleMakeVrtAvergN("B", vertex_map[this->get_Bi()], obj);
    //std::cout << "Calculated Avgerage normal B: "<< BvgN << '\n';
    BavgN = BvgN;
}
void Triangle::triangleMakeCVrtNormal(const Obj &obj) const {
    //std::cout << "Vertex C" << '\n';
    Eigen::RowVector3d CvgN = this->triangleMakeVrtAvergN("C", vertex_map[this->get_Ci()], obj);
    //std::cout << "Calculated Avgerage normal C: " << CvgN << '\n';
    CavgN = CvgN;
}

Eigen::RowVector3d Triangle::triangleMakeVrtAvergN (const std::string &label, const std::vector<double> &vtm, const Obj &obj) const {
    //std::cout << "------ Building averaged vertex normal for triangle " << this->get_index() << " vertex " << label << " ------\n";
    Eigen::RowVector3d sumN(0,0,0);
    double count = 0.0;
    //for (const auto &v: vtm) {
    //    std::cout << v << '\n';
    //}
    for (auto const &i: vtm) {
        // What's inside this for loop will be a vector of double
        //std::cout << "Jesus Chroist" << '\n';
        //std::cout << "Triangle Normal: " << this->triangleNormal(this->dir_to_cam) << '\n';
        //std:: cout << "Parent Triangle Normal: " << obj.triangles()[i].triangleNormal(this->dir_to_cam) << '\n';
        //std::cout << "Before" << '\n';
        double costheta = this->triangleNormal(this->dir_to_cam).dot(obj.triangles()[i].triangleNormal(this->dir_to_cam));
        //std::cout << "Pass cosine calc" << '\n';
        //std::cout << "Cosine Theta: " << costheta << '\n';
        double theta = acos(costheta) * (180/3.141592653589793238);
        //std::cout << "Pass arc cosine calc" << '\n';
        //std::cout << "Theta (After Arccos): " << theta << '\n';
        //std::cout << "Cutoff: " << obj.get_cutoff() << '\n';
        if (theta < obj.get_cutoff()) {
            //std::cout << "EVER IN HERE??" << '\n';
            sumN += obj.triangles()[i].triangleNormal(this->dir_to_cam);
            count += 1.0;
        }
        //std::cout << "For Triangle "<< this->get_index() << " theta is: " << theta << " normal is: " << obj.triangles()[this->get_index()].triangleNormal(this->dir_to_cam) << " Sum is " << sumN << '\n';
    }
    Eigen::RowVector3d AvgN;
    double pd0 = sumN[0] / count, pd1 = sumN[1] /count, pd2 = sumN[2] /count;
    AvgN << pd0, pd1, pd2;
    //std::cout << "Averaged Normal is: " << AvgN << '\n';
    AvgN = -1 * AvgN;
    return AvgN;
}

Eigen::RowVector3d Triangle::triangleSmoothedNormal(const double &b_b, const double &b_g) const {
    //std::cout << "Inside Smoothed Normal" << '\n';
    //std::cout << AavgN << '\n';
    //std::cout << BavgN << '\n';
    //std::cout << CavgN << '\n';
    Eigen::RowVector3d sumN = (1.0 - b_b - b_g) * AavgN + b_b * BavgN + b_g * CavgN;
    sumN.normalize();
    return sumN;
}