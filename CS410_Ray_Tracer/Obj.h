#ifndef OBJ_H
#define OBJ_H
#include <Eigen/Dense>
#include <string>
#include <vector>

class Triangle;
class Obj {
private:
    std::string filename;
    Eigen::Matrix4d TRANS_MATRIX;
    double cutoff;
public:
    Obj (std::string file, Eigen::Matrix4d TM, double cA);
    std::vector<std::string> lines() const;
    std::vector<Eigen::RowVector3d> vertices() const;
    std::vector<Triangle> triangles() const;
    double get_cutoff() const;
    std::vector<std::vector<double>> build_vertex_map() const;
};

#endif /* OBJ_H */
