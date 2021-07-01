#include "Obj.h"
#include "Triangle.h"
#include <fstream>
#include <iostream>
#include <utility>
#include <sstream>

Obj::Obj (std::string f, Eigen::Matrix4d trans_matrix, double cutoffA) : filename(std::move(f)), TRANS_MATRIX(std::move(trans_matrix)),
                                                                         cutoff(cutoffA) {}

std::vector<std::string> Obj::lines() const {
    std::vector<std::string> ls;
    std::string input;
    std::ifstream read;
    read.open(this->filename);
    if (read.fail()) { std::cerr << "ERROR! Could not open file /" << filename
                                 << '\n'; EXIT_FAILURE; }
    while (getline(read, input)) {
        ls.push_back(input);
    }
    read.close();
    return ls;
}

std::vector<Eigen::RowVector3d> Obj::vertices() const {
    std::vector<Eigen::RowVector3d> vs;
    for (const auto &line: this->lines()) {
        std::vector<std::string> tokens;
        std::istringstream iss(line);
        while(iss) {
            std::string temp;
            iss >> temp;
            tokens.push_back(temp);
        }
        if (tokens[0] == "v") {
            Eigen::Matrix<double, 4, 1> temp_vert(stod(tokens[1]), stod(tokens[2]), stod(tokens[3]), 1);
            temp_vert = TRANS_MATRIX * temp_vert;
            Eigen::RowVector3d vert(temp_vert(0,0), temp_vert(1,0), temp_vert(2,0));
            vs.push_back(vert);
        }
    }
    return vs;
}

std::vector<Triangle> Obj::triangles() const {
    double Ns = 0;
    Eigen::RowVector3d Ka;
    Eigen::RowVector3d Kd;
    Eigen::RowVector3d Ks;
    double illumination = 0;
    std::vector<Triangle> tris;
    std::string mtl_filename;
    int i = 0;
    for (const auto &line: this->lines()) {
        std::vector<std::string> tokens;
        std::istringstream iss(line);
        while(iss) {
            std::string temp;
            iss >> temp;
            tokens.push_back(temp);
        }
        if (tokens[0] == "mtllib") {
            mtl_filename = tokens[1];
        }
        else if (tokens[0] == "usemtl") {
            std::string material;
            std::string mtl_input;
            std::ifstream mtl_read;
            mtl_read.open(mtl_filename);
            if (mtl_read.fail()) { std::cerr << "ERROR! Could not open file /" << filename
                                             << '\n'; EXIT_FAILURE; }
            while (getline(mtl_read, mtl_input)) {
                std::vector<std::string> mtl_tokens;
                std::istringstream mtl_iss(mtl_input);
                while(mtl_iss) {
                    std::string mtl_temp;
                    mtl_iss >> mtl_temp;
                    mtl_tokens.push_back(mtl_temp);
                }
                if (mtl_tokens[0] == "newmtl") {
                    material = mtl_tokens[1];
                }
                if ((tokens[1] == material) and (mtl_tokens[0] == "Ns")) {
                    Ns = stod(mtl_tokens[1]);
                }
                else if ((tokens[1] == material) and (mtl_tokens[0] == "Ka")) {
                    Ka << stod(mtl_tokens[1]), stod(mtl_tokens[2]), stod(mtl_tokens[3]);
                }
                else if ((tokens[1] == material) and (mtl_tokens[0] == "Kd")) {
                    Kd << stod(mtl_tokens[1]), stod(mtl_tokens[2]), stod(mtl_tokens[3]);
                }
                else if ((tokens[1] == material) and (mtl_tokens[0] == "Ks")) {
                    Ks << stod(mtl_tokens[1]), stod(mtl_tokens[2]), stod(mtl_tokens[3]);
                }
                else if ((tokens[1] == material) and (mtl_tokens[0] == "illum")) {
                    illumination = stod(mtl_tokens[1]);
                }
            }
            mtl_read.close();
        }
        else if (tokens[0] == "f") {
            /*std::cout << "---------- Index: " << i << " ----------" << '\n';
            std::cout << "Ns: " << Ns << '\n';
            std::cout << "Ka: " << Ka << '\n';
            std::cout << "Kd: " << Kd << '\n';
            std::cout << "Ks: " << Ks << '\n';
            std::cout << "Illum: " << illumination << '\n';*/
            std::vector<Eigen::RowVector3d> verts = this->vertices();
            Triangle newTri(line, verts, Ns, Ka, Kd, Ks, illumination, i);
            tris.push_back(newTri);
            i += 1;
        }
    }
    return tris;
}

double Obj::get_cutoff() const {
    return cutoff;
}

std::vector<std::vector<double>> Obj::build_vertex_map() const {
    //std::cout << "How many times is this running" << '\n';
    std::vector<Triangle> triangles = this->triangles();
    std::vector<std::vector<double>> vtmap(this->vertices().size());
    for (const auto &triangle: triangles) {
        vtmap[triangle.get_Ai()].push_back(triangle.get_index());
        vtmap[triangle.get_Bi()].push_back(triangle.get_index());
        vtmap[triangle.get_Ci()].push_back(triangle.get_index());
    }
    return vtmap;
}