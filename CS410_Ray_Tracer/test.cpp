#include "Camera.h"
#include "Light.h"
#include "Sphere.h"
#include "Ray.h"
#include "Triangle.h"
#include "Obj.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <sstream>

#define PI 3.141592653589793238

using namespace std;
using namespace Eigen;

Ray pixel_ray(const int &i, const int &j, const Camera &c) {
    double px = (i/(c.get_width() - 1) * (c.get_rb() - c.get_lb()) + c.get_lb());
    double py = (j/(c.get_height() - 1) * (c.get_bb() - c.get_tb()) + c.get_tb());
    RowVector3d Lv = c.eye() + (c.get_near() * c.setupW()) + (px * c.setupU()) + (py * c.setupV());
    RowVector3d Dv = Lv - c.eye();
    Matrix<double, 1, 17> TEMP_INIT_MATRIX;
    TEMP_INIT_MATRIX << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    RowVector3d TEMP_INIT_VEC(0, 0, 0);
    Ray r(Lv, Dv, TEMP_INIT_MATRIX, TEMP_INIT_VEC);
    return r;
}

Matrix<double, 1, 17> ray_find(const Ray &ray, const vector<Sphere> &s_objs, const vector<Triangle> &t_objs) {
    for (const auto &s_obj : s_objs) {
        ray.sphere_test(s_obj);
    }
    ray.tri_test(t_objs);
    return ray.get_best_sph();
}

RowVector3d pairwise_product(const RowVector3d &am, const RowVector3d &mats) {
    RowVector3d pp;
    double pp0 = am[0] * mats[0], pp1 = am[1] * mats[1], pp2 = am[2] * mats[2];
    pp << pp0, pp1, pp2;
    return pp;
}

bool shadow(RowVector3d &pt, RowVector3d &lt, const vector<Sphere> &so, const vector<Triangle> &to) {
    Matrix<double, 1, 17> TEMP_INIT_MATRIX;
    TEMP_INIT_MATRIX << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    RowVector3d TEMP_INIT_VEC(0, 0, 0);
    RowVector3d L = lt - pt;
    Ray r(pt, L, TEMP_INIT_MATRIX, TEMP_INIT_VEC);
    double dtl = L.dot(r.get_D());
    for (const auto &sph: so) {
        if ((r.sphere_test(sph)) and (r.get_best_t() < dtl)) {
            return true;
        }
    }
    if ((r.tri_test(to)) and (r.get_best_t() < dtl)) { return true;}
    return false;
}

RowVector3d ray_trace(const Ray &r, RowVector3d &accum, RowVector3d &refatt, const vector<Sphere> &sph_obj, const vector<Light> &l_obj, const vector<Triangle> &tris, const RowVector3d &ambi, const double &ca, int level) {
    Matrix<double, 1, 17> TEMP_INIT_MATRIX;
    TEMP_INIT_MATRIX << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    RowVector3d TEMP_INIT_VEC(0, 0, 0);
    RowVector3d color;
    if ((ray_find(r, sph_obj, tris).sum() != 0) and (r.get_best_obj() == "sphere")) {
        Matrix<double, 1, 17> bs(r.get_best_sph());
        Sphere s(bs[0], bs[1], bs[2], bs[3], bs[4], bs[5], bs[6], bs[7], bs[8], bs[9], bs[10], bs[11], bs[12], bs[13], bs[14], bs[15], bs[16]);
        RowVector3d surface_normal(r.get_best_pt() - s.center());
        surface_normal.normalize();
        color = pairwise_product(ambi, s.ambient_materials());
        for (const auto &l : l_obj) {
            RowVector3d shadow_pt_test = r.get_best_pt();
            RowVector3d shadow_lt_test = l.light_coordinates();
            RowVector3d toL = l.light_coordinates() - r.get_best_pt();
            toL.normalize();
            RowVector3d temp_best_pt = r.get_best_pt();
            double Sn_dot_toL = surface_normal.dot(toL);
            bool sh = shadow(shadow_pt_test, shadow_lt_test, sph_obj, tris);
            if ((Sn_dot_toL > 0.0) and (!sh)) {
                color += pairwise_product(s.diffuse_materials(), l.light_color_values()) * Sn_dot_toL;
                RowVector3d toC = r.get_L() - r.get_best_pt();
                toC.normalize();
                RowVector3d spR = (2 * Sn_dot_toL * surface_normal) - toL;
                spR.normalize();
                double CdR = toC.dot(spR);
                if (CdR > 0.0) {
                    color += (pairwise_product(s.specular_materials(), l.light_color_values())) * pow(CdR, Sphere::get_alpha());
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            accum[i] += refatt[i] * s.opacity_materials()[i] * color[i];
        }
        if (level > 0) {
            //RowVector3d flec(0,0,0);
            RowVector3d Uinv = -1 * r.get_D();
            RowVector3d refR = (2 * surface_normal.dot(Uinv) * surface_normal) - Uinv;
            refR.normalize();
            RowVector3d temp_best_pt = r.get_best_pt();
            RowVector3d reflection_atten = pairwise_product(s.reflective_materials(), refatt);
            ray_trace(Ray(temp_best_pt, refR, TEMP_INIT_MATRIX, TEMP_INIT_VEC), accum, reflection_atten, sph_obj, l_obj, tris, ambi, ca, (level - 1));
            /*ray_trace(Ray(temp_best_pt, refR, TEMP_INIT_MATRIX, TEMP_INIT_VEC), flec, reflection_atten, sph_obj, l_obj, tris, ambi, (level - 1));
            for (int i = 0; i < 3; i++) {
                accum[i] += refatt[i] * s.opacity_materials()[i] * flec[i];
            }*/
        }
        if ((level > 0) and (s.get_refraction_index() != 0.0)) {
            //cout << "Entering Here" << '\n';
            RowVector3d thru(0,0,0);
            vector<RowVector3d> fractR_holder = s.refract_exit(-1*r.get_D(), r.get_best_pt(), s.get_refraction_index(), 1.0);
            Ray fractR(fractR_holder[0], fractR_holder[1], TEMP_INIT_MATRIX, TEMP_INIT_VEC);
            if ((fractR.get_L().sum() != 0) and (fractR.get_D().sum() != 0)) {
                RowVector3d temp_1_vec(1,1,1);
                RowVector3d refact = pairwise_product((temp_1_vec - s.reflective_materials()), refatt);
                //cout << refact << '\n';
                ray_trace(fractR, thru, refact, sph_obj, l_obj, tris, ambi, ca, (level - 1));
                for (int i = 0; i < 3; i++) {
                    //cout << "Refatt[i]: " << refatt[i] << '\n';
                    //cout << "1.0 - s.reflective_materials()[i]: " << (1.0 - s.reflective_materials()[i]) << '\n';
                    //cout << "thru[i]: " <<thru[i] << '\n';
                    accum[i] += refatt[i] * (1.0 - s.reflective_materials()[i]) * thru[i];
                    //cout << accum[i] << '\n';
                    //cout << accum << '\n';
                }
            }
        }
    }
    else if (r.get_best_obj() == "triangle") {
        Triangle tri = r.get_best_tri();
        Eigen::RowVector3d dir_to_cam = -1 * r.get_D();
        Eigen::RowVector3d SurfaceNormal;
        if (tri.get_cutoff() == 0.0) {
            SurfaceNormal = tri.triangleNormal(dir_to_cam);
        }
        else {
            SurfaceNormal = tri.triangleSmoothedNormal(r.get_best_beta(), r.get_best_gamma());
        }
        /*if (ca == 0) {
            SurfaceNormal = tri.triangleNormal(dir_to_cam);
        }
        else {
            SurfaceNormal = tri.triangleSmoothedNormal(r.get_best_beta(), r.get_best_gamma());
            SurfaceNormal = -1*(SurfaceNormal);
        }*/
        //cout << SurfaceNormal << '\n';
        color = pairwise_product(ambi, tri.get_ambient_mats());
        for (const auto &l : l_obj) {
            RowVector3d shadow_pt_test = r.get_best_pt();
            RowVector3d shadow_lt_test = l.light_coordinates();
            RowVector3d toL = l.light_coordinates() - r.get_best_pt();
            toL.normalize();
            //RowVector3d temp_best_pt = r.get_best_pt();
            double Sn_dot_toL = SurfaceNormal.dot(toL);
            bool sh = shadow(shadow_pt_test, shadow_lt_test, sph_obj, tris);
            if ((Sn_dot_toL > 0.0) and (!sh)) {
                color += pairwise_product(tri.get_diffuse_mats(), l.light_color_values()) * Sn_dot_toL;
                RowVector3d toC = r.get_L() - r.get_best_pt();
                toC.normalize();
                RowVector3d spR = (2 * Sn_dot_toL * SurfaceNormal) - toL;
                spR.normalize();
                double CdR = toC.dot(spR);
                if (CdR > 0.0) {
                    color += (pairwise_product(tri.get_specular_mats(), l.light_color_values())) * pow(CdR, tri.get_spec_exponent());
                }
            }
        }
        for (int i = 0; i < 3; i++) {
            accum[i] += refatt[i] * color[i];
        }
        if (level > 0) {
            RowVector3d Uinv = -1 * r.get_D();
            RowVector3d refR = (2 * SurfaceNormal.dot(Uinv) * SurfaceNormal) - Uinv;
            refR.normalize();
            RowVector3d temp_best_pt = r.get_best_pt();
            //RowVector3d reflection_atten = tri.get_illumination() * refatt;
            RowVector3d Kr;
            Kr << 0,0,0;
            if (tri.get_illumination() != 2) {
                Kr = tri.get_specular_mats();
            }
            RowVector3d reflection_atten = pairwise_product(Kr, refatt);
            ray_trace(Ray(temp_best_pt, refR, TEMP_INIT_MATRIX, TEMP_INIT_VEC), accum, reflection_atten, sph_obj, l_obj, tris, ambi, ca, (level - 1));
        }
    }
    //cout << "["<<accum <<"]\n";
    return accum;
}

RowVector3d round_pixels(const RowVector3d &pc) {
    double rp0 = round(pc.x()), rp1 = round(pc.y()), rp2 = round(pc.z());
    RowVector3d rp(rp0, rp1, rp2);
    return rp;
}

int main(int argc __attribute__((unused)), char *argv[]) {
    /* Object vectors to extract and hold data from file for construction */
    // Camera
    vector<double> camera_data;
    // Light Objects
    double ambient_red=0, ambient_green=0, ambient_blue=0;
    vector<Light> LightObjs;
    // Sphere Objects
    vector<Sphere> SphereObjs;
    // Getting the recursion depth
    int depth = 0;
    // Triangle Objects
    vector<Triangle> TriangleObjs;
    // Variable used for assigning cutoff
    double cutoff = 0;
    // Initialize the Transformation matrix with the Identity
    Matrix4d transformation_matrix;
    transformation_matrix << 1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1;

    string input;
    string filename = argv[1];
    ifstream read;
    read.open(filename);
    if (read.fail()) { cerr << "ERROR! Could not open file /" << argv[1]
                            << '\n'; return 0; }
    while (getline(read, input)) {
        vector<string> tokens;
        istringstream iss(input);
        while(iss) {
            string temp;
            iss >> temp;
            tokens.push_back(temp);
        }
        if (tokens[0] == "recursionlevel") {
            depth = stoi(tokens[1]);
        }
        else if (tokens[0] == "camera" || tokens[0] == "bounds" || tokens[0] == "res") {
            for (int i = 1; i < tokens.size()-1; i++) {
                double d = 0;
                stringstream(tokens[i]) >> d;
                camera_data.push_back(d);
            }
        }
        else if (tokens[0] == "ambient") {
            stringstream(tokens[1]) >> ambient_red;
            stringstream(tokens[2]) >> ambient_green;
            stringstream(tokens[3]) >> ambient_blue;
        }
        else if (tokens[0] == "light") {
            double ss_x, ss_y, ss_z, ss_w, ss_e_r, ss_e_g, ss_e_b;
            stringstream(tokens[1]) >> ss_x;
            stringstream(tokens[2]) >> ss_y;
            stringstream(tokens[3]) >> ss_z;
            stringstream(tokens[4]) >> ss_w;
            stringstream(tokens[5]) >> ss_e_r;
            stringstream(tokens[6]) >> ss_e_g;
            stringstream(tokens[7]) >> ss_e_b;
            LightObjs.__emplace_back(Light(ss_x, ss_y, ss_z, ss_w, ss_e_r, ss_e_g, ss_e_b));
        }
        else if (tokens[0] == "sphere") {
            vector<double> sphere_data;
            for (int i = 1; i < tokens.size(); i++) {
                double d = 0;
                stringstream(tokens[i]) >> d;
                sphere_data.push_back(d);
            }
            SphereObjs.__emplace_back(Sphere(sphere_data[0], sphere_data[1], sphere_data[2], sphere_data[3],
                                             sphere_data[4], sphere_data[5], sphere_data[6], sphere_data[7],
                                             sphere_data[8], sphere_data[9], sphere_data[10], sphere_data[11],
                                             sphere_data[12], sphere_data[13], sphere_data[14], sphere_data[15],
                                             sphere_data[16]));
        }
        else if (tokens[0] == "trans") {
            if (tokens[1] == "clear") {
                transformation_matrix << 1,0,0,0,
                        0,1,0,0,
                        0,0,1,0,
                        0,0,0,1;
            }
                /* Section for creation the Translation, Scale, and Rotation Matrices */
            else if (tokens[1] == "move") {
                double tx = stod(tokens[2]);
                double ty = stod(tokens[3]);
                double tz = stod(tokens[4]);
                Matrix4d translation_matrix;
                translation_matrix << 1,0,0,tx,
                        0,1,0,ty,
                        0,0,1,tz,
                        0,0,0,1;
                transformation_matrix = translation_matrix * transformation_matrix;
            }
            else if (tokens[1] == "scale") {
                double sx = stod(tokens[2]);
                double sy = stod(tokens[3]);
                double sz = stod(tokens[4]);
                Matrix4d scale_matrix;
                scale_matrix << sx,0,0,0,
                        0,sy,0,0,
                        0,0,sz,0,
                        0,0,0,1;
                transformation_matrix = scale_matrix * transformation_matrix;
            }
            else if (tokens[1] == "rota") {
                /* Implementation of Axis Angle Rotation */
                double wx = stoi(tokens[2]);
                double wy = stoi(tokens[3]);
                double wz = stoi(tokens[4]);
                double angle = stod(tokens[5]);
                /* 1. Calculate the normalization of the vector Wv */
                double magnitude = sqrt((wx*wx) + (wy*wy) + (wz*wz));
                wx = wx/magnitude, wy = wy/magnitude, wz = wz/magnitude;
                RowVector3d Wv(wx, wy, wz);
                /* 2. Rotate so theta becomes the z-axis */
                RowVector3d Mv(0,0,0);
                for (int i = 0; i < Wv.cols(); i++) {
                    Mv[i] = Wv[i];
                }
                int j;
                double min = 1000;
                for (int i = 0; i < 3; i ++) {
                    if (Mv[i] <= min) {
                        min = Mv[i];
                        j = i;
                    }
                }
                Mv[j] = 1.0;
                RowVector3d Uv(Wv.cross(Mv));
                double ux = Uv.x(), uy = Uv.y(), uz = Uv.z();
                double u_magnitude = sqrt((ux*ux) + (uy*uy) +(uz*uz));
                ux = ux/u_magnitude, uy = uy/u_magnitude, uz = uz/u_magnitude;
                Uv << ux,uy,uz;
                RowVector3d Vv(Wv.cross(Uv));
                Matrix3d R2z;
                R2z.row(0) = Uv, R2z.row(1) = Vv, R2z.row(2) = Wv;
                Matrix3d R2z_transpose(R2z.transpose());
                /* 3. Rotate by theta about the z-axis */
                double cosine = cos(angle * PI/ 180);
                double sine = sin(angle * PI/ 180);
                Matrix3d RMz;
                RMz << cosine, -sine, 0,
                        sine, cosine,  0,
                        0,    0,       1;
                Matrix3d RM_sub(R2z_transpose * RMz * R2z);
                /* 4. Create the final 4D rotation matrix */
                Matrix4d RM_final;
                RM_final << RM_sub(0,0), RM_sub(0,1), RM_sub(0,2), 0,
                        RM_sub(1,0), RM_sub(1,1), RM_sub(1,2), 0,
                        RM_sub(2,0), RM_sub(2,1), RM_sub(2,2), 0,
                        0                   , 0                   , 0                   , 1;
                /*. 5. Multiply the Rotation Matrix by the Transformation matrix */
                transformation_matrix = RM_final * transformation_matrix;
            }
        }
        else if (tokens[0] == "cutoffAngle") {
            cutoff = stod(tokens[1]);
        }
        else if (tokens[0] == "load") {
            Obj obj(tokens[1], transformation_matrix, cutoff);
            //cout << "Cutoff: " << obj.get_cutoff() << '\n';
            //cout << "----- Vertex Map -----\n";
            //cout << obj.build_vertex_map().size() << '\n';
            /*for (const auto &vrt: obj.build_vertex_map()) {
                for (const auto &i: vrt) {
                    cout << i << " ";
                }
                cout << '\n';
            }*/
            vector<vector<double>> vmap = obj.build_vertex_map();
            for (auto &triangle: obj.triangles()) {
                if (cutoff != 0.0) {
                    triangle.set_cutoff(cutoff);
                    triangle.set_vtmap(vmap);
                    triangle.triangleMakeAVrtNormal(obj);
                    triangle.triangleMakeBVrtNormal(obj);
                    triangle.triangleMakeCVrtNormal(obj);
                }
                triangle.set_cutoff(cutoff);
                TriangleObjs.push_back(triangle);
            }
            //tri.triangleMakeAVrtNormal(*this);
            //tri.triangleMakeBVrtNormal(*this);
            //tri.triangleMakeCVrtNormal(*this);
        }
    }
    read.close();
    Camera c(camera_data[0], camera_data[1], camera_data[2],
             camera_data[3], camera_data[4], camera_data[5],  camera_data[6],
             camera_data[7], camera_data[8], camera_data[9], camera_data[10],
             camera_data[11], camera_data[12], camera_data[13], camera_data[14],
             camera_data[15]);
    const RowVector3d ambient(ambient_red, ambient_green, ambient_blue);
    /* Ray Tracing Algorithm */
    /*int count = 1;
    for (const auto &tri: TriangleObjs) {
        cout << "--------- Testing Triangle " << count << " -----------\n";
        cout << "Ai: " << tri.get_Ai() << '\n';
        cout << "Bi: " << tri.get_Bi() << '\n';
        cout << "Ci: " << tri.get_Ci() << '\n';
        cout << "vertA: " << tri.get_vert_A() << '\n';
        cout << "vertB: " << tri.get_vert_B() << '\n';
        cout << "vertC: " << tri.get_vert_C() << '\n';
        cout << "Index: " << tri.get_index() << '\n';
        cout << "Smoothed Normal: " << tri.triangleSmoothedNormal() << '\n';
        count += 1;
    }
*/
    vector<RowVector3d> pixels;
    pixels.reserve(c.get_width() * c.get_height());
    /** ---------------- Note -----------------------
     *  Before we were being tested with the Non uniform resolution image, width was being passed in first
     *  and height was being passed in second. Keep an eye on this, the enumeration is the problem for almost
     *  everything. This current implementation works for square images and images where the width > height.
     *    - May need to see what happens when an image has height > width, mess around with enumeration
     */

    for (int i = 0; i < c.get_height(); i++) {
        for (int j = 0; j < c.get_width(); j++) {
            Ray r = pixel_ray(j, i, c);
            RowVector3d rgb(0, 0, 0);
            RowVector3d reflect(1, 1, 1);
            RowVector3d pix_color = ray_trace(r, rgb, reflect, SphereObjs, LightObjs, TriangleObjs, ambient, cutoff, depth);
            if (pix_color.sum() != 0) {
                pix_color = pix_color * 255;
                double r0 = pix_color[0], g1 = pix_color[1], b2 = pix_color[2];
                for (int p = 0; p < pix_color.size(); p++) {
                    if (pix_color[p] < 0 && p == 0) { r0 = 0; }
                    else if (pix_color[p] > 255 && p == 0) { r0 = 255; }
                    else if (pix_color[p] < 0 && p == 1) { g1 = 0; }
                    else if (pix_color[p] > 255 && p == 1) { g1 = 255; }
                    else if (pix_color[p] < 0 && p == 2) { b2 = 0; }
                    else if (pix_color[p] > 255 && p == 2) { b2 = 255; }
                }
                RowVector3d final_pix(r0, g1, b2);
                pix_color = round_pixels(final_pix);
                pixels.push_back(pix_color);
            }
            else { pixels.push_back(pix_color); }
        }
    }

    ofstream write_ppm(argv[2]);
    write_ppm << "P3\n";
    write_ppm << to_string(int(c.get_width())) + " " + to_string(int(c.get_height())) + " " + "255\n";
    for (const auto &pix : pixels) {
        ostringstream ss;
        ss << pix.x() << " " << pix.y() << " " << pix.z() << " ";
        write_ppm << ss.str();
        ss.clear();
    }
    write_ppm << '\n';
    write_ppm.close();

    return 0;
}
