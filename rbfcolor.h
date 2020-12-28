#pragma once

#include "utils.h"
#include <vector>
#include <string>

using std::vector, std::string;

class rbfcolor {
private:
    int ngrid;
    int res_int[8];
    double res_frac[8];
    int** weightindex;
    double** weightmap;

    // image
    vector<unsigned char> image;    // image content
    unsigned int width, height;     // image width & height
    unsigned int img_area;
    double* ori_r;
    double* ori_g;
    double* ori_b;
    double* res_r;
    double* res_g;
    double* res_b;

    // grid_lab
    unsigned int grid_size;
    LAB* grid_lab;
    // grid_rgb
    double* grid_R;
    double* grid_G;
    double* grid_B;
    
    // functions
    void calc_singlegrid(double ori_r, double ori_g, double ori_b);
    void prepare_grid();
    bool load_image(string filename);
    double calc_LAB_distance(double L1, double A1, double B1, double L2, double A2, double B2);
    LAB calc_single_point(double param, int matrixsize, double* oldpalette_L, double* oldpalette_A, double* oldpalette_B, double* diffpalette_L, double* diffpalette_A, double* diffpalette_B, LAB vsrc);
    void calc_grid_result(int palette_size, double* oldpalette_L, double* oldpalette_A, double* oldpalette_B, double* diffpalette_L, double* diffpalette_A, double* diffpalette_B);
public:
    rbfcolor(int ngrid, string image_path);
    void draw(RGB* input_palette, RGB* output_palette, int palette_size);
    void save_img(string name);
    void gridacc_kmeans(int center_num, RGB* res);
};