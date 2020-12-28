#include "utils.h"
#include <math.h>
#include <stdio.h>

double GaussianKernel(double r, double sigma) {
    // sigma is chosen to be the mean distancebetween all pairs of colors in the original palette
    return exp(-(r * r) / (2 * sigma * sigma));
}

LAB RGB2LAB(RGB rgb) {
    double var_r = rgb.R / 255.0;
    double var_g = rgb.G / 255.0;
    double var_b = rgb.B / 255.0;

    if(var_r > 0.04045) var_r = pow((var_r + 0.055) / 1.055, 2.4);
    else var_r /= 12.92;
    if(var_g > 0.04045) var_g = pow((var_g + 0.055) / 1.055, 2.4);
    else var_g /= 12.92;
    if(var_b > 0.04045) var_b = pow((var_b + 0.055) / 1.055, 2.4);
    else var_b /= 12.92;

    var_r *= 100; var_g *= 100; var_b *= 100;

    double X = var_r * 0.4124 + var_g * 0.3576 + var_b * 0.1805;
    double Y = var_r * 0.2126 + var_g * 0.7152 + var_b * 0.0722;
    double Z = var_r * 0.0193 + var_g * 0.1192 + var_b * 0.9505;

    double var_x = X / 95.047;
    double var_y = Y / 100.0;
    double var_z = Z / 108.883;
    if(var_x > 0.008856) var_x = pow(var_x, 1.0 / 3.0);
    else var_x = 7.787 * var_x + 16.0 / 116.0;
    if(var_y > 0.008856) var_y = pow(var_y, 1.0 / 3.0);
    else var_y = 7.787 * var_y + 16.0 / 116.0;
    if(var_z > 0.008856) var_z = pow(var_z, 1.0 / 3.0);
    else var_z = 7.787 * var_z + 16.0 / 116.0;
    
    LAB lab(116.0 * var_y - 16.0, 500.0 * (var_x - var_y), 200.0 * (var_y - var_z));
    return lab;
}

RGB LAB2RGB(LAB lab) {
    double var_y = (lab.L + 16.0) / 116.0;
    double var_x = lab.A / 500.0 + var_y;
    double var_z = var_y - lab.B / 200.0;

    if(var_y > 0.206893034422) var_y = pow(var_y, 3.0);
    else var_y = (var_y - 16.0 / 116.0) / 7.787;
    if(var_x > 0.206893034422) var_x = pow(var_x, 3.0);
    else var_x = (var_x - 16.0 / 116.0) / 7.787;
    if(var_z > 0.206893034422) var_z = pow(var_z, 3.0);
    else var_z = (var_z - 16.0 / 116.0) / 7.787;

    double X = 95.047 * var_x;
    double Y = 100.0 * var_y;
    double Z = 108.883 * var_z;

    var_x = X / 100.0;
    var_y = Y / 100.0;
    var_z = Z / 100.0;

    double var_r = var_x * 3.2406 + var_y * -1.5372 + var_z * -0.4986;
    double var_g = var_x * -0.9689 + var_y * 1.8758 + var_z * 0.0415;
    double var_b = var_x * 0.0557 + var_y * -0.2040 + var_z * 1.0570;

    if(var_r > 0.0031308) var_r = 1.055 * pow(var_r, 1.0 / 2.4) - 0.055;
    else var_r *= 12.92;
    if(var_g > 0.0031308) var_g = 1.055 * pow(var_g, 1.0 / 2.4) - 0.055;
    else var_g *= 12.92;
    if(var_b > 0.0031308) var_b = 1.055 * pow(var_b, 1.0 / 2.4) - 0.055;
    else var_b *= 12.92;

    RGB rgb(var_r * 255.0, var_g * 255.0, var_b * 255.0);
    return rgb;
}

double determinant(double** src, int size) {
    // Gaussian Elimination
    if(size == 1) return src[0][0];
    #ifdef DEBUG
    for(int i = 0; i < size; i ++) {
        for(int j = 0; j < size; j ++) {
            if(j == size - 1) printf("%.3lf\n", src[i][j]);
            else printf("%.3lf ", src[i][j]);
        }
    }
    printf("\n");
    #endif
    for(int i = 0; i < size - 1; i ++) {
        for(int j = i + 1; j < size; j ++) {
            if(fabs(src[j][i] - 0.0) < 1e-4) continue;
            double ratio = src[j][i] / src[i][i];
            for(int k = i; k < size; k ++) src[j][k] -= src[i][k] * ratio;
        }
    }
    double ans = 1.0;
    for(int i = 0; i < size; i ++) ans *= src[i][i];
    return ans;
}

double adjugate(double** src, int size, int del_row, int del_col) {
    double** assist = new double*[size - 1];
    for(int i = 0; i < size - 1; i ++) assist[i] = new double[size - 1];
    for(int i = 0; i < size; i ++) {
        if(i == del_row) continue;
        for(int j = 0; j < size; j ++) {
            if(j == del_col) continue;
            int assist_i = i > del_row ? i - 1 : i;
            int assist_j = j > del_col ? j - 1 : j;
            assist[assist_i][assist_j] = src[i][j];
        }
    }
    double ans = determinant(assist, size - 1);
    for(int i = 0; i < size - 1; i ++) delete []assist[i];
    delete []assist;
    return ans;
}

void inverse(double** src, double** dst, int size) {
    double** adju = new double*[size];
    for(int i = 0; i < size; i ++) adju[i] = new double[size];
    for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++) adju[j][i] = adjugate(src, size, i, j) * pow(-1, i + j);
    #ifdef DEBUG
    for(int i = 0; i < size; i ++) {
        for(int j = 0; j < size; j ++) {
            if(j == size - 1) printf("%.3lf\n", adju[i][j]);
            else printf("%.3lf ", adju[i][j]);
        }
    }
    printf("\n");
    #endif
    double val = determinant(src, size);
    if(val - 0.0 < 1e-4) {
        printf("Error Occured in Inverse\n");
        return;
    }
    for(int i = 0; i < size; i ++)
    for(int j = 0; j < size; j ++) dst[i][j] = adju[i][j] / val;
    for(int i = 0; i < size; i ++) delete []adju[i];
    delete []adju;
}