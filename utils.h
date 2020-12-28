#pragma once

class RGB {
public:
    double R;
    double G;
    double B;
    RGB(double R = 0.0, double G = 0.0, double B = 0.0) : R(R), G(G), B(B) {}
    void set(RGB rgb) {this->R = rgb.R; this->G = rgb.G; this->B = rgb.B;}
};

class LAB {
public:
    double L;
    double A;
    double B;
    LAB(double L = 0.0, double A = 0.0, double B = 0.0) : L(L), A(A), B(B) {}
    void set(LAB lab) {this->L = lab.L; this->A = lab.A; this->B = lab.B;}
};

LAB RGB2LAB(RGB rgb);
RGB LAB2RGB(LAB lab);

double GaussianKernel(double r, double sigma);

void inverse(double** src, double** dst, int size);