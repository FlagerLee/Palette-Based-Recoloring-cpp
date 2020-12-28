#include "rbfcolor.h"
#include "png/lodepng.h"
#include <math.h>
#include <algorithm>

const int RBF_param_coff = 5;

bool big_change(LAB dir) {
    return fabs(dir.L) > 0.5 || fabs(dir.A) > 0.5 || fabs(dir.B) > 0.5;
}
bool out_boundary(RGB testrgb) {
    double out_threshold = 0.5;
    return testrgb.R < -out_threshold || testrgb.R > 255.0 + out_threshold ||
           testrgb.G < -out_threshold || testrgb.G > 255.0 + out_threshold ||
           testrgb.B < -out_threshold || testrgb.B > 255.0 + out_threshold;
}
double find_boundary(LAB vsrc, LAB dir, double l, double r) {
    double mid;
    for(int iter = 0; iter < 15; iter ++) {
        mid = 0.5 * (l + r);
        LAB lab(vsrc.L + mid * dir.L, vsrc.A + mid * dir.A, vsrc.B + dir.B);
        RGB testrgb = LAB2RGB(lab);
        
        if(out_boundary(testrgb)) r = mid;
        else l = mid;
    }
    return l;
}

rbfcolor::rbfcolor(int ngrid, string image_path) {
    this->ngrid = ngrid;
    this->grid_size = this->ngrid * this->ngrid * this->ngrid;

    if(!load_image(image_path)) exit(0);
}

void rbfcolor::calc_singlegrid(double ori_r, double ori_g, double ori_b) {
    int ntmp = this->ngrid + 1;
    int ntmpsqr = ntmp * ntmp;
    double tmpx_frac, tmpy_frac, tmpz_frac;
    int tmpx_int, tmpy_int, tmpz_int;

    tmpx_frac = ori_r / 255.0 * (double)this->ngrid;
    tmpx_int = floor(tmpx_frac);
    tmpx_frac = tmpx_frac - tmpx_int;
    if(tmpx_int == this->ngrid) {
        tmpx_int = this->ngrid - 1;
        tmpx_frac = 1.0;
    }

    tmpy_frac = ori_g / 255.0 * (double)this->ngrid;
    tmpy_int = floor(tmpy_frac);
    tmpy_frac = tmpy_frac - tmpy_int;
    if(tmpy_int == this->ngrid) {
        tmpy_int = this->ngrid - 1;
        tmpy_frac = 1.0;
    }

    tmpz_frac = ori_b / 255.0 * (double)this->ngrid;
    tmpz_int = floor(tmpz_frac);
    tmpz_frac = tmpz_frac - tmpz_int;
    if(tmpz_int == this->ngrid) {
        tmpz_int = this->ngrid - 1;
        tmpz_frac = 1.0;
    }

    int corner_ind = tmpz_int * ntmpsqr + tmpy_int * ntmp + tmpx_int;

    this->res_int[0] = corner_ind;
    this->res_int[1] = corner_ind + ntmpsqr;
    this->res_int[2] = corner_ind + ntmp;
    this->res_int[3] = corner_ind + ntmpsqr + ntmp;
    this->res_int[4] = corner_ind + 1;
    this->res_int[5] = corner_ind + ntmpsqr + 1;
    this->res_int[6] = corner_ind + ntmp + 1;
    this->res_int[7] = corner_ind + ntmpsqr + ntmp + 1;

    this->res_frac[0] = (1.0 - tmpx_frac) * (1.0 - tmpy_frac) * (1.0 - tmpz_frac);
    this->res_frac[1] = (1.0 - tmpx_frac) * (1.0 - tmpy_frac) * tmpz_frac;
    this->res_frac[2] = (1.0 - tmpx_frac) * tmpy_frac * (1.0 - tmpz_frac);
    this->res_frac[3] = (1.0 - tmpx_frac) * tmpy_frac * tmpz_frac;
    this->res_frac[4] = tmpx_frac * (1.0 - tmpy_frac) * (1.0 - tmpz_frac);
    this->res_frac[5] = tmpx_frac * (1.0 - tmpy_frac) * tmpz_frac;
    this->res_frac[6] = tmpx_frac * tmpy_frac * (1.0 - tmpz_frac);
    this->res_frac[7] = tmpx_frac * tmpy_frac * tmpz_frac;
}

void rbfcolor::prepare_grid() {
    this->grid_size = pow(ngrid + 1, 3);
    this->grid_lab = new LAB[this->grid_size];

    double step = 255.0 / ngrid;
    int tot = 0;
    for(int i = 0; i < this->ngrid + 1; i ++)
        for(int j = 0; j < this->ngrid + 1; j ++)
            for(int k = 0; k < this->ngrid + 1; k ++) {
                RGB rgb(k * step, j * step, i * step);
                this->grid_lab[tot] = RGB2LAB(rgb);
                tot ++;
            }

    this->img_area = this->width * this->height;
    this->weightmap = new double*[this->img_area];
    this->weightindex = new int*[this->img_area];

    for(int i = 0; i < this->img_area; i ++) {
        this->calc_singlegrid(this->ori_r[i], this->ori_g[i], this->ori_b[i]);

        this->weightindex[i] = new int[8];
        for(int j = 0; j < 8; j ++) this->weightindex[i][j] = res_int[j];
        this->weightmap[i] = new double[8];
        for(int j = 0; j < 8; j ++) this->weightmap[i][j] = res_frac[j];
    }

    this->grid_R = new double[this->grid_size];
    this->grid_G = new double[this->grid_size];
    this->grid_B = new double[this->grid_size];
}

bool rbfcolor::load_image(const string filename) {
    unsigned int error = lodepng::decode(this->image, this->width, this->height, filename);
    if(error != 0)
    {
        printf("Decode error %d: ", error);
        printf("%s\n", lodepng_error_text(error));
        return false;
    }
    if(this->width * this->height * 4 != this->image.size()) return false;
    int size = this->width * this->height;
    this->ori_r = new double[size];
    this->ori_g = new double[size];
    this->ori_b = new double[size];
    for(int i = 0; i < size; i ++) {
        this->ori_r[i] = this->image[i * 4];
        this->ori_g[i] = this->image[i * 4 + 1];
        this->ori_b[i] = this->image[i * 4 + 2];
    }
    this->prepare_grid();

    this->res_r = new double[size];
    this->res_g = new double[size];
    this->res_b = new double[size];
    return true;
}

double rbfcolor::calc_LAB_distance(double L1, double A1, double B1, double L2, double A2, double B2) {
    double K1 = 0.045, K2 = 0.015;
    double del_L = L1 - L2;
    double C1 = sqrt(A1 * A1 + B1 * B1);
    double C2 = sqrt(A2 * A2 + B2 * B2);
    double C_AB = C1 - C2;
    double H_AB = (A1 - A2) * (A1 - A2) + (B1 - B2) * (B1 - B2) - C_AB * C_AB;
    return del_L * del_L + C_AB * C_AB / (1.0 + K1 * C1) / (1.0 + K1 * C1) + H_AB / (1.0 + K2 * C1) / (1.0 + K2 * C1);
}

LAB rbfcolor::calc_single_point(double param, int matrixsize, double* oldpalette_L, double* oldpalette_A, double* oldpalette_B, double* diffpalette_L, double* diffpalette_A, double* diffpalette_B, LAB vsrc) {
    double** tmpMat = new double*[matrixsize];
    for(int i = 0; i < matrixsize; i ++) tmpMat[i] = new double[matrixsize];
    for(int u = 0; u < matrixsize; u ++)
        for(int v = 0; v < matrixsize; v ++) {
            double r = this->calc_LAB_distance(oldpalette_L[u], oldpalette_A[u], oldpalette_B[u], oldpalette_L[v], oldpalette_A[v], oldpalette_B[v]);
            tmpMat[u][v] = exp(-r * param);
        }
    
    double* tmpD = new double[matrixsize];
    for(int u = 0; u < matrixsize; u ++) {
        double r = this->calc_LAB_distance(oldpalette_L[u], oldpalette_A[u], oldpalette_B[u], vsrc.L, vsrc.A, vsrc.B);
        tmpD[u] = exp(-r * param);
    }
    #ifdef DEBUG
    for(int i = 0; i < matrixsize; i ++) {
        for(int j = 0; j < matrixsize; j ++) {
            if(j == matrixsize - 1) printf("%.3lf\n", tmpMat[i][j]);
            else printf("%.3lf ", tmpMat[i][j]);
        }
    }
    printf("\n");
    #endif
    inverse(tmpMat, tmpMat, matrixsize);
    #ifdef DEBUG
    for(int i = 0; i < matrixsize; i ++) {
        for(int j = 0; j < matrixsize; j ++) {
            if(j == matrixsize - 1) printf("%.3lf\n", tmpMat[i][j]);
            else printf("%.3lf ", tmpMat[i][j]);
        }
    }
    printf("\n");
    #endif

    double* precompute_inv = new double[matrixsize];
    for(int i = 0; i < matrixsize; i ++) precompute_inv[i] = 0.0;
    for(int i = 0; i < matrixsize; i ++)
    for(int j = 0; j < matrixsize; j ++)
    precompute_inv[i] += tmpMat[i][j] * tmpD[j];

    double delta_L = 0.0;
    double delta_A = 0.0;
    double delta_B = 0.0;

    double scale = 0.0;
    for(int i = 0; i < matrixsize; i ++) {
        scale += std::max(precompute_inv[i], 0.0);
    }

    for(int i = 0; i < matrixsize; i ++) if(precompute_inv[i] > 0.0) {
        delta_L += precompute_inv[i] / scale * diffpalette_L[i];
        delta_A += precompute_inv[i] / scale * diffpalette_A[i];
        delta_B += precompute_inv[i] / scale * diffpalette_B[i];
    }

    LAB lab(vsrc.L + delta_L, vsrc.A + delta_A, vsrc.B + delta_B);

    for(int i = 0; i < matrixsize; i ++) delete []tmpMat[i];
    delete []tmpMat;
    delete []tmpD;
    delete []precompute_inv;
    return lab;
}

void rbfcolor::calc_grid_result(int palette_size, double* oldpalette_L, double* oldpalette_A, double* oldpalette_B, double* diffpalette_L, double* diffpalette_A, double* diffpalette_B) {
    int tot = 0;
    double totr = 0.0;
    for(int u = 0; u < palette_size; u ++)
    for(int v = u + 1; v < palette_size; v ++) {
        double r = calc_LAB_distance(oldpalette_L[u], oldpalette_A[u], oldpalette_B[u], oldpalette_L[v], oldpalette_A[v], oldpalette_B[v]);

        tot = tot + 1;
        totr += sqrt(r);
    }
    if(palette_size > 1) totr /= (double)tot;
    else totr = 1.0;

    double param = (double)RBF_param_coff / (totr * totr);

    for(int i = 0; i < this->grid_size; i ++) {
        LAB vsrc = this->grid_lab[i];

        double* tdiff_L = new double[palette_size];
        double* tdiff_A = new double[palette_size];
        double* tdiff_B = new double[palette_size];

        for(int j = 0; j < palette_size; j ++) {
            LAB dir(diffpalette_L[j], diffpalette_A[j], diffpalette_B[j]);
            if(big_change(dir)) {
                LAB lab(vsrc.L + dir.L, vsrc.A + dir.A, vsrc.B + dir.B);
                RGB pc = LAB2RGB(lab);

                if(out_boundary(pc)) {
                    LAB M(oldpalette_L[j] + dir.L, oldpalette_A[j] + dir.A, oldpalette_B[j] + dir.B);
                    LAB Mdir(vsrc.L - oldpalette_L[j], vsrc.A - oldpalette_A[j], vsrc.B - oldpalette_B[j]);

                    double t1 = find_boundary(M, Mdir, 0.0, 1.0);

                    if(!out_boundary(LAB2RGB(LAB(M.L + Mdir.L, M.A + Mdir.A, M.B + Mdir.B)))) puts("M + Mdir Error");

                    double t2 = find_boundary(LAB(oldpalette_L[j], oldpalette_A[j], oldpalette_B[j]), dir, 1.0, 300.0);

                    tdiff_L[j] = dir.L - (1.0 - t1) * Mdir.L;
                    tdiff_A[j] = dir.A - (1.0 - t1) * Mdir.A;
                    tdiff_B[j] = dir.B - (1.0 - t1) * Mdir.B;

                    if(t1 > 1.0) printf("t1 > 1 at case 1 %lf, Mdir: %lf, %lf, %lf\n", t1, Mdir.L, Mdir.A, Mdir.B);
                    if(t2 < 1.0) printf("t2 < 1 at case 1 %lf\n", t2);
                }
                else {
                    double t1 = find_boundary(vsrc, dir, 1.0, 300.0);
                    double t2 = find_boundary(LAB(oldpalette_L[j], oldpalette_A[j], oldpalette_B[j]), dir, 1.0, 300.0);

                    if(t2 < 1.0) printf("t2 < 1 at case 2 %lf\n", t2);
                    double lambda = std::min(t1 / t2, 1.0);
                    tdiff_L[j] = diffpalette_L[j];
                    tdiff_A[j] = diffpalette_A[j] * lambda;
                    tdiff_B[j] = diffpalette_B[j] * lambda;
                }
            }
            else {
                tdiff_L[j] = diffpalette_L[j];
                tdiff_A[j] = diffpalette_A[j];
                tdiff_B[j] = diffpalette_B[j];
            }
        }

        LAB res = calc_single_point(param, palette_size, oldpalette_L, oldpalette_A, oldpalette_B, tdiff_L, tdiff_A, tdiff_B, vsrc);

        RGB pc = LAB2RGB(res);
        this->grid_R[i] = pc.R;
        this->grid_G[i] = pc.G;
        this->grid_B[i] = pc.B;

        this->grid_R[i] = std::max(0.0, std::min(this->grid_R[i], 255.0));
        this->grid_G[i] = std::max(0.0, std::min(this->grid_G[i], 255.0));
        this->grid_B[i] = std::max(0.0, std::min(this->grid_B[i], 255.0));
    }
}

void rbfcolor::draw(RGB* input_palette, RGB* output_palette, int palette_size) {
    double* oldpalette_L = new double[palette_size + 1];
    double* oldpalette_A = new double[palette_size + 1];
    double* oldpalette_B = new double[palette_size + 1];
    double* newpalette_L = new double[palette_size + 1];
    double* newpalette_A = new double[palette_size + 1];
    double* newpalette_B = new double[palette_size + 1];

    LAB tmplab;

    for(int i = 0; i < palette_size; i ++) {
        tmplab = RGB2LAB(RGB(input_palette[i].R * 255.0, input_palette[i].G * 255.0, input_palette[i].B * 255.0));
        oldpalette_L[i] = tmplab.L; oldpalette_A[i] = tmplab.A; oldpalette_B[i] = tmplab.B;

        tmplab = RGB2LAB(RGB(output_palette[i].R * 255.0, output_palette[i].G * 255.0, output_palette[i].B * 255.0));
        newpalette_L[i] = tmplab.L; newpalette_A[i] = tmplab.A; newpalette_B[i] = tmplab.B;
    }

    oldpalette_L[palette_size] = 0.0;
    oldpalette_A[palette_size] = 0.0;
    oldpalette_B[palette_size] = 0.0;

    double* diffpalette_L = new double[palette_size + 1];
    double* diffpalette_A = new double[palette_size + 1];
    double* diffpalette_B = new double[palette_size + 1];
    
    for(int i = 0; i < palette_size; i ++) {
        diffpalette_L[i] = newpalette_L[i] - oldpalette_L[i];
        diffpalette_A[i] = newpalette_A[i] - oldpalette_A[i];
        diffpalette_B[i] = newpalette_B[i] - oldpalette_B[i];
    }

    diffpalette_L[palette_size] = 0.0;
    diffpalette_A[palette_size] = 0.0;
    diffpalette_B[palette_size] = 0.0;

    calc_grid_result(palette_size, oldpalette_L, oldpalette_A, oldpalette_B, diffpalette_L, diffpalette_A, diffpalette_B);

    for(int i = 0; i < this->img_area; i ++) {
        double tmpR = 0.0;
        double tmpG = 0.0;
        double tmpB = 0.0;
        for(int k = 0; k < 8; k ++) {
            tmpR += this->grid_R[this->weightindex[i][k]] * this->weightmap[i][k];
            tmpG += this->grid_G[this->weightindex[i][k]] * this->weightmap[i][k];
            tmpB += this->grid_B[this->weightindex[i][k]] * this->weightmap[i][k];
        }
        this->res_r[i] = tmpR;
        this->res_g[i] = tmpG;
        this->res_b[i] = tmpB;
    }

    delete []oldpalette_L;
    delete []oldpalette_A;
    delete []oldpalette_B;
    delete []newpalette_L;
    delete []newpalette_A;
    delete []newpalette_B;
    delete []diffpalette_L;
    delete []diffpalette_A;
    delete []diffpalette_B;
}

void rbfcolor::save_img(string name) {
    vector<unsigned char> output_data;
    output_data.clear();

    for(int i = 0; i < this->img_area; i ++) {
        output_data.push_back((unsigned char)((unsigned int)floor(this->res_r[i])));
        output_data.push_back((unsigned char)((unsigned int)floor(this->res_g[i])));
        output_data.push_back((unsigned char)((unsigned int)floor(this->res_b[i])));
        output_data.push_back((unsigned char)255);
    }

    unsigned error = lodepng::encode(name, output_data, this->width, this->height);
}

void rbfcolor::gridacc_kmeans(int center_num, RGB* res) {
    int step_size = 255 / (this->ngrid - 1);
    int* sample_cnt = new int[this->grid_size];
    LAB* sample_sum = new LAB[this->grid_size]; // LAB* actually
    for(int i = 0; i < this->grid_size; i ++) {
        sample_cnt[i] = 0;
        sample_sum[i].L = 0.0;
        sample_sum[i].A = 0.0;
        sample_sum[i].B = 0.0;
    }

    for(int i = 0; i < this->img_area; i ++) {
        LAB p = RGB2LAB(RGB(this->ori_r[i], this->ori_g[i], this->ori_b[i]));
        int bin1 = round(this->ori_r[i] / (double)step_size);
        int bin2 = round(this->ori_g[i] / (double)step_size);
        int bin3 = round(this->ori_b[i] / (double)step_size);

        int bin = bin1 * ngrid * ngrid + bin2 * ngrid + bin3;

        sample_cnt[bin] ++;
        sample_sum[bin].L += p.L;
        sample_sum[bin].A += p.A;
        sample_sum[bin].B += p.B;
    }

    int tot = 0;
    for(int i = 0; i < grid_size; i ++) if(sample_cnt[i] > 0) tot ++;

    LAB* D = new LAB[tot];
    tot = 0;
    for(int i = 0; i < grid_size; i ++) if(sample_cnt[i] > 0) {
        D[tot] = LAB(sample_sum[i].L / (double)sample_cnt[i],
                     sample_sum[i].A / (double)sample_cnt[i],
                     sample_sum[i].B / (double)sample_cnt[i]
                    );
        sample_cnt[tot] = sample_cnt[i];

        tot ++;
    }

    LAB* Center = new LAB[center_num + 1];
    double* pickcnt = new double[tot];
    for(int i = 0; i < tot; i ++) pickcnt[i] = sample_cnt[i];

    for(int i = 0; i < center_num; i ++) {
        int idx = 0;
        for(int j = 0; j < tot; j ++) if(pickcnt[j] > pickcnt[idx]) idx = j;

        Center[i] = D[idx];
        for(int j = 0; j < tot; j ++) {
            double dis = 0.0;
            dis += (D[idx].L - D[j].L) * (D[idx].L - D[j].L);
            dis += (D[idx].A - D[j].A) * (D[idx].A - D[j].A);
            dis += (D[idx].B - D[j].B) * (D[idx].B - D[j].B);
            dis /= 80.0 * 80.0;

            pickcnt[j] *= 1 - exp(-dis);
        }
    }

    Center[center_num] = RGB2LAB(RGB());
    center_num ++;

    int* cnt = new int[center_num];
    LAB* sumD = new LAB[center_num];

    for(int iter = 0; iter < 20; iter ++) {
        for(int i = 0; i < center_num; i ++) {
            sumD[i] = LAB();
            cnt[i] = 0;
        }
        for(int i = 0; i < tot; i ++) {
            int min_id = -1;
            double min_v = 1e10;
            for(int j = 0; j < center_num; j ++) {
                double r = pow(D[i].L - Center[j].L, 2) + pow(D[i].A - Center[j].A, 2) + pow(D[i].B - Center[j].B, 2);
                if(r < min_v) {
                    min_v = r;
                    min_id = j;
                }
            }
            cnt[min_id] += sample_cnt[i];
            sumD[min_id].L += sample_cnt[i] * D[i].L;
            sumD[min_id].A += sample_cnt[i] * D[i].A;
            sumD[min_id].B += sample_cnt[i] * D[i].B;
        }

        for(int i = 0; i < center_num; i ++) if(cnt[i] > 0) {
            Center[i].L = sumD[i].L / (double)cnt[i];
            Center[i].A = sumD[i].A / (double)cnt[i];
            Center[i].B = sumD[i].B / (double)cnt[i];
        }
        // printf("Iteration[%d]:\n", iter);
        // for(int i = 0; i < center_num; i ++) printf("Center[%d]: (LAB)(%lf, %lf, %lf)\n", i, Center[i].L, Center[i].A, Center[i].B);
    }

    center_num --;
    for(int i = 0; i < center_num; i ++) {
        RGB ps = LAB2RGB(Center[i]);
        res[i].R = std::max(0.0, std::min(1.0, ps.R / 255.0));
        res[i].G = std::max(0.0, std::min(1.0, ps.G / 255.0));
        res[i].B = std::max(0.0, std::min(1.0, ps.B / 255.0));
    }

    delete []sample_cnt;
    delete []sample_sum;
    delete []D;
    delete []Center;
    delete []pickcnt;
    delete []cnt;
    delete []sumD;
}