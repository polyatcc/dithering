#pragma  once

#include <vector>
#include <cmath>
#include <random>
#include "dithering.h"

using  namespace std;

random_device rd;
mt19937 mersenne(rd());

double nearest_col(double pixl, int bit) {
    int k = round(pixl);
    int m = k >> (8 - bit);
    int l = 0;
    for (int i = 0; i < (7 / bit)  + 1; i++) {
        l = m + (l << bit);
    }
    l = l >> ((7 / bit + 1) * bit - 8);
    return l;
}

double correction_gamma(double this_p, double gm) {
    this_p = this_p / 255.0;
    this_p = this_p > 1 ? 1 : this_p;
    if (gm == 0) {
        if (this_p < 0.04045) {
            double temp = (255.0 * this_p) / 12.92;
            return temp;
        } else {
            double temp = 255.0 * (pow((200.0 * this_p + 11.0) / 211.0, 2.4));
            return temp;
        }
    }
    return 255.0 * pow(this_p, gm);
}

double gamma_rev(double this_p, double gm) {
    this_p = this_p / 255.0;
    this_p = this_p > 1 ? 1 : this_p;
    if (gm == 0) {
        if (this_p < 0.0031308) {
            return this_p * 12.92 * 255.0;
        } else {
            return 255.0 * ((211.0 * pow(this_p, 0.4166) - 11.0) / 200.0);
        }
    }
    double temp = 255.0 * pow(this_p, 1.0 / gm);
    return temp;
}

void no_dithering(int widht, int height, int bits, double gm, unsigned char** arr_pixels, int &grad) {
    unsigned char* tmp = *arr_pixels;
    for (int i = 0; i < widht * height; i++) {
        double this_p = tmp[i];
        this_p = correction_gamma(this_p, gm);
        this_p = nearest_col(this_p, bits);
        tmp[i] = (unsigned char)gamma_rev(this_p, gm);
    }
}

void dithering8(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double this_p = tmp[i * widht + j];
            this_p = correction_gamma(this_p, gm);
            this_p = nearest_col(this_p + 255.0 * (matrix_diz8[i % 8][j % 8]) - 0.5, bits);
            tmp[i * widht + j] = (unsigned char)(gamma_rev(this_p, gm));
        }
    }
}

void dith_rand(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double this_p = tmp[i * widht + j];
            this_p = correction_gamma(this_p, gm);
            this_p = nearest_col(this_p + 255.0 * pow(-1, j) * (double (mersenne()  % 50) / 100.), bits);
            tmp[i * widht + j] = (unsigned char) gamma_rev(this_p, gm);
        }
    }
}

double div_error(int i, int j, int widht, double gm, vector<double>& err, int bits, unsigned char* tmp) {
    double this_p = tmp[i * widht + j];
    this_p = correction_gamma(this_p, gm);
    this_p = this_p / 255.0 + err[i * widht + j] / 255.0;
    this_p = nearest_col(255.0 * this_p, bits);
    auto arr_err = correction_gamma(double(tmp[i * widht + j]), gm) + err[i * widht + j] - this_p;
    tmp[i * widht + j] = (unsigned char)gamma_rev(this_p, gm);
    return arr_err;
}

void floyd_steinberg(int widht, int height, int bits, double gm, unsigned char** arrp, vector<double>& err, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double arr_err = div_error(i, j, widht, gm, err, bits, tmp);
            if (j < widht - 1) {
                err[i * widht + j + 1] += 7.0 / 16.0 * arr_err;
            }
            if (i < height - 1) {
                err[i * widht + j + widht] += 5.0 / 16.0 * arr_err;
                if (j >= 0) {
                    err[i * widht + j + widht - 1] += 3.0 / 16.0 * arr_err;
                }
                if (j < widht - 1) {
                    err[i * widht + j + widht + 1] += 1.0 / 16.0 * arr_err;
                }
            }
        }
    }
}

void JJN(int widht, int height, int bits, double gm, unsigned char** arrp, vector<double>& err, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double arr_err = div_error(i, j, widht, gm, err, bits, tmp);
            for (int h = 0; h < 3; h++) {
                for (int w = -2; w < 3; w++) {
                    if ((i + h) < height) {
                        if ((0 <= j + w) && (j + w) < widht) {
                            err[i * widht + j + h * widht + w] += arr_err * matrix_jjn[h][w + 2];
                        }
                    }
                }
            }
        }
    }
}

void Sierra(int widht, int height, int bits, double gm, unsigned char** arrp, vector<double>& err, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double arr_err = div_error(i, j, widht, gm, err, bits, tmp);
            for (int h = 0; h < 3; h++) {
                for(int w = -2; w < 3; w++) {
                    if (i + h < height){
                        if ((0 <= j + w) && (j + w < widht)) {
                            err[i * widht + j + h * widht + w] += arr_err * matrix_sierra[h][w + 2];
                        }
                    }
                }
            }
        }
    }
}

void Atkinson(int widht, int height, int bits, double gm, unsigned char** arrp, vector<double>& err, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double arr_err = div_error(i, j, widht, gm, err, bits, tmp);
            for (int h = 0; h < 3; h++) {
                for (int w = -2; w < 3; w++) {
                    if (i + h < height) {
                        if ((0 <= j + w) && (j + w < widht)) {
                            err[w + i * widht + j + h * widht] += matrix_atk[h][w + 2] * arr_err ;
                        }
                    }
                }
            }
        }
    }
}

void halftone(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad) {
    unsigned char* tmp = *arrp;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < widht; j++) {
            double  this_p = tmp[j + i * widht];
            this_p = correction_gamma(this_p, gm);
            this_p = nearest_col(this_p + 255.0 * ((matrix_ht4[i % 4][j % 4] - 0.5)), bits);
            tmp[i * widht + j] = (unsigned char)gamma_rev(this_p, gm);
        }
    }
}

