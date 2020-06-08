#pragma  once

#include <vector>
#include <cmath>
#include <random>

using  namespace std;

random_device rd;
mt19937 mersenne(rd());

const double matrix_diz8[8][8] = {
        {0.0 / 64.0, 48.0 / 64.0, 12.0 / 64.0, 60.0 / 64.0, 3.0 / 64.0, 51.0 / 64.0, 15.0 / 64.0, 63.0 / 64.0},
        {32.0 / 64.0, 16.0 / 64.0, 44.0 / 64.0, 28.0 / 64.0, 35.0 / 64.0, 19.0 / 64.0, 47.0 / 64.0, 31.0 / 64.0},
        {8.0 / 64.0, 56.0 / 64.0, 4.0 / 64.0, 52.0 / 64.0, 11.0 / 64.0, 59.0 / 64.0, 7.0 / 64.0, 55.0 / 64.0},
        {40.0 / 64.0, 24.0 / 64.0, 36.0 / 64.0, 20.0 / 64.0, 43.0 / 64.0, 27.0 / 64.0, 39.0 / 64.0, 23.0 / 64.0},
        {2.0 / 64.0, 50.0 / 64.0, 14.0 / 64.0, 62.0 / 64.0, 1.0 / 64.0, 49.0 / 64.0, 13.0 / 64.0, 61.0 / 64.0},
        {34.0 / 64.0, 18.0 / 64.0, 46.0 / 64.0, 30.0 / 64.0, 33.0 / 64.0, 17.0 / 64.0, 45.0 / 64.0, 29.0 / 64.0},
        {10.0 / 64.0, 58.0 / 64.0, 6.0 / 64.0, 54.0 / 64.0, 9.0 / 64.0, 57.0 / 64.0, 5.0 / 64.0, 53.0 / 64.0},
        {42.0 / 64.0, 26.0 / 64.0, 38.0 / 64.0, 22.0 / 64.0, 41.0 / 64.0, 25.0 / 64.0, 37.0 / 64.0, 21.0 / 64.0},
};
const double matrix_jjn[3][5] = {
        {0.0 / 64.0, 0.0 / 64.0, 0.0 /64.0, 7.0 / 64.0, 5.0 / 64.0},
        {3.0 / 64.0, 5.0 / 64.0, 7.0 / 64.0, 5.0 / 64.0, 3.0 / 64.0},
        {1.0 / 64.0, 3.0 / 64.0, 5.0 / 64.0, 3.0 / 64.0, 1.0 / 64.0},
};
const double matrix_sierra[3][5] = {
        {0./32., 0./32., 0./32., 5./32., 3./32.},
        {2./32., 4./32., 5./32., 4./32., 2./32.},
        {0./32., 2./32., 3./32., 2./32., 0./32.},
};
const double matrix_atk[3][5] = {
        {0.0 / 8.0, 0.0 / 8.0, 0.0 / 8.0, 1.0 / 8.0, 1.0 / 8.0},
        {0.0 / 8.0, 1.0 /8.0, 1.0 / 8.0, 1.0 / 8.0, 0.0 / 8.0},
        {0.0 / 8.0, 0.0 / 8.0, 1.0 / 8.0, 0.0 / 8.0, 0.0 / 8.0},
};
const double matrix_ht4[4][4] = {
        {7./16., 13./16., 11./16., 4./16.},
        {12./16., 16./16., 14./16., 8./16.},
        {10./16., 15./16., 6./16., 2./16.},
        {5./16., 9./16., 3./16., 1./16.},
};

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

