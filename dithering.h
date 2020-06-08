#ifndef KGIG_3_LAB_DITHERING_H
#define KGIG_3_LAB_DITHERING_H

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

double nearest_col(double pixl, int bit);
double correction_gamma(double this_p, double gm);
double gamma_rev(double this_p, double gm);
void no_dithering(int widht, int height, int bits, double gm, unsigned char** arr_pixels, int &grad);
void dithering8(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad);
void dith_rand(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad);
double div_error(int i, int j, int widht, double gm, std::vector<double>& err, int bits, unsigned char* tmp);
void floyd_steinberg(int widht, int height, int bits, double gm, unsigned char** arrp, std::vector<double>& err, int &grad);
void JJN(int widht, int height, int bits, double gm, unsigned char** arrp, std::vector<double>& err, int &grad);
void Sierra(int widht, int height, int bits, double gm, unsigned char** arrp, std::vector<double>& err, int &grad);
void Atkinson(int widht, int height, int bits, double gm, unsigned char** arrp, std::vector<double>& err, int &grad);
void halftone(int widht, int height, int bits, double gm, unsigned char** arrp, int &grad);



#endif //KGIG_3_LAB_DITHERING_H
