#include <iostream>
#include <math.h>
#include <vector>
#include "dithering.h"

using namespace std;

int a, widht, height, sg;

int main(int argc, char** argv) {
    if (argc != 7) {
        cerr << "input data can't be accepted \n";
        return 1;
    }
    int grad = atoi(argv[3]);
    int dithering = atoi(argv[4]);
    int bits = atoi(argv[5]);
    double gm = atof(argv[6]);

    // открытие
    FILE * input_file; 
    input_file = fopen(argv[1], "rb");
    if (input_file == NULL) { //не удалось открыть файл
        cerr << "the file was not found \n";
        return 1;
    }
    if (fscanf(input_file, "P%i%i%i%i\n", &a, &widht, &height, &sg) != 4) {
        cerr << "invalid file \n"; //нет 4х аргументов файла
        fclose(input_file);
        return 1;
    }
    if (a != 5) {
        cerr << "unsupported type of file \n";
        fclose(input_file);
        return 1;
    }
    auto* arr_pixels = (unsigned char*)malloc(widht * height * sizeof(unsigned char));
    if (arr_pixels == NULL) {
        cerr << "allocation memory failed \n";
        fclose(input_file);
        free(arr_pixels);
        return 1;
    }

    // градиент
    if (grad == 0) {
        if (fread(arr_pixels, sizeof(unsigned char), widht * height, input_file) != widht * height) {
            cerr << "incorrect size of data in input file \n"; //некоректный размер массива данных
            fclose(input_file);
            free(arr_pixels);
            return 1;
        }
    } else if (grad == 1) {
        for (int i = 0; i < height; i++) {
            for (int j = 0; j < widht; j++) {
                double k = j / double(widht);
                arr_pixels[i * widht + j] = round(255 * k);
            }
        }
    } else {
        free(arr_pixels);
        fclose(input_file);
        cerr << "invalid grad";
        return 1;
    }

    // дизеринг
    switch (dithering) {
        case 0: {
            no_dithering(widht, height, bits, gm, &arr_pixels, grad);
            break;
        }
        case 1 : {
            dithering8(widht, height, bits, gm, &arr_pixels, grad);
            break;
        }
        case 2: {
            dith_rand(widht, height, bits, gm, &arr_pixels, grad);
            break;
        }
        case 3: {
            vector<double> err(widht * height, 0);
            floyd_steinberg(widht, height, bits, gm, &arr_pixels, err, grad);
            break;
        }
        case 4: {
            vector<double> err4(widht * height, 0);
            JJN(widht, height, bits, gm, &arr_pixels, err4, grad);
            break;
        }
        case 5: {
            vector<double> err5(widht * height, 0);
            Sierra(widht, height, bits, gm, &arr_pixels, err5, grad);
            break;
        }
        case 6: {
            vector<double> err6(widht * height, 0);
            Atkinson(widht, height, bits, gm, &arr_pixels, err6, grad);
            break;
        }
        case 7: {
            halftone(widht, height, bits, gm, &arr_pixels, grad);
            break;
        }
        default: {
            cerr << "invalid dithering";
        }
    }


    //закрытие
    FILE* output_file = fopen(argv[2], "wb");
    if (output_file == nullptr) {
        cerr << "the output file was not found \n"; //не удалось открыть файл
        fclose(input_file);
        free(arr_pixels);
        return 1;
    }
    if (fprintf(output_file, "P%i\n%i %i\n%i\n", a, widht, height, 255) == 4) {
        cerr << "invalid output data";
        fclose(input_file);
        free(arr_pixels);
        return 1;
    }
    if (fwrite(arr_pixels, sizeof(unsigned char), widht * height, output_file) != widht * height) {
        cerr << "invalid output data";
        fclose(input_file);
        free(arr_pixels);
        return 1;
    }
    fclose(input_file);
    free(arr_pixels);
    fclose(output_file);
    return 0;
}