#include "cp.h"
#include <cmath> 

void correlate(int ny, int nx, const float* data, float* result) {
    //double mean;
    double sum, sum_sqr;
    double ** matrix = new double*[ny];
    for (int i = 0; i < ny; ++i)
        matrix[i] = new double[nx];
    // first optimization
    for ( int i = 0; i < ny; i++ ) {
        sum = 0;
        for ( int j = 0; j < nx; j++ ) {
            sum += data[j + i * nx];
        }
        sum = sum / nx;
        sum_sqr = 0;
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j] = data[j + i * nx] - sum;
            sum_sqr += matrix[i][j] * matrix[i][j];
        }
        sum_sqr = std::sqrt(sum_sqr);
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j] /= sum_sqr;            
        }
    }
    double sab, a, b, r;
    for ( int j = 0; j < ny; j++ ) {
        for ( int i = j; i < ny; i++ ) {
            sab = 0.0;
            for (int x = 0; x < nx; x++) {
                a = matrix[i][x];
                b = matrix[j][x];
                sab += a * b;
            }
            r = nx * sab ;
            r /= nx;
            result[i + j * ny] = r;
        }
    }
    
    for (int i = 0; i < ny; ++i)
        delete [] matrix[i];
    delete [] matrix;

}
