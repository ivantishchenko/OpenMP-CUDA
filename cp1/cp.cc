#include "cp.h"
#include <cmath> 

void correlate(int ny, int nx, const float* data, float* result) {
    //double mean;
    double sum;
    double ** matrix = new double*[ny];
        for (int i = 0; i < ny; ++i)
    matrix[i] = new double[nx];
    
    for ( int i = 0; i < ny; i++ ) {
        sum = 0;
        for ( int j = 0; j < nx; j++ ) {
            sum += data[j + i * nx];
        }
        sum = sum / nx;
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j] = data[j + i * nx] - sum;
        }
    }
    
    for ( int j = 0; j < ny; j++ ) {
        for ( int i = j; i < ny; i++ ) {
            double sab = 0.0;
            double saa = 0;
            double sbb = 0;
            for (int x = 0; x < nx; x++) {
                double a = matrix[i][x];
                double b = matrix[j][x];
                sab += a * b;
                saa += a * a;
                sbb += b * b;
            }
            double r = nx * sab ;
            r /= std::sqrt(nx * saa);
            r /= std::sqrt(nx * sbb);
            result[i + j * ny] = r;
        }
    }
    
    for (int i = 0; i < ny; ++i)
        delete [] matrix[i];
    delete [] matrix;

}
