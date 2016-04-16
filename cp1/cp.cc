#include "cp.h"
#include <cmath> 

void correlate(int ny, int nx, const float* data, float* result) {
    double mean;
    double** matrix = new double*[ny];
        for (int i = 0; i < ny; ++i)
    matrix[i] = new double[nx];
    
    for ( int i = 0; i < ny; i++ ) {
        mean = 0;
        for ( int j = 0; j < nx; j++ ) {
            mean += data[j + i * nx];
        }
        mean /= nx;
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j]= data[j + i * nx] - mean;
        }
    }
    
    //for ( int i = 0; i < ny; i++ ) {
    //    for ( int j = 0; j < nx; j++ ) {
    //        matrix[i][j] = matrix[i][j] / sum;
    //    }
    //}
    
    for ( int j = 0; j < ny; j++ ) {
        for ( int i = j; i < ny; i++ ) {
            double sa = 0.0;
            double sb = 0.0;
            double sab = 0.0;
            double saa = 0.0;
            double sbb = 0.0;
            for (int x = 0; x < nx; x++) {
                double a = matrix[i][x];
                double b = matrix[j][x];

                sab += a * b;
                saa += a * a;
                sbb += b * b;
            }
            double r = nx * sab - sa * sb;
            r /= std::sqrt(nx * saa - sa * sa);
            r /= std::sqrt(nx * sbb - sb * sb);
            result[i + j * ny] = r;
        }
    }
    
    for (int i = 0; i < ny; ++i)
        delete [] matrix[i];
    delete [] matrix;

}
