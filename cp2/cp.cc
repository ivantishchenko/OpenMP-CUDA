#include "cp.h"
#include <cmath> 
#include <vector>

void correlate(int ny, int nx, const float* data, float* result) {
    //double mean;
    //double sum, sum_sqr;
    
    std::vector<std::vector<double>> matrix(ny, std::vector<double>(nx));
   
    // first optimization
    for ( int i = 0; i < ny; i++ ) {
        double sum = 0;
        for ( int j = 0; j < nx; j++ ) {
            sum += data[j + i * nx];
        }
        sum = sum / nx;
        double sum_sqr = 0;
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j] = data[j + i * nx] - sum;
            sum_sqr += matrix[i][j] * matrix[i][j];
        }
        sum_sqr = std::sqrt(sum_sqr);
        for ( int j = 0; j < nx; j++ ) {
            matrix[i][j] /= sum_sqr;            
        }
    }
  
    #pragma omp parallel for schedule(dynamic,1)
    for ( int j = 0; j < ny; j++ ) {
        for ( int i = j; i < ny; i++ ) {
            double sab = 0.0;
            for (int x = 0; x < nx; x++) {
                double a = matrix[i][x];
                double b = matrix[j][x];
                sab += a * b;
            }
            double r = nx * sab;
            r /= nx;
            result[i + j * ny] = r;
        }
    }

}
