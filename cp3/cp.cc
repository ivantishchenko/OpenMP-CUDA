#include "cp.h"
#include <cmath> 
#include <vector>
#include <iostream>
#include "../common/vector.h"

void correlate(int ny, int nx, const float* data, float* result) {
    //double mean;
    //double sum, sum_sqr;
/*
    std::cout << std::endl;
   for (int i = 0; i < ny; i++) { 
        for (int j = 0; j < nx ; j++) { 
               if ( j == nx - 1) std::cout << data [i* nx + j] << std::endl; 
               else std::cout << data [i* nx + j] << " ";  
        } 
    }   
*/


    int nx_full = nx / 4;
    int _nx = std::ceil(nx / 4.0);
    double4_t * x = double4_alloc(_nx * ny); 
     //std::cout << " NX " << _nx << std::endl;     
    for (int i = 0; i < ny; i++) { 
        for (int j = 0, k = 0; j < nx_full ; j++, k += 4) { 
            double4_t tmp = { data[k + 0 + i * nx], data[k + 1 + i * nx], data[k + 2 + i * nx], data[k + 3 + i * nx] };
            x[i * _nx  + j] = tmp; 
            
        }  
        if ( nx % 4 > 0) {
        double4_t tmp = double4_0;
            for ( int a = 0; a < nx % 4; a++ ) {
               int t = (nx - nx % 4 + a);
                tmp[a] = data[t + i * nx];/*
                std::cout << " T " << tmp[a];
                std::cout << "(" << i << "|" << t << ")"; */
            }
             //std::cout << std::endl;
            x[i * _nx  + _nx - 1] = tmp;
        }
  
    } 
    /*
    std::cout << " AFTER PADDING" << std::endl;
        for (int i = 0; i < ny; i++) { 
        for (int j = 0; j < _nx ; j++) { 
            for ( int k = 0; k < 4; k++ ) {
               if ( k == 3) std::cout << x[i * _nx  + j][k] << std::endl; 
               else std::cout << x[i * _nx  + j][k] << " ";  
            }
        } 
    }   
    */
        // OPTIMIZATION
     for ( int i = 0; i < ny; i++ ) {
        double4_t sum = double4_0;
        for ( int j = 0; j < _nx; j++ ) {
            sum += x[i * _nx  + j];
        }
        double4_t sum_sqr = double4_0;
        double sum_one = (sum[0] + sum[1] + sum[2] + sum[3]) / nx;
        for ( int j = 0; j < _nx; j++ ) {
            if ( j == _nx - 1 ) {
                for ( int p = 0; p < 4; p++) {
                    if ( x[i * _nx  + j][p] != 0 ) x[i * _nx  + j][p] = x[i * _nx  + j][p] - sum_one;
                }
            } 
            else {
                x[i * _nx  + j] = x[i * _nx  + j] - sum_one;
            }
              
            sum_sqr += x[i * _nx  + j] * x[i * _nx  + j];
        }
        
        double sum_sqr_one = std::sqrt(sum_sqr[0] + sum_sqr[1] + sum_sqr[2] + sum_sqr[3]);
        for ( int j = 0; j < _nx; j++ ) {
            x[i * _nx  + j] /= sum_sqr_one;            
        }
    }
    /*
     std::cout << " AFTER OPTIMIZATION" << std::endl;
           // OUTPUT
    for (int i = 0; i < ny; i++) { 
        for (int j = 0; j < _nx ; j++) { 
            for ( int k = 0; k < 4; k++ ) {
               if ( k == 3) std::cout << x[i * _nx  + j][k] << std::endl; 
               else std::cout << x[i * _nx  + j][k] << " ";  
            }
        } 
    }   */

  // COREALTION itself
    #pragma omp parallel for schedule(dynamic)
    for ( int j = 0; j < ny; j++ ) {
        for ( int i = j; i < ny; i++ ) {
            double4_t sab = double4_0;
            for (int l = 0; l < _nx; l++) {
                double4_t a = x[i * _nx  + l];
                double4_t b = x[j * _nx  + l];
                sab += a * b;
            }
            result[i + j * ny] = sab[0] + sab[1] + sab[2] + sab[3];
        }
    }
   // std::cout << std::endl; 
    free(x); 

}
