#include "cp.h"
#include <cmath> 
#include <vector>
#include <iostream>
#include "../common/vector.h"
#define N 9

void correlate(int ny, int nx, const float* data, float* result) {

    int nx_full = nx / 8;
    int _nx = std::ceil(nx / 8.0);
    const int WIN_SIZE = 3;
    float8_t * x = float8_alloc(_nx * (ny + WIN_SIZE)); 
     //std::cout << " NX " << _nx << std::endl;    
    for (int i = 0; i < ny; i++) { 
        for (int j = 0, k = 0; j < nx_full ; j++, k += 8) { 
            float8_t tmp = { data[k + 0 + i * nx], data[k + 1 + i * nx], 
            data[k + 2 + i * nx], data[k + 3 + i * nx], data[k + 4 + i * nx], data[k + 5 + i * nx],
            data[k + 6 + i * nx], data[k + 7 + i * nx]};
            x[i * _nx  + j] = tmp; 
        }  
        if ( nx % 8 > 0) {
        float8_t tmp = float8_0;
            for ( int a = 0; a < nx % 8; a++ ) {
               int t = (nx - nx % 8 + a);
                tmp[a] = data[t + i * nx];
     
            }
             //std::cout << std::endl;
            x[i * _nx  + _nx - 1] = tmp;
        }
  
    } 
  
        // OPTIMIZATION
     for ( int i = 0; i < ny; i++ ) {
        float8_t sum = float8_0;
        for ( int j = 0; j < _nx; j++ ) {
            sum += x[i * _nx  + j];
        }
        float8_t sum_sqr = float8_0;
        float sum_one = (sum[0] + sum[1] + sum[2] + sum[3] + sum[4] + sum[5] + sum[6] + sum[7]) / nx;
        //CHANGES
        int j = 0;
        for ( j = 0; j < nx_full; j++ ) {
            
            x[i * _nx  + j] = x[i * _nx  + j] - sum_one;
            sum_sqr += x[i * _nx  + j] * x[i * _nx  + j];
            
           
        }
        
        //if ( j == _nx - 1 ) {
            for ( int a = 0; a < nx % 8; a++ ) {
               // std::cout << " T = " << nx - nx % 4 + a << std::endl; 
                x[i * _nx  + j][a] = x[i * _nx  + j][a] - sum_one;
                sum_sqr[a] += x[i * _nx  + j][a] * x[i * _nx  + j][a];
            }
       // }
        //std::cout << " T = "<< std::endl; 
            
        
        float sum_sqr_one = std::sqrt(sum_sqr[0] + sum_sqr[1] + sum_sqr[2] + sum_sqr[3] +sum_sqr[4] + sum_sqr[5] +sum_sqr[6] +sum_sqr[7]);
        for ( int j = 0; j < _nx; j++ ) {
            x[i * _nx  + j] /= sum_sqr_one;            
        }
    }


    
  // COREALTION itself
    #pragma omp parallel for schedule(dynamic)
    for ( int j = 0; j < ny; j += WIN_SIZE ) {
        for ( int i = j; i < ny; i+= WIN_SIZE ) {
            float8_t sab[N] = {float8_0, float8_0, float8_0, float8_0, float8_0, float8_0, float8_0, float8_0, float8_0}; 
       
       
            int i_rem = WIN_SIZE;
            int j_rem = WIN_SIZE;
            if ( i + WIN_SIZE > ny ) {
                i_rem = ny % WIN_SIZE;
            }
            
            if ( j + WIN_SIZE > ny ) {
                j_rem = ny % WIN_SIZE;
            }
            
            for ( int k = 0; k < _nx; k++) {
  
                sab[0] += x[i * _nx  + k] * x[j * _nx  + k];
                sab[1] += x[(i) * _nx  + k] * x[(j + 1) * _nx  + k];
                sab[2] += x[(i) * _nx  + k] * x[(j + 2) * _nx  + k];
                sab[3] += x[(i + 1) * _nx  + k] * x[(j) * _nx  + k];
                sab[4] += x[(i + 1) * _nx  + k] * x[(j + 1) * _nx  + k];
                sab[5] += x[(i + 1) * _nx  + k] * x[(j + 2) * _nx  + k];
                sab[6] += x[(i + 2) * _nx  + k] * x[(j) * _nx  + k];
                sab[7] += x[(i + 2) * _nx  + k] * x[(j + 1) * _nx  + k];
                sab[8] += x[(i + 2) * _nx  + k] * x[(j + 2) * _nx  + k];

            }


            for ( int v = 0; v < j_rem; v++ ) {
                for ( int u = 0; u < i_rem; u++ ) {
                    float tmp = 0.0;
                    for ( int p = 0; p < 8; p++ ) {
                        tmp += sab[u * WIN_SIZE + v][p];
                    }
                    
                    result[(i + u) + (j + v) * ny] = tmp;

                }
            }
                        
            
        }
    }    
       
  
    /*
    std::cout << std::endl; 
      for (int i = 0; i < ny; i++) { 
        for (int j = 0; j < nx ; j++) { 
               if ( j == nx - 1) std::cout << result[i * nx  + j] << std::endl; 
               else std::cout << result[i * nx  + j] << " ";  
        } 
    }*/
       
   // std::cout << std::endl; 
    free(x); 

}
