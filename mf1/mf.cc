#include "mf.h"
#include <algorithm> 
#include <iostream>


void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    std::vector<float> window_elements; 
    int limit_y, limit_x, beg_x, beg_y;
    std::vector<float>::iterator max_half;
   
    for ( int i = 0; i < ny; i++ ) {
        limit_y = std::min(ny, i + hy + 1);
        beg_y = std::max(0, i - hy);
        
        for ( int j = 0; j < nx; j++ ) {
            
            window_elements.clear();
           
            limit_x = std::min(nx, j + hx + 1);
            beg_x = std::max(0, j - hx);
            

            for ( int k = beg_y; k < limit_y; k++ ) {
                for ( int l = beg_x; l < limit_x; l++ ) {
                    window_elements.push_back(in[l + k * nx]);
                }
            }

            std::nth_element(window_elements.begin(), window_elements.begin() + window_elements.size()/2, window_elements.end());
            
            if ( window_elements.size() % 2 == 1 ) {
                 out[j + i * nx] = window_elements[window_elements.size()/2];
            }
            else {
                 max_half = std::max_element(window_elements.begin(), window_elements.begin() + window_elements.size()/2);
                 out[j + i * nx] = (window_elements[window_elements.size()/2] + *max_half) / 2 ;
            }
        }
    }   
}