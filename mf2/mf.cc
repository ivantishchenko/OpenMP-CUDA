#include "mf.h"
#include <algorithm> 

void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    #pragma omp parallel for
    for ( int i = 0; i < ny; i++ ) {
        int limit_y = std::min(ny, i + hy + 1);
        int beg_y = std::max(0, i - hy);
        int row_matr = i * nx;
        for ( int j = 0; j < nx; j++ ) {
            
            std::vector<float> window_elements;
           
            int limit_x = std::min(nx, j + hx + 1);
            int beg_x = std::max(0, j - hx);
            

            for ( int k = beg_y; k < limit_y; k++ ) {
                int row_window = k * nx;
                for ( int l = beg_x; l < limit_x; l++ ) {
                    window_elements.push_back(in[l + row_window]);
                }
            }

            int win_size_half = window_elements.size()/2;

            std::nth_element(window_elements.begin(), window_elements.begin() + win_size_half, window_elements.end());
            
            if ( window_elements.size() % 2 == 1 ) {
                 out[j + row_matr] = window_elements[win_size_half];
            }
            else {
                 std::vector<float>::iterator max_half = std::max_element(window_elements.begin(), window_elements.begin() + win_size_half);
                 out[j + row_matr] = (window_elements[win_size_half] + *max_half) / 2 ;
            }
        }
    }
       
}