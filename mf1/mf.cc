#include "mf.h"
#include <algorithm> 
#include <iostream>

// in[x + y*nx] and out[x + y*nx].
void mf(int ny, int nx, int hy, int hx, const float* in, float* out) {
    // limit is our boundaries of a window, beg is beging boundaries
    
    /*
    std::vector<int> v{5, 6, 4, 3, 2, 7};
    std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
    
    if ( v.size() % 2 == 1 ) {
        std::cout << "The median is " << v[v.size()/2] << '\n';
    }
    else {
        std::cout << "The median is " << ( (v[v.size()/2] + v[v.size()/2 - 1]) / 2 )  << '\n';
    }*/
    
   // std::cout << "IN" << std::endl;
   // for (int i = 0; i < ny; i++ ) {
  //      for ( int j = 0; j < nx; j++) {
  //          if (j == nx - 1)  std::cout << in[j + i * nx] << std::endl;
  //          else std::cout << in[j + i * nx] << " ";
   //     }
   // }
    std::vector<float> window_elements; 
    int limit_y, limit_x, beg_x, beg_y;
    for ( int i = 0; i < ny; i++ ) {
        limit_y = std::min(ny - 1, i + hy) + 1;
        beg_y = std::max(0, i - hy);
        
        for ( int j = 0; j < nx; j++ ) {
            
            window_elements.clear();
           
            limit_x = std::min(nx - 1, j + hx) + 1;
            beg_x = std::max(0, j - hx);
            
           // std::cout << "ITR " << i <<" "<< j << " END " << limit_y <<std::endl;
            
            for ( int k = beg_y; k < limit_y; k++ ) {
                for ( int l = beg_x; l < limit_x; l++ ) {
                    window_elements.push_back(in[l + k * nx]);
                }
            }
            
            //std::cout << "VECTOR " << i  << " "<< hy << "  " << hx << std::endl;
            
           // for (auto i = window_elements.begin(); i != window_elements.end(); ++i)
            //std::cout << *i << ' ';
          //  std::cout << std::endl;
            
            
            std::nth_element(window_elements.begin(), window_elements.begin() + window_elements.size()/2, window_elements.end());
            
            if ( window_elements.size() % 2 == 1 ) {
                 out[j + i * nx] = window_elements[window_elements.size()/2];
            }
            else {
                 auto max_half = std::max_element(window_elements.begin(), window_elements.begin() + window_elements.size()/2);
                 out[j + i * nx] = (window_elements[window_elements.size()/2] + *max_half) / 2 ;
            }
            //std::cout << in[j + i * nx];       
        }
    }
    
          // std::cout << "OUT" << std::endl;
        //    for (int i = 0; i < ny; i++ ) {
        //        for ( int j = 0; j < nx; j++) {
         //           if (j == nx - 1)  std::cout << out[j + i * nx] << std::endl;
         //           else std::cout << out[j + i * nx] << " ";
          //      }
         //   }
    
}