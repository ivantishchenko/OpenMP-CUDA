#include "so.h"
#include <algorithm>
#include <omp.h>
#include <iostream>

void psort(int n, data_t* data) {
    //#pragma omp parallel num_threads(p)
    //omp_get_max_threads()
    //omp_get_thread_num()
    
    int nthrds = omp_get_max_threads();
    int parts = n / nthrds; 
    #pragma omp parallel for schedule(dynamic) 
    for ( int i = 0; i < nthrds; i++) {
        int beg = parts * i;
        int end = beg + parts;
        std::sort(data + beg, data + end);  
    }
    
    while ( nthrds > 1 ) {
        #pragma omp parallel num_threads(nthrds) 
        {
            int parts = n / nthrds; 
            int id = omp_get_thread_num();
            int beg = parts * id;
            int end = beg + parts;
            int mid = beg + parts / 2;
            std::inplace_merge(data + beg, data + mid, data + end);
        }
        nthrds = nthrds / 2;
    }
    std::inplace_merge(data, data + n / 2, data + n);
}
