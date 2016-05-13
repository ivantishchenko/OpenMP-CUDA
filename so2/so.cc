#include "so.h"
#include <algorithm>
#include <omp.h>
#include <iostream>
 void quickSort ( data_t * arr, int left, int right, int nthrds) {
     
     
        if (nthrds < 2) {
            std::sort(arr + left, arr + right); 
            return;
        }
        
      int i = left, j = right;

      data_t tmp;

      data_t pivot = arr[(left + right) / 2];

    

      /* partition */

      while (i <= j) {

            while (arr[i] < pivot)
                  i++;

            while (arr[j] > pivot)
                  j--;

            if (i <= j) {

                  tmp = arr[i];

                  arr[i] = arr[j];

                  arr[j] = tmp;

                  i++;

                  j--;

            }

      }

 

      /* recursion */
      
      if (left < j) {
          #pragma omp task
            quickSort(arr, left, j, nthrds / 2);
      }
            

      if (i < right) {
          #pragma omp task
          quickSort(arr, i, right, nthrds / 2);
      }

            
 }
 void psort(int n, data_t* data) {
     int nthrds = omp_get_max_threads();
#pragma omp parallel 
#pragma omp single 
   {
       quickSort(data, 0, n, nthrds);
   }
   /*
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
    std::inplace_merge(data, data + n / 2, data + n);*/
 }