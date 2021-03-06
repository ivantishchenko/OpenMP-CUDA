#include "so.h"
#include <algorithm>
#include <omp.h>
#include <iostream>
 void quickSort ( data_t * arr, int left, int right, int nthrds) {
     
     
        if (nthrds < 2) {
            std::sort(arr + left, arr + right+1); 
            return;
        }
        
      int i = left, j = right;

      data_t tmp;

        data_t pivot;
        int pivot_index = (left + right) / 2; 
        int min_i = std::max(pivot_index - 3, 0);
        int max_i = std::min(pivot_index + 3, right);
        
        if ( (arr[min_i] > arr[pivot_index] && arr[min_i] < arr[max_i] )|| (arr[min_i] > arr[max_i] && arr[min_i] < arr[pivot_index])) {
            pivot = arr[min_i];
        }
        else if ((arr[max_i] > arr[pivot_index] && arr[max_i] < arr[min_i] )|| (arr[max_i] > arr[min_i] && arr[max_i] < arr[pivot_index])) {
            pivot = arr[max_i];
        }
        else {
            pivot = arr[(left + right) / 2];
        }
        
     

    

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
       quickSort(data, 0, n - 1, nthrds);
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