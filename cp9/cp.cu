#include "cp.h"
#include <cuda_runtime.h>
#include <iostream>

template<typename T>
cudaError_t allocDev(T*& d_p, size_t elements) {
    return cudaMalloc((void**)&d_p, elements * sizeof(T));
}

template<typename T>
cudaError_t allocHost(T*& d_p, size_t elements) {
    return cudaHostAlloc((void**)&d_p, elements * sizeof(T), cudaHostAllocMapped);
}

__global__ void medianKernel(float * output, const float * input, int ny, int nx) {
    int x = threadIdx.x + blockIdx.x * blockDim.x; 
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    if (x >= ny || y >= ny) return;
    
    //int index = x + nx * y; 
    //int res = index + 1;   
    //output[index] = res; 
    /// CORELATE
 
 
    if (y > x) {
        return;
    }
    float sab, a, b;

            sab = 0.0;
            for (int j = 0; j < nx; j++) {
                a = input[x * nx + j];
                b = input[y * nx + j];
                sab += a * b;
                //printf("a = %d ; b = %d\n", x, y);
            }
            
            output[x + y * ny] = sab;
        
            //printf("%d  ", sab);
    //CORELaTE
    
    
}

void correlate(int ny, int nx, const float* data, float* result) {
    const size_t N = nx * ny;      
    float * input_CPU = nullptr;
   // double * output_CPU = nullptr;  
    float * input_GPU = nullptr;  
    float * output_GPU = nullptr;  
    
    // ALOC CPU
    cudaError_t res_in_host = allocHost(input_CPU, N);      // modifies d, input is not bytes
    //cudaError_t res_out_host = allocHost(output_CPU, N);      // modifies d, input is not bytes
    // ALOC GPU
    cudaError_t res_in_dev = allocDev(input_GPU, N);     
    cudaError_t res_out_dev = allocDev(output_GPU, ny * ny); 
    //Check for errors
  //  if (res_in_host == cudaSuccess && res_out_host == cudaSuccess ) std::cout << "Allocated host memory" << std::endl;
   // if (res_in_dev == cudaSuccess && res_out_dev == cudaSuccess ) std::cout << "Allocated host memory" << std::endl;

    //OPTIMAZATION
    
    float sum, sum_sqr;

    // first optimization
    for ( int i = 0; i < ny; i++ ) {
        sum = 0;
        for ( int j = 0; j < nx; j++ ) {
            sum += data[j + i * nx];
        }
        sum = sum / nx;
        sum_sqr = 0;
        for ( int j = 0; j < nx; j++ ) {
            input_CPU[j + i * nx] = data[j + i * nx] - sum;
            sum_sqr += input_CPU[j + i * nx] * input_CPU[j + i * nx];
        }
        sum_sqr = std::sqrt(sum_sqr);
        for ( int j = 0; j < nx; j++ ) {
            input_CPU[j + i * nx] /= sum_sqr;            
        }
    }
    
  
        // CPU -> GPU
    cudaMemcpy(input_GPU, input_CPU, N * sizeof(float), cudaMemcpyHostToDevice);
    
    // block size
    dim3 dimBlock(8, 8);
    dim3 dimGrid((ny + dimBlock.x - 1) / dimBlock.x, (ny + dimBlock.y - 1) / dimBlock.y);
    //launch KERNEL
    medianKernel<<<dimGrid, dimBlock>>>(output_GPU, input_GPU, ny, nx);
    // GPU -> CPU
    cudaMemcpy(result, output_GPU, ny * ny * sizeof(float), cudaMemcpyDeviceToHost);
    
    
    
    ///std::cout << std::endl << "HELLO WORLD    ["<<output_CPU[3] << "] :from the GPU" << std::endl;
    
    cudaFree(input_CPU);
    //cudaFree(output_CPU);
    cudaFree(input_GPU);
    cudaFree(output_GPU);
}