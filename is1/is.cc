#include "is.h"
#include <cstdio>
#include "vector.h"

Result segment(int ny, int nx, const float* data) {
printf("\n");
    double4_t * sums_matrix = double4_alloc((nx + 1) * (ny + 1); 
    double4_t sum = double4_0;
    // padding zzero
    for ( int i = 0; i < nx + 1; i++ ) {
        sums_matrix[i] = double4_0;
    }
    for ( int j = 0; j < ny + 1; j++ ) {
        sums_matrix[nx * j] = double4_0;
    }
    
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {

            double4_t v;
            for ( int c = 0; c < 3; ++c )
                v[c] = data[c + 3 * x + 3 * nx * y];
            v[3] = 0;

            printf("[%f,    %f,    %f] ", v[0], v[1], v[2]);
            // precalc happens here
           // sum += v;
           // sums_matrix[ (x + 1) + (y + 1) * nx] = sum;
            
            
            sums_matrix[ (x + 1) + (y + 1) * nx] = v + sums_matrix[ (x + 1) + y * nx] +  sums_matrix[ x + (y + 1) * nx] - sums_matrix[ x + y * nx];
            
        }
        printf("\n");
    }
    
    // DEBUG AND MAKE SURE SUMS_MATR is CORRECT
    
    for (int y0=0; y0<=ny-h; y0++) {
        for (int x0=0; x0<=nx-w; x0++) {
            y1 = y0 + h;
            x1 = x0 + w;
            vXc = S00[y1*snx + x1] - S00[y1*snx + x0] - S00[y0*snx + x1] + S00[y0*snx + x0];
            vYc = vPc - vXc;
            hXY4 = vXc * vXc * divX + vYc * vYc * divY;
            hXY = hXY4[0]+hXY4[1]+hXY4[2];
            if (hXY > max_hXY) {
                max_hXY = hXY;
                tx0 = x0;
                ty0 = y0;
                tx1 = x1;
                ty1 = y1;
                asm ("#dummy");
            }
        }
    }
    
    Result result { ny/3, nx/3, 2*ny/3, 2*nx/3, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f} };
    free(sums_matrix);
    return result;
}
