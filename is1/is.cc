#include "is.h"
#include <cstdio>
#include "vector.h"

Result segment(int ny, int nx, const float* data) {
//printf("\n");
    double4_t * sums_matrix = double4_alloc((nx + 1) * (ny + 1)); 
    // SUMS_MATRIX
    for ( int i = 0; i < nx + 1; i++ ) {
        sums_matrix[i] = double4_0;
    }
    for ( int j = 0; j < ny + 1; j++ ) {
        sums_matrix[(nx + 1) * j] = double4_0;
    }
    /*
      for (int y = 0; y < ny; ++y) {
            double4_t sum = double4_0;
            for (int x = 0; x < nx; ++x) {
                for ( int c = 0; c < 3; ++c )
                    sum[c] += data[c + 3 * x + 3 * nx * y];
                    
                sums_matrix[ (x + 1) + (y + 1) * (nx + 1)] = sum + sums_matrix[ (x + 1) + y  * (nx + 1)];
            }
      }
    */
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {

            double4_t v;
            for ( int c = 0; c < 3; ++c )
                v[c] = data[c + 3 * x + 3 * nx * y];
            v[3] = 0;

            //printf("[%f,    %f,    %f] ", v[0], v[1], v[2]);
            // precalc happens here
           // sum += v;
           // sums_matrix[ (x + 1) + (y + 1) * nx] = sum;
            
            
            sums_matrix[ (x + 1) + (y + 1) * (nx + 1)] = v + sums_matrix[ (x + 1) + y * (nx + 1)] +  sums_matrix[ x + (y + 1) * (nx + 1)] - sums_matrix[ x + y * (nx + 1)];
            
        }
        //printf("\n");
    }
   //printf("ALL SUMS\n");
    
    //printf("[%f,    %f,    %f] ", vPc[0], vPc[1], vPc[2]);
    
    //printf("\n");
    /*
    for ( int i = 0; i < ny + 1; i++ ) {
        for ( int j = 0; j < nx + 1; j++ ) {


                
   
                

                printf("[%f,    %f,    %f] ", sums_matrix[j + i * (nx + 1)][0], sums_matrix[j + i * (nx + 1) ][1], sums_matrix[j + i * (nx + 1)][2]);
                if (j == nx) printf("\n");
        }

    }*/
    double4_t vPc = sums_matrix[nx + ny * (nx + 1)];
    int snx = nx + 1;
    int tx0 = 0;
    int ty0 = 0; 
    int tx1 = 0; 
    int ty1 = 0;
    double max_hXY = 0;
    double4_t vPc_sqr = vPc * vPc;
    // DEBUG AND MAKE SURE SUMS_MATR is CORRECT
    #pragma omp parallel for schedule(dynamic)
    for ( int h = 1; h <= ny; h++ ) {
          int lx0 = 0;
          int ly0 = 0; 
          int lx1 = 0; 
          int ly1 = 0;
          double lmax_hXY = 0;
        for ( int w = 1; w <= nx; w++ ) {

           if ( ny * nx != w * h)  {
            double divX = 1.0 / (w * h);
            double divY = 1.0 / (ny * nx  - w * h);
                
                double divXdivY = divX + divY;
                double4_t vPc2divY = 2 * vPc * divY;
                double4_t vPcdivY = vPc_sqr * divY;
                
                for ( int y0 = 0; y0 <= ny - h; y0++ ) {
                    for ( int x0 = 0; x0 <= nx - w; x0++ ) {
                        int y1 = y0 + h;
                        int x1 = x0 + w;

                        // dirty fix here
                        double4_t vXc = sums_matrix[y1 * snx + x1] - sums_matrix[y1 * snx + x0] - sums_matrix[y0 * snx + x1] + sums_matrix[y0 * snx + x0];
                        //double4_t vYc = vPc - vXc;
                        double4_t hXY4 = vXc*(vXc*(divXdivY)-(vPc2divY)) + (vPcdivY);

                        //double4_t hXY4 = vXc * vXc * divX + vYc * vYc * divY;
                        double hXY = hXY4[0] + hXY4[1] + hXY4[2];
                        if ( hXY > lmax_hXY ) {
                            //double4_t tmp_out = vYc * divY;
                            //double4_t tmp_inn = vXc * divX;
                                lmax_hXY = hXY;
                               // outer = tmp_out;
                               // inner = tmp_inn;
                                lx0 = x0;
                                ly0 = y0;
                                lx1 = x1;
                                ly1 = y1;
                                asm ("#nopnopnop");
                 
                        }
                   
                        
                    }
                }
            }

        }
        if ( lmax_hXY > max_hXY ) {
                #pragma omp critical 
                            {
                                if ( lmax_hXY > max_hXY ) {
                                max_hXY = lmax_hXY;
                               // outer = tmp_out;
                               // inner = tmp_inn;
                                tx0 = lx0;
                                ty0 = ly0;
                                tx1 = lx1;
                                ty1 = ly1;
                                asm ("#nopnopnop");
                                
                                }
                            }
        }
        
    }
  int w = tx1 - tx0;
  int h = ty1 - ty0;
  double divX = 1.0 / (w * h);
  double divY = 1.0 / (ny * nx  - w * h);
  double4_t vXc = sums_matrix[ty1 * snx + tx1] - sums_matrix[ty1 * snx + tx0] - sums_matrix[ty0 * snx + tx1] + sums_matrix[ty0 * snx + tx0];
  double4_t vYc = vPc - vXc;
  double4_t inner = vXc * divX; 
  double4_t outer = vYc * divY;

    //Result result { ty0, tx0, ty1, tx1, {0.322511, 0.404981, 0.809137}, {(float)inner[0], (float)inner[1], (float)inner[2]} };
    Result result { ty0, tx0, ty1, tx1, {(float)outer[0], (float)outer[1], (float)outer[2]}, {(float)inner[0], (float)inner[1], (float)inner[2]} };
    free(sums_matrix);
    return result;
}
