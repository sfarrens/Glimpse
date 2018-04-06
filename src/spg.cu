/*! Copyright CEA, 2015-2016
 * author : Francois Lanusse < francois.lanusse@gmail.com >
 *
 * This software is a computer program whose purpose is to reconstruct mass maps
 * from weak gravitational lensing.
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 *
 */
#include "spg.cuh"

#include <cuda.h>

// Constant memory array used to store the pre-conditionning matrix
__constant__ float PP[SIZE_P];
__constant__ float P[SIZE_P];


void spg_store_matrix ( unsigned int nz, float * p, float *pp )
{

    cudaMemcpyToSymbol ( P,   p, nz*nz*sizeof ( float ) );
    cudaMemcpyToSymbol ( PP, pp, nz*nz*sizeof ( float ) );
}

__global__ void spg_l1_kernel (unsigned long nx, float *pt_u, float *pt_x, float * pt_w, int niter )
{

    extern __shared__ float cache[];
    
    unsigned long x = threadIdx.x;
    unsigned long z = threadIdx.y;
    unsigned long x_offset = blockIdx.x * blockDim.x;
    unsigned long nz = blockDim.y;
    unsigned long idx = z * blockDim.x + x;
    float bb=0.0005;
    float g = 0;
    float gg0 = 0;
    float epsilon = 1e-3;
    float epsilon_0 = 1e-5;

    // Load part of the matrix into shared memory
    float u = 0;
    float w = 0;
    float gold = 0;
    float uold = 0;

    // Initialize the loop by computing the transform of input vector
    double px = 0;
    if ( x + x_offset  < nx ) {
        cache[idx] = pt_x[z * nx + x + x_offset];
        u          = pt_u[z * nx + x + x_offset];
        w          = pt_w[z * nx + x + x_offset];
    }else{
        cache[idx] = 0;
    }

    __syncthreads();

    // Compute matrix product
    for ( unsigned int z1=0; z1 < nz; z1++ ) {
        px += cache[z1 * blockDim.x + x] * P[z1 * nz + z];
    }
    
    __syncthreads();
    
    // Compute l2 norm of initial gradient
    cache[idx] = px * px;
    __syncthreads();
    for ( unsigned int s=nz/2; s>0; s>>=1 ) {
        if ( z < s ) {
            cache[idx] += cache[ ( z+s ) * blockDim.x + x];
        }

        __syncthreads();
    }
    gg0 = sqrtf(cache[x]);

    // Main loop of minimisation algorithm
    for ( int iter=0; iter < niter; iter++ ) {
        // Initialize gradient with A^tx
        g =-px;

        // Fill the cache for performing matrix multiplication
        __syncthreads();
        cache[z * blockDim.x + x] = u;

        __syncthreads();

        // Compute matrix product, which gives g[z]
        for ( unsigned int z1=0; z1 < nz; z1++ ) {
            g += PP[z1 * nz + z] * cache[z1 * blockDim.x + x];
        }

        __syncthreads();
        // Compute BB steps //////////////////////////////////
        if ( iter>0 ) {

            cache[idx] = ( g - gold ) * ( u - uold );
            __syncthreads();
            for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                if ( z < s ) {
                    cache[idx] += cache[ ( z + s ) * blockDim.x + x];
                }

                __syncthreads();
            }

            float sy = cache[x];
            __syncthreads();

            if ( iter % 2 ) {
                // Reduction on yy accross threads
                cache[idx] = ( g - gold ) * ( g - gold );
                __syncthreads();
                for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                    if ( z < s ) {
                        cache[idx] += cache[ ( z+s ) * blockDim.x + x];
                    }

                    __syncthreads();
                }
                bb = cache[x] == 0 ? 0 : sy/cache[x];
            } else {
                // Reduction on ss accross threads
                cache[idx] = ( u - uold ) * ( u - uold );
                __syncthreads();
                for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                    if ( z < s ) {
                        cache[idx] += cache[ ( z+s ) * blockDim.x + x];
                    }

                    __syncthreads();
                }
                bb = sy == 0 ? 0 : cache[x]/sy;
            }
        }
        ////////////////////////////////////////////////
        // Compute optimality check
        __syncthreads();
        float o;
        if(x + x_offset  >= nx || sqrtf(gg0) < epsilon_0){
            cache[idx] = 0;
        }else{
            if ( u <=  -w + 1e-5*w ) {
                o = fminf ( 0.0f, g );
            } else if ( u >= w - 1e-5*w ) {
                    o = fmaxf ( 0.0f, g );
            } else {
                    o = g;
            }
            cache[idx] = o*o/gg0; 
        }
        
        // Reduction across redshifts
        __syncthreads();
        for ( unsigned int s=nz/2; s>0; s>>=1 ) {
            if ( z < s ) {
                cache[idx] += cache[ ( z+s ) * blockDim.x + x];
            }

            __syncthreads();
        }
        
        if(x + x_offset  >= nx){
            cache[idx] = 0;
        }else if(z == 0)
            cache[idx] = sqrtf(cache[idx]);

            __syncthreads();
        // Update variable only if optimality is not reached
        if(cache[x] > epsilon) {
            uold = u;
            gold = g;
            
            // Small check to prevent all hell from breaking loose
            if ( fabsf ( bb ) > 1.0 ) {
                bb=0.001;
            }
            u = u - bb*g;

            // Apply projection
            u = u - copysignf (fdimf(fabsf ( u ), w ), u);
        }
        
        // Reduction across lines of sights
        for ( unsigned int s=blockDim.x/2; s>0; s>>=1 ) {
            if ( z == 0 && x < s ) {
                cache[x] = fmaxf(cache[x], cache[x + s]);
            }
            __syncthreads();
        }

        // Decide to exit loop if optimality is reached for every l.o.s.
        if ( cache[0] <= epsilon ) {
            break;
        }
    }

    // Compute transform of resulting array
    __syncthreads();
    cache[idx] = u - copysignf (fdimf(fabsf ( u ), w ), u);

    __syncthreads();
    px = 0;
    // Compute matrix product
    for ( unsigned int z1=0; z1 < nz; z1++ ) {
        px += cache[z1 * blockDim.x + x] * P[z1 * nz + z];
    }

    if ( x + x_offset  < nx ) {
        // Save the variables back to the main memory
        pt_u[z * nx + x + x_offset] = u;
        pt_x[z * nx + x + x_offset] = px;
    }
}


void spg_l1 ( unsigned long nx, unsigned int nz,  float *pt_u, float *pt_x, float *pt_w, int niter)
{
    dim3 dimBlock ( 1024 / nz, nz );
    dim3 dimGrid ( nx / dimBlock.x + ( nx % dimBlock.x > 0 ? 1 : 0 ) );

    spg_l1_kernel<<< dimGrid, dimBlock, dimBlock.x * dimBlock.y * sizeof ( float ) >>> (nx, pt_u, pt_x, pt_w, niter);
}

__global__ void spg_pos_kernel (unsigned long nx, float *pt_u, float *pt_x, int niter )
{

    extern __shared__ float cache[];

    unsigned int x = threadIdx.x;
    unsigned int z = threadIdx.y;
    unsigned int x_offset = blockIdx.x * blockDim.x;
    unsigned int nz = blockDim.y;
    unsigned int idx = z * blockDim.x + x;
    float bb=0.0005;
    float g = 0;
    float gg0 = 0;
    float epsilon = 1e-4;
    float epsilon_0 = 1e-5;
    float u = 0;
    float gold = 0;
    float uold = 0;

    // Initialize the loop by computing the transform of input
    // Load part of the matrix into shared memory vector
    double px = 0;
    if ( x + x_offset  < nx ) {
        cache[idx] = pt_x[z * nx + x + x_offset];
        u = pt_u[z * nx + x + x_offset];
    }else{
        cache[idx] = 0;
    }

    __syncthreads();

    // Compute matrix product
    for ( unsigned int z1=0; z1 < nz; z1++ ) {
        px += cache[z1 * blockDim.x + x] * P[z1 * nz + z];
    }
    
    __syncthreads();
    
    // Compute l2 norm of initial gradient
    cache[idx] = px * px;
    __syncthreads();
    for ( unsigned int s=nz/2; s>0; s>>=1 ) {
        if ( z < s ) {
            cache[idx] += cache[ ( z+s ) * blockDim.x + x];
        }

        __syncthreads();
    }
    gg0 = sqrtf(cache[x]);

    // Main loop of minimisation algorithm
    for ( int iter=0; iter < niter; iter++ ) {
        // Initialize gradient with A^tx
        g =-px;

        // Fill the cache for performing matrix multiplication
        __syncthreads();
        cache[z * blockDim.x + x] = u;

        __syncthreads();

        // Compute matrix product, which gives g[z]
        for ( unsigned int z1=0; z1 < nz; z1++ ) {
            g += PP[z1 * nz + z] * cache[z1 * blockDim.x + x];
        }

        __syncthreads();
        // Compute BB steps //////////////////////////////////
        if ( iter>0 ) {

            cache[idx] = ( g - gold ) * ( u - uold );
            __syncthreads();
            for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                if ( z < s ) {
                    cache[idx] += cache[ ( z + s ) * blockDim.x + x];
                }

                __syncthreads();
            }

            float sy = cache[x];
            __syncthreads();

            if ( iter % 2 ) {
                // Reduction on yy accross threads
                cache[idx] = ( g - gold ) * ( g - gold );
                __syncthreads();
                for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                    if ( z < s ) {
                        cache[idx] += cache[ ( z+s ) * blockDim.x + x];
                    }

                    __syncthreads();
                }
                bb = cache[x] == 0 ? 0 : sy/cache[x];
            } else {
                // Reduction on ss accross threads
                cache[idx] = ( u - uold ) * ( u - uold );
                __syncthreads();
                for ( unsigned int s=nz/2; s>0; s>>=1 ) {
                    if ( z < s ) {
                        cache[idx] += cache[ ( z+s ) * blockDim.x + x];
                    }

                    __syncthreads();
                }
                bb = sy == 0 ? 0 : cache[x]/sy;
            }
        }
        ////////////////////////////////////////////////
        // Compute optimality check
        __syncthreads();
        if(x + x_offset  >= nx || sqrtf(gg0) < epsilon_0){
            cache[idx] = 0;
        } else if ( u ==  0) {
            cache[idx] = fmaxf(g, 0) * fmaxf(g, 0)/gg0;
        }else {
            cache[idx] = g*g/gg0;
        }

        // Reduction across redshifts
        __syncthreads();
        for ( unsigned int s=nz/2; s>0; s>>=1 ) {
            if ( z < s ) {
                cache[idx] += cache[ ( z+s ) * blockDim.x + x];
            }

            __syncthreads();
        }
        
        if(x + x_offset  >= nx){
            cache[idx] = 0;
        }else if(z == 0)
            cache[idx] = sqrtf(cache[idx]);

            __syncthreads();
        // Update variable only if optimality is not reached
        if(cache[x] > epsilon) {
            uold = u;
            gold = g;
            
            // Small check to prevent all hell from breaking loose
            if ( fabsf ( bb ) > 1.0 ) {
                bb=0.001;
            }
            u = u - bb*g;

            // Apply projection
            u = u - fmaxf ( u, 0.f );
        }
        
        // Reduction across lines of sights
        for ( unsigned int s=blockDim.x/2; s>0; s>>=1 ) {
            if ( z == 0 && x < s ) {
                cache[x] = fmaxf(cache[x], cache[x + s]);
            }
            __syncthreads();
        }

        // Decide to exit loop if optimality is reached for every l.o.s.
        if ( cache[0] <= epsilon ) {
            break;
        }
    }

    // Compute transform of resulting array
    __syncthreads();
    cache[idx] =  u - fmaxf ( u, 0.f );

    __syncthreads();
    px = 0;
    // Compute matrix product
    for ( unsigned int z1=0; z1 < nz; z1++ ) {
        px += cache[z1 * blockDim.x + x] * P[z1 * nz + z];
    }

    if ( x + x_offset  < nx ) {
        // Save the variables back to the main memory
        pt_u[z * nx + x + x_offset] = u;
        pt_x[z * nx + x + x_offset] = px;
    }
}

void spg_pos ( unsigned long nx, unsigned int nz,  float *pt_u, float *pt_x, int niter )
{
    dim3 dimBlock ( 1024/nz, nz );
    dim3 dimGrid ( nx/dimBlock.x + ( nx % dimBlock.x > 0 ? 1 : 0 ) );

    spg_pos_kernel<<< dimGrid, dimBlock, dimBlock.x * dimBlock.y * sizeof ( float ) >>> (nx, pt_u, pt_x, niter );
}
