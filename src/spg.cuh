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

#include <cuda_runtime.h>
#include <cuda.h>
#include <iostream>

// Maximum size of the p matrix, corresponding to a 64*64 matrix 
#define SIZE_P          4096


#define checkCudaErrors(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
    if (code != cudaSuccess) {
        std::cerr << "CUDA Error : " << cudaGetErrorString(code) << " " << file << " " << line << std::endl;
        if (abort) {
            exit(code);
        }
    }
}

void spg_store_matrix( unsigned int nz, float *p, float *pp );
void spg_l1( unsigned long nx, unsigned int nz, float *pt_u, float *pt_x, float *pt_w, int niter );
void spg_pos(unsigned long nx, unsigned int nz, float *pt_u, float *pt_x, int niter);
