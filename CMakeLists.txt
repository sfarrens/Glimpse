cmake_minimum_required(VERSION 2.8)

# Paths to external libraries
set(NICAEA_URL http://www.cosmostat.org/wp-content/uploads/2014/12/nicaea_2.5.tgz)
set(NICAEA_MD5 27c4d6d46f234133fc69130e6460c2d4)

set(SPARSE2D_URL https://github.com/CosmoStat/Sparse2D.git)

set(NFFT_URL https://www-user.tu-chemnitz.de/~potts/nfft/download/nfft-3.3.0.tar.gz)
set(NFFT_MD5 8efe6a43a9565f61ebf836a314b9aeb2)

# Adding customized cmake module
list(APPEND CMAKE_MODULE_PATH  "${CMAKE_SOURCE_DIR}/cmake/Modules/")

include(FindPkgConfig)
include(ExternalProject)

project(glimpse)

# Include external libraries
find_package(GSL    REQUIRED)
find_package(CCFits REQUIRED)
find_package(FFTW   REQUIRED)
find_package(NFFT   REQUIRED)

# Build boost manually on Mac OS X systems otherwise just use the system version
if(APPLE)
	include(BuildBoost)
else(APPLE)
	find_package(Boost  COMPONENTS program_options REQUIRED)
endif(APPLE)

# Include external projects
## NICAEA
ExternalProject_Add(NICAEA
        PREFIX NICAEA
        URL ${NICAEA_URL}
        URL_MD5 ${NICAEA_MD5}
        PATCH_COMMAND patch -p1 < ${PROJECT_SOURCE_DIR}/nicaea.patch)
        set(NICAEA_LIBRARY_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/lib)
        set(NICAEA_INCLUDE_DIR ${CMAKE_BINARY_DIR}/NICAEA/src/NICAEA/include)
        set(NICAEA_LIBRARIES -lnicaea)

## Sparse2D
ExternalProject_Add(SPARSE2D
        PREFIX SPARSE2D
        GIT_REPOSITORY ${SPARSE2D_URL}
        CMAKE_ARGS     -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/extern)
        set(SPARSE2D_LIBRARY_DIR ${CMAKE_BINARY_DIR}/extern/lib)
        set(SPARSE2D_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/include)
        set(SPARSE2D_LIBRARIES -lsparse2d)
        
include_directories(${SPARSE2D_INCLUDE_DIR}
		    ${GSL_INCLUDE_DIR}
		    ${CCFITS_INCLUDE_DIR}
		    ${CFITSIO_INCLUDE_DIR}
		    ${Boost_INCLUDE_DIR}
		    ${NICAEA_INCLUDE_DIR}
		    ${NFFT_INCLUDE_DIR}
		    ${FFTW_INCLUDE_DIR})
		    
link_directories(${SPARSE2D_LIBRARY_DIR}
		    ${GSL_LIBRARY_DIR}
		    ${CCFITS_LIBRARY_DIR}
		    ${CFITSIO_LIBRARY_DIR}
		    ${Boost_LIBRARY_DIR}
		    ${NICAEA_LIBRARY_DIR}
		    ${NFFT_LIBRARY_DIR}
		    ${FFTW_LIBRARY_DIR})

set(GLIMPSE_SRC src/survey.cpp
		src/redshift_distribution.cpp
		src/field.cpp
		src/surface_reconstruction.cpp
		src/starlet_2d.cpp
		src/wavelet_transform.cpp)

if(APPLE)
  set(CMAKE_CXX_FLAGS "-DMACOS -O3 -fpermissive")
  set(CMAKE_C_FLAGS "-O3 -fpermissive")
else(APPLE)
  set(CMAKE_CXX_FLAGS "-O3 -fpermissive -fopenmp")
  set(CMAKE_C_FLAGS "-O3 -fpermissive")
endif(APPLE)
  
add_executable(glimpse src/glimpse.cpp ${GLIMPSE_SRC})
add_dependencies(glimpse NICAEA SPARSE2D)
target_link_libraries(glimpse ${NICAEA_LIBRARIES}
			      ${SPARSE2D_LIBRARIES}
			      ${Boost_LIBRARIES} 
			      ${GSL_LIBRARIES}
			      ${NFFT_LIBRARIES}
			      ${FFTW_LIBRARIES}
			      ${CCFITS_LIBRARY}
			      ${CFITSIO_LIBRARY})
