#MakeLists.txt for CUDA
PROJECT(UK-TI) 
cmake_minimum_required(VERSION 2.8)
FIND_PACKAGE(CUDA REQUIRED)
if(CUDA_FOUND)
	message(STATUS "Building with CUDA Version:" ${CUDA_VERSION})
        INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
        list(APPEND SUMO_LIBRARIES ${CUDA_LIBRARIES} ${CUDA_DRIVER_LIBRARY})
        add_definitions(-DUSECUDA)
else()
        message(STATUS "CUDA disabled.")
endif()

option(BUILD_WITH_OPENMP "OpenMP" OFF)
if(BUILD_WITH_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DUSEOPENMP)
    endif()
endif()

set(GENCODE -gencode=arch=compute_30,code=sm_30
-gencode=arch=compute_35,code=compute_35)
set(GENCODE ${GENCODE} -gencode=arch=compute_20,code=sm_20)
set(GENCODE ${GENCODE} -gencode=arch=compute_10,code=sm_10)

# debug options
set(CUDA_NVCC_FLAGS_DEBUG "-G;-g")
set(CUDA_NVCC_FLAGS_RELEASE "-O2")

SET(CUDA_SOURCE
src/GPUKriging.cu
src/KrigingSolver.cpp
src/main.cpp
)

SET(CUDA_HEADER
src/KrigingSolver.h
src/GPUKriging.h
)

cuda_add_executable(
UKTI
${CUDA_SOURCE}
${CUDA_HEADER}
)
