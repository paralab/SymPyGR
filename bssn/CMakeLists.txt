cmake_minimum_required(VERSION 2.8)
project(bssn)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

option(WITH_CUDA " build with cuda " ON)
option(BSSN_ENABLE_CUDA "enable cuda execution" ON)


if(BSSN_ENABLE_CUDA)
    add_definitions(-DBSSN_ENABLE_CUDA)
endif()


set(INC_FILES
        include/derivs.h
        include/rhs.h
        include/dendro.h
        include/bssn_profiler.h
        include/profile_param.h
        include/block.h
        include/grUtils.h
        )

set(SRC_FILES
        src/derivs.cpp
        src/rhs.cpp
        src/bssn_profiler.cpp
        src/profile_param.cpp
        src/block.cpp
        src/grUtils.cpp
        )

if(WITH_CUDA)
    set(INC_FILES ${INC_FILES} include/cudaUtils.h)
    set(SRC_FILES ${SRC_FILES} src/cudaUtils.cpp)
endif()

add_library(bssnUtils ${INC_FILE} ${SRC_FILES})
target_include_directories(bssnUtils PRIVATE include)
target_include_directories(bssnUtils PRIVATE utils)


if(WITH_CUDA)
    find_package(CUDA REQUIRED)
    if(CUDA_FOUND)
        set(CUDA_PROPAGATE_HOST_FLAGS OFF)
        set(CUDA_SEPARABLE_COMPILATION ON)

        #list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
        #list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
        list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
        #list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)

        #list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_60,code=compute_60)
        #list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_61,code=compute_61)
        list( APPEND CUDA_NVCC_FLAGS -use_fast_math )
        list( APPEND CUDA_NVCC_FLAGS -std=c++11)
        list( APPEND CUDA_NVCC_FLAGS -maxrregcount 64)
    endif()
endif()

if(WITH_CUDA)

    set(CUDA_GR_INCLUDE_H
            cuda_gr/include/bssn_rhs_deriv_mem_cuda.h
            cuda_gr/include/params_cu.h
            cuda_gr/include/profile_gpu.h
        )



    set(CUDA_GR_SRC_CPP
            cuda_gr/src/bssn_rhs_deriv_mem_cuda.cpp
            cuda_gr/src/params_cu.cpp
            cuda_gr/src/profile_gpu.cpp
            )

endif()

set(CUDA_GR_INC
        cuda_gr/include/cudaUtils.cuh
        cuda_gr/include/derivs.cuh
        #cuda_gr/include/rhs.cuh
        cuda_gr/include/derivs_bssn.cuh
        cuda_gr/include/rhs_bssn.cuh
        cuda_gr/include/rhs_cuda.cuh
        )

set(CUDA_GR_SRC
        cuda_gr/src/cudaUtils.cu
        cuda_gr/src/derivs.cu
        #cuda_gr/src/rhs.cu
        cuda_gr/src/derivs_bssn.cu
        cuda_gr/src/rhs_bssn.cu
        cuda_gr/src/rhs_cuda.cu
        )

if(WITH_CUDA)
    cuda_add_library(bssn_cuda ${CUDA_GR_INC} ${CUDA_GR_SRC})
    target_include_directories(bssn_cuda PRIVATE cuda_gr/include)
    target_include_directories(bssn_cuda PRIVATE include)
endif()


set(SOURCE_FILES cuda_gr/include/gpuBSSNExample.h cuda_gr/src/gpuBSSNExample.cpp ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP} ${INC_FILES} ${SRC_FILES} )
add_executable(cudaBSSN ${SOURCE_FILES})
target_include_directories(cudaBSSN PRIVATE include)
target_include_directories(cudaBSSN PRIVATE utils)
target_include_directories(cudaBSSN PRIVATE cuda_gr/include)
target_link_libraries(cudaBSSN bssnUtils m)

if(WITH_CUDA)
    target_include_directories(cudaBSSN PRIVATE ${CUDA_INCLUDE_DIRS})
    target_link_libraries(cudaBSSN bssn_cuda m)
    target_link_libraries(cudaBSSN ${CUDA_LIBS})
endif()
