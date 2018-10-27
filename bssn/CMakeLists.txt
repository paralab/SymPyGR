cmake_minimum_required(VERSION 2.8)
project(bssn) 

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(CUDA REQUIRED)
if(CUDA_FOUND)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(CUDA_SEPARABLE_COMPILATION ON)

    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)

    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_60,code=compute_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_61,code=compute_61)
    list( APPEND CUDA_NVCC_FLAGS -use_fast_math )
    list( APPEND CUDA_NVCC_FLAGS -std=c++11)
endif()

option(ENABLE_CUDA_TEST "enable cuda test" ON)

set(INC_FILES
        cuda_gr/include/block.h
        cuda_gr/include/cudaBSSN.h
        cuda_gr/include/dataGeneration.h
        cuda_gr/include/mergeSort.h
        cuda_gr/include/profileGPU.h
        cuda_gr/include/profilerClass.h
        cuda_gr/include/rhsMethods.h
        )

set(SRC_FILES
        cuda_gr/src/cudaBSSN.cpp
        cuda_gr/src/dataGeneration.cpp
        cuda_gr/src/mergeSort.cpp
        cuda_gr/src/profileGPU.cpp
        cuda_gr/src/profilerClass.cpp
        cuda_gr/src/rhsMethods.cpp
        )

set(OTHER
        include/def.h
        include/utils.h
        src/utils.cpp
        )

set(CUDA_GR_INC
        cuda_gr/cuda_include/bssnEqns.cuh
        cuda_gr/cuda_include/cudaBCS.cuh
        cuda_gr/cuda_include/cudaDerivs.cuh
        cuda_gr/cuda_include/cudaRHS.cuh
        cuda_gr/cuda_include/deviceDerivs.cuh
        cuda_gr/cuda_include/getOutput.cuh
        )

set(CUDA_GR_SRC
        cuda_gr/cuda_src/bssnEqns.cu
        cuda_gr/cuda_src/cudaBCS.cu
        cuda_gr/cuda_src/cudaDerivs.cu
        cuda_gr/cuda_src/cudaRHS.cu
        cuda_gr/cuda_src/deviceDerivs.cu
        cuda_gr/cuda_src/deviceDerivsShared.cu
        cuda_gr/cuda_src/getOutput.cu
        )

set(SOURCE_FILES 
        ${INC_FILES}
        ${SRC_FILES}
        ${OTHER}
        )

if(ENABLE_CUDA_TEST)
        add_definitions(-DENABLE_CUDA_TEST)
        set(INC 
                include/bssn_profiler.h
                include/computeBSSN.h
                include/derivs.h
                include/profile_param.h
                include/rhs.h
                )

        set(SRC 
                src/bssn_profiler.cpp
                src/derivs.cpp
                src/profile_param.cpp
                src/rhs.cpp
                )

        set(SOURCE_FILES 
                ${INC_FILES}
                ${SRC_FILES}
                ${OTHER}
                ${SRC}
                ${INC}
                )
endif()

cuda_add_library(cuda_lib ${CUDA_GR_SRC} ${CUDA_GR_INC})
target_include_directories(cuda_lib PRIVATE cuda_gr/cuda_include)
target_include_directories(cuda_lib PRIVATE cuda_gr/include)
target_include_directories(cuda_lib PRIVATE cuda_gr/utils)
target_include_directories(cuda_lib PRIVATE include)

add_executable(cudaBSSN ${SOURCE_FILES})
target_include_directories(cudaBSSN PRIVATE cuda_gr/cuda_include)
target_include_directories(cudaBSSN PRIVATE cuda_gr/include)
target_include_directories(cudaBSSN PRIVATE cuda_gr/utils)
target_include_directories(cudaBSSN PRIVATE include)
target_include_directories(cudaBSSN PRIVATE utils)

target_link_libraries(cudaBSSN cuda_lib)