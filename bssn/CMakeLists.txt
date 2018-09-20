cmake_minimum_required(VERSION 2.8)
project(bssn) 



# CUDA PACKAGE
find_package(CUDA REQUIRED)

if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
    SET(ALL_CUDA_LIBS ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY})
    SET(LIBS ${LIBS} ${ALL_CUDA_LIBS})
    message(STATUS "CUDA_LIBRARIES: ${CUDA_INCLUDE_DIRS} ${ALL_CUDA_LIBS}")
    set(CUDA_PROPAGATE_HOST_FLAGS ON)
    set(CUDA_SEPARABLE_COMPILATION OFF)

    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)

    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_60,code=compute_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_61,code=compute_61)
    list( APPEND CUDA_NVCC_FLAGS -use_fast_math )
    list( APPEND CUDA_NVCC_FLAGS -lineinfo )

    # list( APPEND CUDA_NVCC_FLAGS --maxrregcount 16 )
    # list( APPEND CUDA_NVCC_FLAGS --ptxas-options=-v )

endif()


find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

include_directories(${PROJECT_BINARY_DIR}
        include
        utils
        cuda_gr/utils
        )

set(INC_FILES
        include/computeBSSN.h
        include/derivs.h
        include/rhs.h
        include/def.h
        include/bssn_profiler.h
        include/profile_param.h
        include/utils.h

        include/rhs_cuda.cuh
        include/derivs_cuda.cuh
        include/bssneqn_solve.cuh
        include/cuda_bssneqs.cuh

        cuda_gr/include/dataGeneration.h
        cuda_gr/include/block.h
        cuda_gr/include/cudaBSSN.h
        cuda_gr/include/rhsMethods.h
        # cuda_gr/include/cuda_rhs.cuh
        cuda_gr/include/mergeSort.h
        # cuda_gr/include/intermediateVariables.h
        )

set(SRC_FILES
        src/derivs.cpp
        src/bssn_profiler.cpp
        src/rhs.cpp
        src/profile_param.cpp
        src/utils.cpp

        cuda_src/bssneqn_solve.cu
        cuda_src/device_derivs.cu
        cuda_src/rhs_cuda.cu

        cuda_gr/src/dataGeneration.cpp
        cuda_gr/src/cudaBSSN.cpp
        cuda_gr/src/rhsMethods.cpp
        # cuda_gr/src/cuda_rhs.cu
        cuda_gr/src/mergeSort.cpp
        # cuda_gr/src/intermediateVariables.cpp
        )


add_executable(computeBSSN ${INC_FILES} ${SRC_FILES})
target_link_libraries(computeBSSN ${LIBS} )
target_include_directories(computeBSSN PRIVATE cuda_gr/include)
