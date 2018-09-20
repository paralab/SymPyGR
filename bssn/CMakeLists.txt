cmake_minimum_required(VERSION 2.8)
project(bssn) 



# CUDA PACKAGE
find_package(CUDA REQUIRED)

if(CUDA_FOUND)
    include_directories(${CUDA_INCLUDE_DIRS})
    SET(ALL_CUDA_LIBS ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY})
    SET(LIBS ${LIBS} ${ALL_CUDA_LIBS})
    message(STATUS "CUDA_LIBRARIES: ${CUDA_INCLUDE_DIRS} ${ALL_CUDA_LIBS}")
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(CUDA_SEPARABLE_COMPILATION ON)

    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_30,code=sm_30)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_50,code=sm_50)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_60,code=sm_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode arch=compute_61,code=sm_61)

    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_60,code=compute_60)
    list( APPEND CUDA_NVCC_FLAGS -gencode=arch=compute_61,code=compute_61)
    list( APPEND CUDA_NVCC_FLAGS -use_fast_math )
    list( APPEND CUDA_NVCC_FLAGS -lineinfo )
    list( APPEND CUDA_NVCC_FLAGS -std=c++11)

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
        cuda_gr/cuda_include
        )

set(INC_FILES
        include/computeBSSN.h
        include/derivs.h
        include/rhs.h
        include/def.h
        include/bssn_profiler.h
        include/profile_param.h
        include/utils.h

        # cuda_gr/include/dataGeneration.h
        # cuda_gr/include/block.h
        # cuda_gr/include/cudaBSSN.h
        # cuda_gr/include/rhsMethods.h
        # cuda_gr/include/cuda_rhs.cuh
        # cuda_gr/include/mergeSort.h
        # cuda_gr/include/intermediateVariables.h
        )

set(SRC_FILES
        src/derivs.cpp
        src/bssn_profiler.cpp
        src/rhs.cpp
        src/profile_param.cpp
        src/utils.cpp
        src/computeBSSN.cpp

        # cuda_src/bssneqn_solve.cu
        # cuda_src/device_derivs.cu
        # cuda_src/rhs_cuda.cu
        
        # cuda_gr/src/dataGeneration.cpp
        # cuda_gr/src/cudaBSSN.cpp
        # cuda_gr/src/rhsMethods.cpp
        # cuda_gr/src/cuda_rhs.cu
        # cuda_gr/src/mergeSort.cpp
        # cuda_gr/src/intermediateVariables.cpp
        )


# set(CUDA_GR_INC
        # cuda_gr/cuda_include/cudaRHS.cuh
        # cuda_gr/cuda_include/cudaDerivs.cuh
        # cuda_gr/cuda_include/deviceDerivs.cuh
        # cuda_gr/cuda_include/bssnEqns.cuh
        # cuda_gr/cuda_include/getOutput.cuh
        # cuda_gr/cuda_include/cudaBCS.cuh
        
        # )

# set(CUDA_GR_SRC
        # cuda_gr/cuda_src/cudaRHS.cu
        # cuda_gr/cuda_src/cudaDerivs.cu
        # cuda_gr/cuda_src/deviceDerivs.cu
        # cuda_gr/cuda_src/bssnEqns.cu
        # cuda_gr/cuda_src/getOutput.cu
        # cuda_gr/cuda_src/cudaBCS.cu
        # )

# cuda_add_library(bssn_cuda ${CUDA_GR_INC} ${CUDA_GR_SRC})
# target_include_directories(bssn_cuda PRIVATE cuda_gr/cuda_include)
# target_include_directories(bssn_cuda PRIVATE cuda_gr/include)
# target_include_directories(bssn_cuda PRIVATE include)


# add_executable(computeBSSN ${CUDA_GR_INCLUDE_H} ${CUDA_GR_SRC_CPP} ${INC_FILES} ${SRC_FILES})
add_executable(computeBSSN ${INC_FILES} ${SRC_FILES})
# target_link_libraries(computeBSSN bssn_cuda ${LIBS} )
# target_include_directories(bssn_cuda PRIVATE cuda_gr/cuda_include)
target_include_directories(computeBSSN PRIVATE cuda_gr/include)
