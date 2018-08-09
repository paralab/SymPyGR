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
        )

set(INC_FILES
        include/computeBSSN.h
        include/derivs.h
        include/rhs.h
        include/def.h
        include/derivs_cuda.h
        include/rhs_cuda.h
        include/bssn_profiler.h
        include/profile_param.h
	    include/bssneqn_solve.h
        include/cuda_bssneqs.h
        include/test_support.h
        include/utils.h
        )

set(SRC_FILES
        src/derivs.cpp
	    src/computeBSSN.cpp
        src/bssn_profiler.cpp
	    src/rhs.cpp
        src/profile_param.cpp
        src/test_support.cpp
        src/utils.cpp
        )

if(CUDA_FOUND)
    #compile cuda files and add the compiled object files to your normal source files
    FILE(GLOB_RECURSE CUDA_SRC  *.cu)
    cuda_compile(cuda_objs ${CUDA_SRC})
    SET(SRC ${SRC} ${cuda_objs})
endif()

add_executable(computeBSSN ${INC_FILES} ${SRC_FILES} ${SRC} src/computeBSSN.cpp)
target_link_libraries(computeBSSN ${LIBS} )
