cmake_minimum_required(VERSION 2.8)
project(bssn)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


include_directories(${PROJECT_BINARY_DIR}
        include
        utils
        )


set(INC_FILES
        include/computeBSSN.h
        include/derivs.h
        include/rhs.h
        include/def.h
        include/bssn_profiler.h
        include/profile_param.h
        )

set(SRC_FILES
        src/derivs.cpp
        src/rhs.cpp
        src/bssn_profiler.cpp
        src/profile_param.cpp
        )

add_executable(computeBSSN ${INC_FILES} ${SRC_FILES} src/computeBSSN.cpp)