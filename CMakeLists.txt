cmake_minimum_required(VERSION 3.0)
project(CS28515Proj1)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
include(eigen)
include_directories(${EIGEN_INCLUDE_DIRS})
include_directories(./build/eigen/src/eigen)

set(CMAKE_CXX_STANDARD 11)

add_definitions(-Wall -std=c++0x)


SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(libnumerical)
ADD_SUBDIRECTORY(libpdesolvers)
ADD_SUBDIRECTORY(libfem)

option(enable-test "Compile all project unit tests" ON)

SET(enable-test OFF)
if (enable-test)
    ADD_SUBDIRECTORY(test)
    FIND_PACKAGE(GTest)
endif ()


# Run test cases:

# Run python graph to make pictures.