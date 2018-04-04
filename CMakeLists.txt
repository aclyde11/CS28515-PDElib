cmake_minimum_required(VERSION 3.7)
project(CS28515Proj1)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp utility.cpp utility.h numvec.cpp numvec.h tridiag.cpp tridiag.h)
add_executable(CS28515Proj1 ${SOURCE_FILES})