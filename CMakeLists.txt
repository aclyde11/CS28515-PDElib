cmake_minimum_required(VERSION 2.8)
project(CS28515Proj1)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ADD_SUBDIRECTORY(main)
ADD_SUBDIRECTORY(libnumerical)

option (enable-test "Compile all project unit tests" ON)

if (enable-test)
    ADD_SUBDIRECTORY(test)
    FIND_PACKAGE(GTest)
endif()
