cmake_minimum_required(VERSION 3.1)

file(GLOB MY_NMERIC "*.h" "*.cpp")
ADD_LIBRARY(femlib STATIC ${MY_NMERIC})
include_directories(${EIGEN3_INCLUDE_DIR})
set(CMAKE_CXX_STANDARD 11)


target_include_directories(femlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(femlib PUBLIC pdesolverlib numericallib)
add_dependencies(femlib eigen)