cmake_minimum_required(VERSION 3.1)

file(GLOB MY_NMERIC "*.h" "*.cpp")
ADD_LIBRARY(pdesolverlib STATIC ${MY_NMERIC})


target_include_directories(pdesolverlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_dependencies(pdesolverlib eigen)
target_link_libraries(pdesolverlib numericallib)