file(GLOB MY_NMERIC "*.h" "*.cpp")
ADD_LIBRARY(numericallib STATIC ${MY_NMERIC})
target_include_directories(numericallib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )