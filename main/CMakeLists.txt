file(GLOB MY_MAIN "*.h" "*.cpp")
ADD_EXECUTABLE(heat_equations ${MY_MAIN})
TARGET_LINK_LIBRARIES(heat_equations numericallib)