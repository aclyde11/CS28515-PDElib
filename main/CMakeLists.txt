file(GLOB MY_MAIN "*.h" "*.cpp")
ADD_EXECUTABLE(main ${MY_MAIN})
TARGET_LINK_LIBRARIES(main numericallib pdesolverlib)