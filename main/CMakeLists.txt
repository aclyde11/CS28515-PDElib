file(GLOB PROJ1_FILE "proj1.cpp")
ADD_EXECUTABLE(proj1 ${PROJ1_FILE})
TARGET_LINK_LIBRARIES(proj1 numericallib pdesolverlib)

file(GLOB PROJ2_FILE "proj2.cpp")
ADD_EXECUTABLE(proj2 ${PROJ2_FILE})
TARGET_LINK_LIBRARIES(proj2 numericallib pdesolverlib)