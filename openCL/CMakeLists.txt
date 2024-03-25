cmake_minimum_required(VERSION 2.6)

#INCLUDE(FindOpenMP)

PROJECT(nbcl)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./)

ADD_EXECUTABLE(nbcl main.cpp nbody.cpp gen.cpp sphere.cpp parse_args.cpp clsetup.cpp glsupport.cpp)
target_link_libraries(nbcl gsl gslcblas cfitsio GL GLU glut GLEW OpenCL)
