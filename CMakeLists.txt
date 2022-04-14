cmake_minimum_required(VERSION 3.4)
project(pyPDFs)



# Fortran
enable_language(Fortran)
#set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -g -fPIC")	# -O0 is required to produce grids for SIH, otherwise -O3 should work
set(CMAKE_Fortran_MODULE_DIRECTORY "./include/mod")

# C++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -fpermissive -Wno-return-type -fPIC ")
set(CMAKE_CXX_COMPILER g++)
#set(CMAKE_BUILD_TYPE DEBUG)


file(GLOB_RECURSE LIBF123L_SOURCE_FILES "Theory/f123/src/*.f")

add_library(f123L SHARED ${LIBF123L_SOURCE_FILES})
target_link_libraries(f123L)
set_target_properties(f123L PROPERTIES LIBRARY_OUTPUT_DIRECTORY ./libs)


