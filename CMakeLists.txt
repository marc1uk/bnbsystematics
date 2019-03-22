set(CMAKE_C_COMPILER "$ENV{GCC_FQ_DIR}/bin/gcc")
set(CMAKE_CXX_COMPILER "$ENV{GCC_FQ_DIR}/bin/g++")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -std=c++1y -g -rdynamic " )
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{G4LIB} $ENV{ROOTSYS}/etc/cmake)

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(BNBSys)

#----------------------------------------------------------------------------
# Setup root
#
find_package(ROOT REQUIRED)
find_package(CLHEP REQUIRED) 

add_subdirectory(Base)
add_subdirectory(Calculators)
add_subdirectory(src)
