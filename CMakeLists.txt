cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PDAQ)
include_directories(${MyProject_SOURCE_DIR})
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - RcmaOOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS RIO Net)

#message(${ROOT_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIRS})


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
#include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#---Create a shared library with geneated dictionary
#add_library(re_tab re_tab.cpp re_tab_main.cpp)
#add_library(ntuple SHARED ntuple.cpp H2E.cpp )
add_library(PDAQ SHARED MParManager.cc MPar.cc MParContainer.cc FTGeo.cc)
target_link_libraries(PDAQ ${ROOT_LIBRARIES})

#---Create  a main program using the library

add_executable(Main main.cpp)
target_link_libraries(Main PDAQ)

