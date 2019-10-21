cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(PDAQ)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


include_directories(${MyProject_SOURCE_DIR})
# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - RcmaOOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS})
#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS Core Gpad Hist Physics Matrix RIO Tree Cint Graf)

#include(RootNewMacros)
include_directories(${ROOT_INCLUDE_DIRS})


#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
message(STATUS "Root Use File ${ROOT_USE_FILE}")
include(${ROOT_USE_FILE})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#---Create a shared library with geneated dictionary
#add_library(re_tab re_tab.cpp re_tab_main.cpp)
#add_library(ntuple SHARED ntuple.cpp H2E.cpp )

set(HDRS MParManager.h MPar.h MParContainer.h FTGeo.h SttRawHit.h SttEvent.h EmcEvent.h SttDetector.h panda_subsystem.h panda_subsystem_sb.h panda_subsystem_emc.h panda_subsystem_stt.h EmcHit.h SttHit.h panda_stt_cal.h Stt_Cal_Event.h SttTrackHit.h SttTrackEvent.h panda_stt_track.h
MLookup.h MLookupContainer.h MLookupManager.h FTdef.h panda_subsystem_sci.h SciEvent.h SciHit.h PDAQ_Trackfilter.h)

set(SRCS MParManager.cc MPar.cc MParContainer.cc FTGeo.cc  SttRawHit.cc SttEvent.cc EmcEvent.cc SttDetector.cc panda_subsystem.cc panda_subsystem_sb.cc panda_subsystem_emc.cc panda_subsystem_stt.cc EmcHit.cc SttHit.cc  panda_stt_cal.cc Stt_Cal_Event.cc
SttTrackHit.cc SttTrackEvent.cc panda_stt_track.cc MLookup.cc MLookupContainer.cc MLookupManager.cc panda_subsystem_sci.cc SciEvent.cc SciHit.cc PDAQ_Trackfilter.cc)

set(CMAKE_INSTALL_LIBDIR "lib")
ROOT_GENERATE_DICTIONARY(G__PDAQ_cc
	 ${HDRS}
	 LINKDEF LinkDef.h)
#ROOT_GENERATE_ROOTMAP(PDAQ LINKDEF LinkDef.h)
#add_custom_command(
	#OUTPUT rootmap
	#COMMAND rlibmap -o libPDAQ.rootmap -l libPDAQ -c ../LinkDef.h
	#DEPENDS LinkDef.h rlibmap)
#add_custom_target(libPDAQ.rootmap ALL DEPENDS rootmap PDAQ)

add_library(PDAQ SHARED ${SRCS} G__PDAQ_cc.cxx)

#file(GLOB HDRS "*.h")


target_link_libraries(PDAQ ${ROOT_LIBRARIES})

#---Create  a main program using the library

add_executable(pdaq_unpacker PDAQ_RawDecoder_EMC_STT.cc)
target_link_libraries(pdaq_unpacker PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_unpacker_HADES PDAQ_RawDecoder_HADES.cc)
target_link_libraries(pdaq_unpacker_HADES PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_calibrater pdaq_calibirator.cc PDAQ_Stt_Calibirator.C)
target_link_libraries(pdaq_calibrater PDAQ ${ROOT_LIBRARIES})

#add_executable(pdaq_cluster_finder PDAQ_Cluster_Finder.cc)
#target_link_libraries(pdaq_cluster_finder PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_cluster_finder_cosy PDAQ_Cluster_Finder_Cosy.cc)
target_link_libraries(pdaq_cluster_finder_cosy PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_drift_cal PDAQ_Drift_Cal.cc)
target_link_libraries(pdaq_drift_cal PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_rad_cal PDAQ_Rad_Cal.cc)
target_link_libraries(pdaq_rad_cal PDAQ ${ROOT_LIBRARIES})

add_executable(pdaq_spatial_resolution PDAQ_Spl_Res.cc)
target_link_libraries(pdaq_spatial_resolution PDAQ ${ROOT_LIBRARIES})
