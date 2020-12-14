# CMakeLists.txt for FAS package. It creates a library with dictionary and a main program
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(PicoDst)

# You need to tell CMake where to find the ROOT installation. This can be done
# in a number of ways:
#   - ROOT built with classic configure/make use the provided
#   $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix
#   for ROOT

find_package(Git)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

set(FAIRSOFT_PATH $ENV{FAIRROOTPATH}/share/fairbase/cmake/modules)
list(APPEND CMAKE_MODULE_PATH ${FAIRSOFT_PATH})
find_package(FairRoot QUIET COMPONENTS Base)

if(NOT FAIRROOT_FOUND)
  message(STATUS "FairRoot package NOT found.")
endif()

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})


#set (CMAKE_CXX_STANDARD 11)
add_definitions(${ROOT_CXX_FLAGS})

set(CMAKE_BUILD_TYPE Debug)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -Wall -pthread")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -Wall")

set(PICO_DST_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/converter"
  "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOT_INLCUDE_DIRS}
)

include_directories(${PICO_DST_INCLUDE_DIRECTORIES})

set(PICO_DST_INCLUDE_LIBRARIES
  ${ROOT_LIBRARIES}
)

set(PicoDst_h_files
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstBaseEvent.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstMCEvent.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstRecoEvent.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstBaseTrack.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstMCTrack.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstRecoTrack.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstFHCal.h"
)

set(PicoDst_cxx_files
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstBaseEvent.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstMCEvent.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstRecoEvent.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstBaseTrack.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstMCTrack.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstRecoTrack.cxx"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/PicoDstFHCal.cxx"
)

set(PicoDst_LinkDef
	"${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDst.LinkDef.h"
)

#---Generate dictionary
ROOT_GENERATE_DICTIONARY(G__PicoDst
  ${PicoDst_h_files}
  LINKDEF ${PicoDst_LinkDef}
)

#---Compile library
add_library(PicoDst SHARED ${PicoDst_cxx_files} G__PicoDst.cxx)
target_link_libraries(PicoDst ${PICO_DST_INCLUDE_LIBRARIES})

set(ENV{PICO_DST_BIN} "${CMAKE_CURRENT_BINARY_DIR}")
set(ENV{PICO_DST_SRC} "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(ENV{PICO_DST_INC} "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(PICO_DST_BIN "${CMAKE_CURRENT_BINARY_DIR}")
set(PICO_DST_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(PICO_DST_INC "${CMAKE_CURRENT_SOURCE_DIR}/include")
message(STATUS "PICO_DST_BIN: " $ENV{PICO_DST_BIN})
message(STATUS "PICO_DST_INC: " $ENV{PICO_DST_INC})
message(STATUS "PICO_DST_SRC: " $ENV{PICO_DST_SRC})
set(BASH $ENV{SHELL})
file( WRITE  ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "#!${BASH}\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export PICO_DST_BIN=$ENV{PICO_DST_BIN}\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export PICO_DST_SRC=$ENV{PICO_DST_SRC}\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export PICO_DST_INC=$ENV{PICO_DST_INC}\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export LD_LIBRARY_PATH=$PICO_DST_BIN:$LD_LIBRARY_PATH\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export PATH=$PICO_DST_BIN:$PATH\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "export CPATH=$PICO_DST_INC:$CPATH\n")


# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  RESULT_VARIABLE GIT_BRANCH_ERROR_CODE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if( GIT_BRANCH_ERROR_CODE )
  set(GIT_BRANCH 0)
endif()

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  RESULT_VARIABLE GIT_COMMIT_HASH_ERROR_CODE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if( GIT_COMMIT_HASH_ERROR_CODE )
  set(GIT_COMMIT_HASH 0)
endif()

# Get the current version
execute_process(
  COMMAND git describe --tags --match "*"
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_TAG_VERSION
  RESULT_VARIABLE GIT_TAG_VERSION_ERROR_CODE
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if( GIT_TAG_VERSION_ERROR_CODE )
  set(GIT_TAG_VERSION 0.0)
endif()

#---If one wants to input those definitions manually, use this commands:
#add_definitions("-DGIT_COMMIT_HASH=${GIT_COMMIT_HASH}")
#add_definitions("-DGIT_BRANCH=${GIT_BRANCH}")
#add_definitions("-DGIT_TAG_VERSION=${GIT_TAG_VERSION}")
#add_definitions("-DPICO_DST_BIN_DIR=${PICO_DST_BIN}")
#add_definitions("-DPICO_DST_SRC_DIR=${PICO_DST_SRC}")
#add_definitions("-DPICO_DST_INC_DIR=${PICO_DST_INC}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/include/PicoDstVersion.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/generated/PicoDstVersion.h
)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/generated)

set_target_properties (PicoDst 
                        PROPERTIES 
                            VERSION ${GIT_TAG_VERSION} SOVERSION ${GIT_BRANCH}-${GIT_COMMIT_HASH}
			    PUBLIC_HEADER ${CMAKE_CURRENT_BINARY_DIR}/generated/PicoDstVersion.h 
)

add_executable(PicoDst-config "${CMAKE_CURRENT_SOURCE_DIR}/bin/PicoDstConfig.cpp")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "\n")
file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "alias PicoDst-config='${PICO_DST_BIN}/PicoDst-config'\n")

if( FAIRROOT_FOUND )
#---Compile executable
if(DEFINED ENV{VMCWORKDIR})
set(MPDROOTDIR $ENV{VMCWORKDIR})
message(STATUS "MpdRoot environment is found: ${MPDROOTDIR}")
set(PICODSTCONVERTER_INCLUDE_DIRECTORIES
  ${FAIRROOT_INCLUDE_DIR}
  ${MPDROOTDIR}/tpc
  ${MPDROOTDIR}/kalman
  ${MPDROOTDIR}/lhetrack
  ${MPDROOTDIR}/mcstack
  ${MPDROOTDIR}/generators
  ${MPDROOTDIR}/strawendcap
  ${MPDROOTDIR}/etof
  ${MPDROOTDIR}/tof
  ${MPDROOTDIR}/sft
  ${MPDROOTDIR}/sts
  ${MPDROOTDIR}/zdc
  ${MPDROOTDIR}/mpddst
  ${MPDROOTDIR}/mpdbase
  ${MPDROOTDIR}/mpdpid
)

#---Add Needed libraries
find_library(BASE NAMES Base PATHS ${FAIRROOT_LIBRARY_DIR})
find_library(MCSTACK NAMES MCStack PATHS $ENV{FAIRLIBDIR})
if (MCSTACK)
  message(STATUS "MCStack library is found")
  add_definitions("-D_OLD_MCSTACK_")
endif()
if (NOT MCSTACK)
  find_library(MPDMCSTACK NAMES MpdMCStack PATHS $ENV{FAIRLIBDIR})
  message(STATUS "MpdMCStack library is found")
  add_definitions("-D_NEW_MCSTACK_")
endif()
find_library(PASSIVE NAMES Passive PATHS $ENV{FAIRLIBDIR})
find_library(MPDFIELD NAMES MpdField PATHS $ENV{FAIRLIBDIR})
find_library(MPDBASE NAMES MpdBase PATHS $ENV{FAIRLIBDIR})
find_library(MPDDST NAMES MpdDst PATHS $ENV{FAIRLIBDIR})
find_library(MPDGENERALGENERATOR NAMES MpdGeneralGenerator PATHS $ENV{FAIRLIBDIR})
find_library(KALMAN NAMES Kalman PATHS $ENV{FAIRLIBDIR})
find_library(CLUSTER NAMES Cluster PATHS $ENV{FAIRLIBDIR})
find_library(TPC NAMES tpc PATHS $ENV{FAIRLIBDIR})
find_library(TOF NAMES Tof PATHS $ENV{FAIRLIBDIR})
find_library(ETOF NAMES Etof PATHS $ENV{FAIRLIBDIR})
find_library(STRAWENDCAP NAMES Strawendcap PATHS $ENV{FAIRLIBDIR})
find_library(ZDC NAMES Zdc PATHS $ENV{FAIRLIBDIR})
find_library(FFD NAMES Ffd PATHS $ENV{FAIRLIBDIR})
find_library(CPC NAMES Cpc PATHS $ENV{FAIRLIBDIR})
find_library(STRAWECT NAMES StrawECT PATHS $ENV{FAIRLIBDIR})
find_library(LHETrack NAMES LHETrack PATHS $ENV{FAIRLIBDIR})
find_library(STS NAMES Sts PATHS $ENV{FAIRLIBDIR})
find_library(MPDPID NAMES MpdPid PATHS $ENV{FAIRLIBDIR})

if(MCSTACK)
 set(PICODSTCONVERTER_INCLUDE_LIBRARIES
   PicoDst ${BASE}
   ${MCSTACK} ${PASSIVE} ${KALMAN} ${CLUSTER} ${TPC} ${STS} ${ZDC} ${FFD}
   ${STRAWECT} ${STRAWENDCAP} ${LHETRACK} ${TOF}
   ${MPDFIELD} ${MPDBASE} ${MPDDST} ${MPDGENERALGENERATOR} ${MPDPID}
 )
elseif(MPDMCSTACK)
 set(PICODSTCONVERTER_INCLUDE_LIBRARIES
   PicoDst ${BASE}
   ${MPDMCSTACK} ${PASSIVE} ${KALMAN} ${CLUSTER} ${TPC} ${STS} ${ZDC} ${FFD}
   ${STRAWECT} ${STRAWENDCAP} ${LHETRACK} ${TOF}
   ${MPDFIELD} ${MPDBASE} ${MPDDST} ${MPDGENERALGENERATOR} ${MPDPID}
 )
elseif()
 message(FATAL_ERROR "No MCSTACK nor MPDMCSTACK libraries was found!")
endif()

add_executable(PicoDstConverter "${CMAKE_CURRENT_SOURCE_DIR}/converter/PicoDstConverter.cpp")
target_include_directories(PicoDstConverter PUBLIC ${PICODSTCONVERTER_INCLUDE_DIRECTORIES})
target_link_libraries(PicoDstConverter ${PICODSTCONVERTER_INCLUDE_LIBRARIES})

file( APPEND ${CMAKE_CURRENT_BINARY_DIR}/setPicoDst.sh "alias PicoDstConverter='${PICO_DST_BIN}/PicoDstConverter'\n")

else(DEFINED ENV{VMCWORKDIR})
message(WARNING "MpdRoot environment is NOT found. Skip corresponding step.")
endif(DEFINED ENV{VMCWORKDIR})

endif()
