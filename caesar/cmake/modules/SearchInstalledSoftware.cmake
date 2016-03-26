#---Check for installed packages depending on the build options/components eamnbled -
#include(ExternalProject)
#include(FindPackageHandleStandardArgs)

#-------------------------
# -- Check for ROOT-------
#-------------------------
message(STATUS "Looking for ROOT")

if(NOT DEFINED ENV{ROOTSYS})
	message(SEND_ERROR "ROOTSYS variable not defined!")
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake)
list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/cmake/modules)

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED MODULE COMPONENTS MathCore RIO Hist Tree Net FITSIO PyROOT RInterface)
#find_package(ROOT REQUIRED)
#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE}) ### NOT WORKING!!
#include(ROOTUseFile) 
#include_directories(${ROOT_INCLUDE_DIR})

message (STATUS "ROOT HEADERS: ${ROOT_INCLUDE_DIRS}, LIBS: ${ROOT_LIBRARIES}")
#--------------------------------

#-------------------------
# -- Check for BOOST -----
#-------------------------
find_package(Boost REQUIRED COMPONENTS filesystem system)
message (STATUS "BOOST HEADERS: ${Boost_INCLUDE_DIRS}, LIBS: ${Boost_LIBRARIES}")
#-------------------------

#-------------------------
# -- Check for OPENCV -----
#-------------------------
find_package(OpenCV REQUIRED)
message (STATUS "OPENCV HEADERS: ${OpenCV_INCLUDE_DIRS}, LIBS: ${OpenCV_LIBS}, ${OpenCV_LIB_COMPONENTS}")
find_package(OpenCV2 REQUIRED)
message (STATUS "OPENCV2 HEADERS: ${OpenCV2_INCLUDE_DIRS}, LIBS: ${OpenCV2_LIBRARIES}")
#-------------------------

#-----------------------------
# -- Check for R PROJECT -----
#-----------------------------
message (STATUS "Looking for R")
find_package (R REQUIRED COMPONENTS base RInside Rcpp rrcovHD truncnorm FNN akima)
message (STATUS "R_INCLUDE_DIR: ${R_INCLUDE_DIR}")
message (STATUS "R_LIBRARIES: ${R_LIBRARIES}")
#-------------------------


#-----------------------------
# -- Check for Python    -----
#-----------------------------
message (STATUS "Looking for Python")
find_package (PythonLibs REQUIRED)
message (STATUS "PYTHON_LIBRARIES: ${PYTHON_LIBRARIES}, PYTHON_INCLUDE_DIRS: ${PYTHON_INCLUDE_DIRS}")



