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
find_package(Boost REQUIRED COMPONENTS filesystem system regex thread log log_setup)
add_definitions(-DBOOST_LOG_USE_NATIVE_SYSLOG -DBOOST_LOG_DYN_LINK)
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


#==================================
#==   Check for Log4Cxx         ===
#==================================
MESSAGE(STATUS "Looking for Log4Cxx")
FIND_PACKAGE(Log4Cxx REQUIRED)
MESSAGE(STATUS "LOG4CXX_INCLUDE_DIR: ${LOG4CXX_INCLUDE_DIRS}")
MESSAGE(STATUS "LOG4CXX_LIBRARIES: ${LOG4CXX_LIBRARIES}")


#==================================
#==   Check for JSONCPP         ===
#==================================
MESSAGE(STATUS "Looking for JSONCPP")

if(NOT DEFINED ENV{JSONCPP_ROOT})
	MESSAGE(SEND_ERROR "JSONCPP_ROOT variable not defined!")
endif()

SET (JSONCPP_ROOT $ENV{JSONCPP_ROOT})
MESSAGE(STATUS "JSONCPP_ROOT: ${JSONCPP_ROOT}")

FIND_PATH (JSONCPP_INCLUDE_DIR
	NAMES json/json.h
  HINTS
 ${JSONCPP_ROOT}/include
)

FIND_LIBRARY (JSONCPP_LIBRARIES NAMES jsoncpp HINTS ${JSONCPP_ROOT}/lib)

MARK_AS_ADVANCED (JSONCPP_INCLUDE_DIR JSONCPP_LIBRARIES)
MESSAGE(STATUS "JSONCPP_INCLUDE_DIR: ${JSONCPP_INCLUDE_DIR}")
MESSAGE(STATUS "JSONCPP_LIBRARIES: ${JSONCPP_LIBRARIES}")


#===========================================
#==   Check for Tango Framework          ===
#===========================================
option(ENABLE_SERVER "Enable Caesar server in building" OFF)
if(ENABLE_SERVER)
	## Define a preprocessor option
	add_definitions(-DBUILD_CAESAR_SERVER)
	

	#================================
	#==   Check for MSGPACK       ===
	#================================
	MESSAGE(STATUS "Looking for MessagePack lib")
	FIND_PACKAGE(MsgPack REQUIRED)
	IF (NOT MSGPACK_FOUND)
		MESSAGE(SEND_ERROR "MessagePack not found!")
	endif()
	MESSAGE(STATUS "MSGPACK_INCLUDE_DIR: ${MSGPACK_INCLUDE_DIRS}")
	MESSAGE(STATUS "MSGPACK_LIBRARIES: ${MSGPACK_LIBRARIES}")
	

	#================================
	#==   Check for OMNIORB       ===
	#================================
	MESSAGE(STATUS "Looking for omniORB")
	if(NOT DEFINED ENV{OMNI_ROOT})
		MESSAGE(SEND_ERROR "OMNI_ROOT variable not defined!")
	endif()

	SET (OMNI_ROOT $ENV{OMNI_ROOT})
	MESSAGE(STATUS "OMNI_ROOT: ${OMNI_ROOT}")

	FIND_PATH (OMNIORB_INCLUDE_DIR
		NAMES omniconfig.h
  	HINTS
  	${OMNI_ROOT}/include
	)

	FIND_LIBRARY (OMNIORB_LIB1 NAMES omniORB4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB2 NAMES COS4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB3 NAMES omniDynamic4 HINTS ${OMNI_ROOT}/lib)
	FIND_LIBRARY (OMNIORB_LIB4 NAMES omnithread HINTS ${OMNI_ROOT}/lib)
	list(APPEND OMNIORB_LIBRARIES ${OMNIORB_LIB1} ${OMNIORB_LIB2} ${OMNIORB_LIB3} ${OMNIORB_LIB4})

	MARK_AS_ADVANCED (OMNIORB_INCLUDE_DIR OMNIORB_LIBRARIES)
	MESSAGE(STATUS "OMNIORB_INCLUDE_DIR: ${OMNIORB_INCLUDE_DIR}")
	MESSAGE(STATUS "OMNIORB_LIBRARIES: ${OMNIORB_LIBRARIES}")

	#================================
	#==   Check for ZMQ           ===
	#================================
	MESSAGE(STATUS "Looking for ZMQ")
	if(NOT DEFINED ENV{ZMQ_ROOT})
		MESSAGE(SEND_ERROR "ZMQ_ROOT variable not defined!")
	endif()

	SET (ZMQ_ROOT $ENV{ZMQ_ROOT})
	MESSAGE(STATUS "ZMQ_ROOT: ${ZMQ_ROOT}")

	FIND_PATH (ZMQ_INCLUDE_DIR
		NAMES zmq.h
  	HINTS
  	${ZMQ_ROOT}/include
	)

	FIND_LIBRARY (ZMQ_LIBRARIES NAMES zmq HINTS ${ZMQ_ROOT}/lib)

	MARK_AS_ADVANCED (ZMQ_INCLUDE_DIR ZMQ_LIBRARIES)
	MESSAGE(STATUS "ZMQ_INCLUDE_DIR: ${ZMQ_INCLUDE_DIR}")
	MESSAGE(STATUS "ZMQ_LIBRARIES: ${ZMQ_LIBRARIES}")


	#==================================	
	#==   Check for TANGO           ===
	#==================================
	MESSAGE(STATUS "Looking for TANGO Framework")
	
	## Define Tango preprocessor flag
  add_definitions(-DUSE_TANGO)
	add_definitions(-DAPPENDERS_HAVE_LEVEL_THRESHOLD=1)

	if(NOT DEFINED ENV{TANGO_ROOT})
		MESSAGE(SEND_ERROR "TANGO_ROOT variable not defined!")
	endif()

	SET (TANGO_ROOT $ENV{TANGO_ROOT})
	MESSAGE(STATUS "TANGO_ROOT: ${TANGO_ROOT}")

	FIND_PATH (TANGO_INCLUDE_DIR
		NAMES tango.h
  	HINTS
  	${TANGO_ROOT}/include/tango
	)

	FIND_LIBRARY (TANGO_LIB1 NAMES tango HINTS ${TANGO_ROOT}/lib)
	FIND_LIBRARY (TANGO_LIB2 NAMES log4tango HINTS ${TANGO_ROOT}/lib)
	list(APPEND TANGO_LIBRARIES ${TANGO_LIB1} ${TANGO_LIB2})

	MARK_AS_ADVANCED (TANGO_INCLUDE_DIR TANGO_LIBRARIES)
	MESSAGE(STATUS "TANGO_INCLUDE_DIR: ${TANGO_INCLUDE_DIR}")
	MESSAGE(STATUS "TANGO_LIBRARIES: ${TANGO_LIBRARIES}")

endif() ### close if enable Tango

