cmake_minimum_required(VERSION 3.1)
#project(VBL LANGUAGES CXX VERSION 0.0.2 DESCRIPTION "virtual biology lab")
project(VBL LANGUAGES CXX)

#this is needed for proper hdf5 recognizion
ENABLE_LANGUAGE(C)

# I try to follow the tutorial at:
# https://cmake.org/cmake-tutorial/
set(VBL_VERSION_MAJOR 4)
set(VBL_VERSION_MINOR 0)
# link against tumorcode?
option(USE_TUMORCODE
       "Enable oxygen interaction with tumorcode"
       ON)
if(USE_TUMORCODE)
  set(VBL_USE_TUMORCODE 1)
#   find_package(HDF5 REQUIRED COMPONENTS CXX)
#   message(STATUS "HDF5_CXX_LIBRARIES      = ${HDF5_CXX_LIBRARIES}")
#   message(STATUS "HDF5_C_LIBRARIES    	= ${HDF5_C_LIBRARIES}")
#   message(STATUS "HDF5_LIBRARIES    	= ${HDF5_LIBRARIES}")
#   message(STATUS "HDF5_LIBRARY_DIRS	= ${HDF5_LIBRARY_DIRS}")
#   message(STATUS "HDF5_INCLUDE_DIRS	= ${HDF5_INCLUDE_DIRS}")
#   include_directories(${HDF5_INCLUDE_DIRS})
else()
  set(VBL_USE_TUMORCODE 0)
endif()

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/include/vbl_build_info.h.in"
  "${PROJECT_BINARY_DIR}/include/vbl_build_info.h"
  )
 
# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")
# c++ 11 for newer cmake
set(CMAKE_CXX_STANDARD 11)
message("compiler identified: ${CMAKE_CXX_COMPILER_VERSION}")
# if (${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
#   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gcc-name=gcc-4.9.2 -std=c++0x")
# endif()

#### chose Release a default
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()

###############  set up spack stuff  ##############

if(FALSE)
# enable @rpath in the install name for any shared library being built
# note: it is planned that a future version of CMake will enable this by default
set(CMAKE_MACOSX_RPATH 1)

# Always use full RPATH
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# http://www.kitware.com/blog/home/post/510

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
   #this one is needed for the project itself
   #SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/milotti_mts/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# Include all the transitive dependencies determined by Spack.
# If we're not running with Spack, this does nothing...
include_directories($ENV{SPACK_TRANSITIVE_INCLUDE_PATH})

endif()
###############  end spack stuff  ##############
# Activate concurrency?
option(ACTIVATE_CONCURRENCY
       "Enable concurrency"
       ON)
       

option(ENABLE_RELEASE_DEBUG_OUTPUT "Enables additional output in release mode" ON)
if(ENABLE_RELEASE_DEBUG_OUTPUT)
add_definitions(-DENABLE_RELEASE_DEBUG_OUTPUT)
endif()
       
# c++ 11 for newer cmake
set(CMAKE_CXX_STANDARD 11)

set (ADDITIONAL_INCLUDE_DIRS CACHE PATH "Additional dirs to search for included files")
set (ADDITIONAL_LIBRARY_DIRS CACHE PATH "Additional dirs to search for linked libraries")
include_directories(${ADDITIONAL_INCLUDE_DIRS})
link_directories(${ADDITIONAL_LIBRARY_DIRS})

#always use new c++11 the old way
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

### it is NOT elegant, but we do not know how to handle gmp on mac, 
### we set the compiler flags manually there
### check if apple
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(APPLE TRUE)
    message("Your are using an apple PC.")
endif()
### compiler
message("CMAKE_CXX_COMPILER:        ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_COMPILER_ID:     ${CMAKE_CXX_COMPILER_ID}")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("using GNU build system")

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -frounding-math")
  if(${APPLE})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG -I/opt/local/include -I/usr/local/include -L/opt/local/lib -lgmpxx -lmpfr -lgmp -lmpfi ")
  endif()

endif()
### Boost
find_package(Boost REQUIRED COMPONENTS "thread")
set(myLibs "${myLibs};${Boost_THREAD_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "random")
set(myLibs "${myLibs};${Boost_RANDOM_LIBRARY}")
find_package(Boost REQUIRED COMPONENTS "system")
set(myLibs "${myLibs};${Boost_SYSTEM_LIBRARY}")

### OPENMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
  message(STATUS "OpenMP_CXX_FLAGS        = ${OpenMP_CXX_FLAGS}")
  message(STATUS "OpenMP_C_FLAGS          = ${OpenMP_C_FLAGS}")
  message(STATUS "OpenMP_EXE_LINKER_FLAGS = ${OpenMP_EXE_LINKER_FLAGS}")
  set (CMAKE_C_FLAGS          "${OpenMP_C_FLAGS} ${CMAKE_C_FLAGS}")
  set (CMAKE_CXX_FLAGS        "${OpenMP_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  #set (CMAKE_C_FLAGS  "${OpenMP_C_FLAGS}   -DOPENMP_ENABLED=1 ${CMAKE_C_FLAGS}  ")
  #set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DOPENMP_ENABLED=1 ${CMAKE_CXX_FLAGS}")
  #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  message(FATAL_ERROR "OpenMP stuff is not found, but labeled as found")
endif()
### NOTE THIS IS A HACK AND SHOULD BE DONE PROPERLY IN FUTURE!!!
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
### last one is for ubuntu 16
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/localdisk/thierry/local/lib/CGAL;/localdisk/thierry/local/lib64/CGAL;/usr/lib/x86_64-linux-gnu/cmake/CGAL")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

### CGAL and related

message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};/localdisk/thierry/local/lib64/CGAL;")
message("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
  message("Found CGAL_LIBRARY: ${CGAL_LIBRARY}")
  message("Found CGAL_INCLUDE_DIRS: ${CGAL_INCLUDE_DIRS}")
  include_directories(${CGAL_INCLUDE_DIRS})
  set(myLibs "${myLibs};${CGAL_LIBRARY}" )
endif()

find_package(GMP REQUIRED)
if(GMP_FOUND)
  message("Found GMP_LIBRARIES: ${GMP_LIBRARIES}")
  message("Found GMP_INCLUDE_DIR: ${GMP_INCLUDE_DIR}")
  include_directories(${GMP_INCLUDE_DIR})
  set(myLibs "${myLibs};${GMP_LIBRARIES}" )
endif()

find_package(Boost REQUIRED COMPONENTS "system")
message("Boost_SYSTEM_LIBRARY: ${Boost_SYSTEM_LIBRARY}")
set(myLibs "${myLibs};${Boost_SYSTEM_LIBRARY}")


find_package(MPFR REQUIRED)
if(MPFR_FOUND)
  message("Found MPFR_LIBRARIES: ${MPFR_LIBRARIES}")
  message("Found MPFR_INCLUDE_DIR: ${MPFR_INCLUDE_DIR}")
  include_directories(${MPFR_INCLUDE_DIR})
  set(myLibs "${myLibs};${MPFR_LIBRARIES}" )
endif()

if( ACTIVATE_CONCURRENCY )
  find_package( TBB REQUIRED )
  message("TBB_LIBRARIES: ${TBB_LIBRARIES}") # variable used by CGAL
  message("TBB_LIBRARY_RELEASE: ${TBB_RELEASE_LIBRARY}") 
  message("TBB_MALLOC_RELEASE_LIBRARY: ${TBB_MALLOC_RELEASE_LIBRARY}")
  if( TBB_FOUND )
    message("tbb found!")
    message("TBB_USE_FILE: ${TBB_USE_FILE}")
    include(${TBB_USE_FILE})
    set(myLibs "${myLibs};${TBB_RELEASE_LIBRARY};${TBB_MALLOC_RELEASE_LIBRARY}")
  endif()
endif()

### cmake debug output CMAKE_EXE_LINKER_FLAGS
message("set content of mylib: ${myLibs}")
message("set comiler flags: ${CMAKE_CXX_FLAGS}")
message("CMAKE_EXE_LINKER_FLAGS: ${CMAKE_EXE_LINKER_FLAGS}")

message("PROJECT_BINARY_DIR ${PROJECT_BINARY_DIR}")
message("PROJECT_BINARY_DIR ${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include/vbl")

#set(myLibs "${myLibs};/localdisk/thierry/tc_install/lib/libkrebs_d.so")
add_subdirectory(src)
add_subdirectory(tools)
install(DIRECTORY "include/vbl" DESTINATION include FILES_MATCHING PATTERN "*.h")
install(FILES "${PROJECT_SOURCE_DIR}/include/vbl.h" DESTINATION include)
