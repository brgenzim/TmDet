#/*
# *  CMakeLists.txt
# *  Thu Nov 04 10:39:23 2023
# *  Copyright (c)  2001-2024 Gabor E. Tusnady
# *  <tusi@enzim.hu>
# *  All rights reserved.
# */


# Project related settings

  PROJECT(TmDet)

  SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
  SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
  SET( ${PROJECT_NAME}_VERSION_MICRO 0 )
  SET( ${PROJECT_NAME}_VERSION_BUILD 0 )

# CMAKE settings

  CMAKE_MINIMUM_REQUIRED( VERSION 3.22 )
  MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)

# Global settings

  SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

# Build types

  SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "Choose the type of build, options are: ${CMAKE_CONFIGURATION_TYPES}." FORCE)
  IF(NOT CMAKE_BUILD_TYPE)
  ENDIF(NOT CMAKE_BUILD_TYPE)

  SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Path settings

  SET( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )
#  SET( CMAKE_INSTALL_PREFIX /usr/local )
#  SET( CMAKE_INSTALL_PREFIX_LIB /usr/local/lib )
#  SET( CMAKE_INSTALL_PREFIX_BIN /usr/local/bin )
#  SET( CMAKE_INSTALL_PREFIX_INCLUDE /usr/local/include )
  SET( CMAKE_INSTALL_PREFIX . )
  SET( CMAKE_INSTALL_PREFIX_LIB ${PROJECT_SOURCE_DIR}/lib )
  SET( CMAKE_INSTALL_PREFIX_BIN ${PROJECT_SOURCE_DIR}/bin )
  SET( CMAKE_INSTALL_PREFIX_INCLUDE ${PROJECT_SOURCE_DIR}/include )


  MESSAGE( "-- Binaries:  " ${CMAKE_INSTALL_PREFIX_BIN} )
  MESSAGE( "-- Libraries: " ${CMAKE_INSTALL_PREFIX_LIB} )
  MESSAGE( "-- Includes:  " ${CMAKE_INSTALL_PREFIX_INCLUDE} )

# Debug definitions
  IF( CMAKE_BUILD_TYPE STREQUAL "DEBUG" OR
      CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
   ADD_DEFINITIONS( "-DDEBUG" )
  ENDIF()

# json support
#find_package(nlohmann_json 3.2.0 REQUIRED)

# gemmi support
find_package(gemmi 0.6.4 CONFIG REQUIRED)

# Include files (system)
  INCLUDE_DIRECTORIES( AFTER SYSTEM 
    "${PROJECT_SOURCE_DIR}/src/Utils"
    "${PROJECT_SOURCE_DIR}/src/TmdetLib"
    "${PROJECT_SOURCE_DIR}/src/Types"
    "${PROJECT_SOURCE_DIR}/../contrib/pugixml/src/"
  )

link_directories(${CMAKE_INSTALL_PREFIX_LIB})

# Compiler settings
  ADD_DEFINITIONS( "-Wall -pedantic -ggdb3 -DDEBUG -Wno-variadic-macros -std=c++20")

# Subdirectories
  ADD_SUBDIRECTORY ( src )
