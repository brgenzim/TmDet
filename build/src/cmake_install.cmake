# Install script for directory: /home/tusi/works/pdbtm_4.0/TmDet/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ".")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "DEBUG")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/lib/libTmDetLib.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/lib" TYPE STATIC_LIBRARY FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/libTmDetLib.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/test")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/test")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/test2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/test2")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/test2")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/fragment_cif")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/fragment_cif")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/dssp")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/dssp")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/domdist")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/tusi/works/pdbtm_4.0/TmDet/bin" TYPE EXECUTABLE FILES "/home/tusi/works/pdbtm_4.0/TmDet/build/src/domdist_af")
  if(EXISTS "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af"
         OLD_RPATH "/home/tusi/works/pdbtm_4.0/TmDet/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/home/tusi/works/pdbtm_4.0/TmDet/bin/domdist_af")
    endif()
  endif()
endif()

