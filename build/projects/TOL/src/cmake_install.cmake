# Install script for directory: /home/schmidt/codes/moris/projects/TOL/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
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
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/TOL" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/TOL/src/tools.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Debug.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Enums.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_FunctionFactory.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Function.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Function_Levelset.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Geometry.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Interpolation.hpp"
    "/home/schmidt/codes/moris/projects/TOL/src/cl_Pairing.hpp"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/TOL" TYPE STATIC_LIBRARY FILES "/home/schmidt/codes/moris/build/projects/TOL/src/lib/libTOL.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/TOL/TOLTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/TOL/TOLTargets.cmake"
         "/home/schmidt/codes/moris/build/projects/TOL/src/CMakeFiles/Export/share/TOL/TOLTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/TOL/TOLTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/TOL/TOLTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/TOL" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/TOL/src/CMakeFiles/Export/share/TOL/TOLTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/TOL" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/TOL/src/CMakeFiles/Export/share/TOL/TOLTargets-noconfig.cmake")
  endif()
endif()

