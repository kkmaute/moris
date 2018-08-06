# Install script for directory: /home/schmidt/codes/moris/projects/MRS/CON/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/MRS/CON" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Array.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Bi_Map.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Bitbool.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Bitset.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_BoostBitset.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Cell.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Dist_Map.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Map.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Param_List.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/cl_Tuple.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/containers.hpp"
    "/home/schmidt/codes/moris/projects/MRS/CON/src/fn_zip.hpp"
    )
endif()

