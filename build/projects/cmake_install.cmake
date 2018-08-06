# Install script for directory: /home/schmidt/codes/moris/projects

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

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/schmidt/codes/moris/build/projects/ALG/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/COM/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/DLA/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/GEN/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/IOS/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/MOD/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/FEM/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/TOL/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/MRS/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/LNA/cmake_install.cmake")
  include("/home/schmidt/codes/moris/build/projects/mains/cmake_install.cmake")

endif()

