# Install script for directory: /home/schmidt/codes/moris/projects/FEM/MSI/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/MSI" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Adof.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Dof_Manager.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Dof_Type_Enums.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Model_Solver_Interface.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Equation_Object.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Node_Obj.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_Pdof_Host.hpp"
    "/home/schmidt/codes/moris/projects/FEM/MSI/src/cl_MSI_Solver_Interface.hpp"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/MSI" TYPE STATIC_LIBRARY FILES "/home/schmidt/codes/moris/build/projects/FEM/MSI/src/lib/libMSI.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/MSI/MSITargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/MSI/MSITargets.cmake"
         "/home/schmidt/codes/moris/build/projects/FEM/MSI/src/CMakeFiles/Export/share/MSI/MSITargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/MSI/MSITargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/MSI/MSITargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/MSI" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/FEM/MSI/src/CMakeFiles/Export/share/MSI/MSITargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/MSI" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/FEM/MSI/src/CMakeFiles/Export/share/MSI/MSITargets-noconfig.cmake")
  endif()
endif()

