# Install script for directory: /home/schmidt/codes/moris/projects/DLA/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/DLA" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Communicator_Epetra.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_DistLinAlg_Enums.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver_Amesos2.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver_Amesos.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver_Aztec.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver_PETSc.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Linear_Solver_Trilinos.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Map_Class.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Map_Epetra.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Map_PETSc.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_MatrixPETSc.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Matrix_Vector_Factory.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Model_Solver_Interface_Solver.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Solver_Factory.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Solver_Input.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Solver_Input_Test.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Sparse_Matrix_EpetraFECrs.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Sparse_Matrix.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Vector_Epetra.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_Vector.hpp"
    "/home/schmidt/codes/moris/projects/DLA/src/cl_VectorPETSc.hpp"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/DLA" TYPE STATIC_LIBRARY FILES "/home/schmidt/codes/moris/build/projects/DLA/src/lib/libDLA.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/DLA/DLATargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/DLA/DLATargets.cmake"
         "/home/schmidt/codes/moris/build/projects/DLA/src/CMakeFiles/Export/share/DLA/DLATargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/DLA/DLATargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/DLA/DLATargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/DLA" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/DLA/src/CMakeFiles/Export/share/DLA/DLATargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/DLA" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/DLA/src/CMakeFiles/Export/share/DLA/DLATargets-noconfig.cmake")
  endif()
endif()

