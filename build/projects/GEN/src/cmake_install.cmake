# Install script for directory: /home/schmidt/codes/moris/projects/GEN/src

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/GEN" TYPE FILE FILES
    "/home/schmidt/codes/moris/projects/GEN/src/geomeng.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/cl_GeometryEngineFactory.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/cl_GeometryEngine.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/cl_GeometryEngine_Levelset.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Core.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Data.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Gen.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Mesh_Data.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Triangle_File.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/cl_Ge_SDF_Triangle.hpp"
    "/home/schmidt/codes/moris/projects/GEN/src/SDF/GeUtilities.hpp"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/GEN" TYPE STATIC_LIBRARY FILES "/home/schmidt/codes/moris/build/projects/GEN/src/lib/libGEN.a")
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/GEN/GENTargets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/GEN/GENTargets.cmake"
         "/home/schmidt/codes/moris/build/projects/GEN/src/CMakeFiles/Export/share/GEN/GENTargets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/GEN/GENTargets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/share/GEN/GENTargets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/GEN" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/GEN/src/CMakeFiles/Export/share/GEN/GENTargets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/share/GEN" TYPE FILE FILES "/home/schmidt/codes/moris/build/projects/GEN/src/CMakeFiles/Export/share/GEN/GENTargets-noconfig.cmake")
  endif()
endif()

