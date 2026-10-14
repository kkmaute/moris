# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Tools Source

# Package version.
set(${TOL}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_Debug.hpp cl_Enums.hpp cl_Geometry.hpp cl_TOL_Memory_Map.hpp fn_TOL_Capacities.hpp)

# Tools Library

# Sources.
set(LIB_SOURCES cl_Debug.cpp cl_Geometry.cpp cl_TOL_Memory_Map.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${ALG}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${TOL}-lib TOL_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${TOL_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${TOL}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${TOL_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${TOL_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${TOL})

foreach(TPL IN LISTS TOL_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${TOL_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${TOL_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${TOL_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    TOL
    TARGET
    ${TOL_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${TOL_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${TOL}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${TOL}-lib DESTINATION ${${TOL}_LIB_INSTALL_DIR} EXPORT ${TOL}Targets)

# Export library into *Targets file
install(EXPORT ${TOL}Targets DESTINATION ${${TOL}_CONFIG_INSTALL_DIR})
