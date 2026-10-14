# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

#  Model Source

# Package version.
set(${MDL}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_MDL_Model.hpp
            #cl_MDL_Mesh_Model_Helper.hpp
)

# fem/mdl Library

# Library sources.
set(LIB_SOURCES cl_MDL_Model.cpp
                #cl_MDL_Mesh_Model_Helper.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES
    ${INT}-lib
    ${MTK}-lib
    ${NLA}-lib
    ${TSA}-lib
    ${SDF}-lib
    ${DLA}-lib
    ${PRM}-lib
    ${VIS}-lib
)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${MDL}-lib MDL_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${MDL_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${MDL}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${MDL_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MDL_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${MDL})

foreach(TPL IN LISTS MDL_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${MDL_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${MDL_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${MDL_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    MDL
    TARGET
    ${MDL_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${MDL_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${MDL}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${MDL}-lib DESTINATION ${${MDL}_LIB_INSTALL_DIR} EXPORT ${MDL}Targets)

# Export library into *Targets file
install(EXPORT ${MDL}Targets DESTINATION ${${MDL}_CONFIG_INSTALL_DIR})
