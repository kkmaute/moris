# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Visualization

# Package version.
set(${VIS}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_VIS_Output_Manager.hpp cl_VIS_Visualization_Mesh.hpp cl_VIS_Vertex_Visualization.hpp
    cl_VIS_Cell_Cluster_Visualization.hpp cl_VIS_Factory.hpp cl_VIS_Side_Cluster_Visualization.hpp
)

# Visualization Library

# Library sources.
set(LIB_SOURCES
    cl_VIS_Output_Manager.cpp cl_VIS_Visualization_Mesh.cpp cl_VIS_Cell_Cluster_Visualization.cpp
    cl_VIS_Factory.cpp cl_VIS_Side_Cluster_Visualization.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${COM}-lib ${CNT}-lib ${HMR}-lib ${MTK}-lib ${ENM}-lib
                     #${MDL}-lib
)

set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

# Include dependency directories.
foreach(GEN_INCLUDES IN LISTS GEN_INCLUDES)
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${VIS}-lib VIS_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${VIS_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${VIS}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${VIS_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${VIS_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${VIS})

foreach(TPL IN LISTS VIS_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${VIS_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${VIS_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${VIS_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    VIS
    TARGET
    ${VIS_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${VIS_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${VIS}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${VIS}-lib DESTINATION ${${VIS}_LIB_INSTALL_DIR} EXPORT ${VIS}Targets)

# Export library into *Targets file
install(EXPORT ${VIS}Targets DESTINATION ${${VIS}_CONFIG_INSTALL_DIR})
