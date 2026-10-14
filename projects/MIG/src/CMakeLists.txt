# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# MIG Source

# Package version.
set(${MIG}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_MIG_Periodic_2D.hpp cl_MIG_Periodic_3D.hpp cl_MIG_Mesh_Editor.hpp cl_MIG.cpp)

include_directories(${MORIS_PACKAGE_DIR}/include)

set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

# Include dependency directories.
foreach(GEN_INCLUDES IN LISTS GEN_INCLUDES)
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# MIG Library

# Library sources.
set(LIB_SOURCES cl_MIG_Periodic_2D.cpp cl_MIG_Periodic_3D.cpp cl_MIG_Mesh_Editor.cpp cl_MIG.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${GEN}-lib ${COM}-lib ${LINALG}-lib ${HMR}-lib ${TOL}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${MIG}-lib MIG_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${MIG_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${MIG}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${MIG_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MIG_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${MIG})

foreach(TPL IN LISTS MIG_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${MIG_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${MIG_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${MIG_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    MIG
    TARGET
    ${MIG_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${MIG_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${MIG}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${MIG}-lib DESTINATION ${${MIG}_LIB_INSTALL_DIR} EXPORT ${MIG}Targets)

# Export library into *Targets file
install(EXPORT ${MIG}Targets DESTINATION ${${MIG}_CONFIG_INSTALL_DIR})
