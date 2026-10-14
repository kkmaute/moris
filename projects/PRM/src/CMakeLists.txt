# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Parameters Source

# Package version.
set(${PRM}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    fn_PRM_FEM_Parameters.hpp
    fn_PRM_GEN_Parameters.hpp
    fn_PRM_HMR_Parameters.hpp
    fn_PRM_MSI_Parameters.hpp
    fn_PRM_OPT_Parameters.hpp
    fn_PRM_SOL_Parameters.hpp
    fn_PRM_VIS_Parameters.hpp
    fn_PRM_XTK_Parameters.hpp
    fn_PRM_MORIS_GENERAL_Parameters.hpp
    fn_PRM_MIG_Parameters.hpp
    fn_PRM_WRK_Parameters.hpp
    fn_PRM_STK_Parameters.hpp
)

# Parameters Library

# Library sources.
set(LIB_SOURCES fn_PRM_FEM_Parameters.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${ENM}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${PRM}-lib PRM_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${PRM_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${PRM}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${PRM_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${PRM_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${PRM})

# Third-party dependencies.
foreach(TPL IN LISTS PRM_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${PRM_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${PRM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${PRM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    PRM
    TARGET
    ${PRM_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${PRM_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${HEADER_INSTALL_DIR}/${PRM})

# Install and export targets.
install(
    TARGETS ${PRM}-lib
    EXPORT ${PRM}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${PRM}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${PRM}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${PRM}
)

# Install exported target definitions.
install(EXPORT ${PRM}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${PRM})
