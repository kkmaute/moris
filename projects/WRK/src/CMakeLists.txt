# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Parameters Source

# Package version.
set(${WRK}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_WRK_Performer_Manager.hpp
    cl_WRK_GEN_Performer.hpp
    cl_WRK_Performer.hpp
    cl_WRK_Workflow.hpp
    cl_WRK_Workflow_Factory.hpp
    cl_WRK_Workflow_HMR_XTK.hpp
    cl_WRK_Workflow_STK_XTK.hpp
    cl_WRK_Workflow_STK_FEM.hpp
    cl_WRK_perform_refinement.hpp
    cl_WRK_perform_remeshing.hpp
    cl_WRK_DataBase_Performer.hpp
    cl_WRK_Reinitialize_Performer.hpp
)

set(GEN_INCLUDES "ADV" "ADV/fields" "PDV")

# Include dependency directories.
foreach(GEN_INCLUDES IN LISTS GEN_INCLUDES)
    include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
endforeach()

# Parameters Library

# Library sources.
set(LIB_SOURCES
    cl_WRK_Performer_Manager.cpp
    cl_WRK_Performer.cpp
    cl_WRK_GEN_Performer.cpp
    cl_WRK_Workflow_Factory.cpp
    cl_WRK_Workflow_HMR_XTK.cpp
    cl_WRK_Workflow_STK_XTK.cpp
    cl_WRK_Workflow_STK_FEM.cpp
    fn_WRK_Workflow_Main_Interface.cpp
    cl_WRK_perform_refinement.cpp
    cl_WRK_perform_remeshing.cpp
    cl_WRK_DataBase_Performer.cpp
    cl_WRK_Reinitialize_Performer.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${INT}-lib
    ${XTK}-lib
    ${HMR}-lib
    ${GEN}-lib
    ${MDL}-lib
    ${OPT}-lib
    ${MIG}-lib
)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${WRK}-lib WRK_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${WRK_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${WRK}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${WRK_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${WRK_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${WRK})

# Third-party dependencies.
foreach(TPL IN LISTS WRK_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${WRK_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${WRK_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${WRK_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    WRK
    TARGET
    ${WRK_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${WRK_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${HEADER_INSTALL_DIR}/${WRK})

# Install and export targets.
install(
    TARGETS ${WRK}-lib
    EXPORT ${WRK}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${WRK}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${WRK}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${WRK}
)

# Install exported target definitions.
install(EXPORT ${WRK}Targets DESTINATION ${CONFIG_INSTALL_DIR}/${WRK})
