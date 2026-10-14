# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Model Solver Interface Source

# Package version.
set(${MSI}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_MSI_Adof.hpp
    cl_MSI_Dof_Manager.hpp
    cl_MSI_Model_Solver_Interface.hpp
    cl_MSI_Equation_Model.hpp
    cl_MSI_Equation_Object.hpp
    cl_MSI_Equation_Set.hpp
    cl_MSI_Pdof_Host.hpp
    cl_MSI_Solver_Interface.hpp
    cl_MSI_Node_Proxy.hpp
    cl_MSI_Element_Proxy.hpp
    cl_MSI_Multigrid.hpp
    cl_MSI_Design_Variable_Interface.hpp
    fn_MSI_get_mesh_index_for_dof_type.hpp
)

# Model Solver Interface Library

# Library sources.
set(LIB_SOURCES
    cl_MSI_Dof_Manager.cpp
    cl_MSI_Pdof_Host.cpp
    cl_MSI_Model_Solver_Interface.cpp
    cl_MSI_Equation_Object.cpp
    cl_MSI_Equation_Set.cpp
    cl_MSI_Equation_Model.cpp
    cl_MSI_Multigrid.cpp
    cl_MSI_Solver_Interface.cpp
    cl_MSI_Design_Variable_Interface.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${COM}-lib ${CNT}-lib ${DLA}-lib ${PRM}-lib ${ENM}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${MSI}-lib MSI_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${MSI_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${MSI}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${MSI_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MSI_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${MSI})

foreach(TPL IN LISTS MSI_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${MSI_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${MSI_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${MSI_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    MSI
    TARGET
    ${MSI_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${MSI_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${MSI}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${MSI}-lib DESTINATION ${${MSI}_LIB_INSTALL_DIR} EXPORT ${MSI}Targets)

# Export library into *Targets file
install(EXPORT ${MSI}Targets DESTINATION ${${MSI}_CONFIG_INSTALL_DIR})
