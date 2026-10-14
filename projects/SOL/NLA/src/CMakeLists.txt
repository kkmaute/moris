# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# NonLinear Algebra Source

# Package version.
set(${NLA}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    #cl_NLA_Arc_Length.hpp
    cl_NLA_Convergence.hpp
    cl_NLA_Solver_Relaxation.hpp
    cl_NLA_Solver_Pseudo_Time_Control.hpp
    cl_NLA_Solver_Load_Control.hpp
    cl_NLA_Solver_Nonconformal_Remapping.hpp
    cl_NLA_Newton_Solver.hpp
    cl_NLA_NLBGS.hpp
    cl_NLA_Nonlinear_Algorithm.hpp
    cl_NLA_Nonlinear_Problem.hpp
    cl_NLA_Nonlinear_Solver_Factory.hpp
    cl_NLA_Nonlinear_Solver.hpp
)

# NonLinear Library

# Library sources.
set(LIB_SOURCES
    #cl_NLA_Arc_Length.cpp
    cl_NLA_Newton_Solver.cpp
    cl_NLA_Nonlinear_Solver_Factory.cpp
    cl_NLA_Nonlinear_Algorithm.cpp
    cl_NLA_Nonlinear_Solver.cpp
    cl_NLA_Nonlinear_Problem.cpp
    cl_NLA_Convergence.cpp
    cl_NLA_Solver_Relaxation.cpp
    cl_NLA_Solver_Pseudo_Time_Control.cpp
    cl_NLA_Solver_Load_Control.cpp
    cl_NLA_Solver_Nonconformal_Remapping.cpp
    cl_NLA_NLBGS.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${DLA}-lib ${COM}-lib ${MSI}-lib ${SOL_CORE}-lib ${ENM}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${NLA}-lib NLA_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${NLA_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${NLA}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${NLA_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${NLA_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${NLA})

foreach(TPL IN LISTS NLA_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${NLA_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${NLA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${NLA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    NLA
    TARGET
    ${NLA_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${NLA_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${NLA}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${NLA}-lib DESTINATION ${${NLA}_LIB_INSTALL_DIR} EXPORT ${NLA}Targets)

# Export library into *Targets file
install(EXPORT ${NLA}Targets DESTINATION ${${NLA}_CONFIG_INSTALL_DIR})
