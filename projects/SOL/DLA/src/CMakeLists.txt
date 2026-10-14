# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Distributed Linear Algebra Source

# Package version.
set(${DLA}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_DLA_Eigen_Solver.hpp
    cl_DLA_Linear_Solver_Amesos.hpp
    cl_DLA_Linear_Solver_Aztec.hpp
    cl_DLA_Linear_Solver_ML.hpp
    cl_DLA_Linear_Solver_Belos.hpp
    cl_DLA_Linear_Solver.hpp
    cl_DLA_Linear_Problem.hpp
    cl_DLA_Linear_System_Trilinos.hpp
    cl_DLA_Solver_Factory.hpp
    cl_DLA_Solver_Interface.hpp
    cl_DLA_Linear_Solver_Algorithm.hpp
    cl_DLA_Linear_Solver_Algorithm_Petsc.hpp
    cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp
    cl_DLA_Preconditioner.hpp
    cl_DLA_Preconditioner_Trilinos.hpp
    cl_DLA_Geometric_Multigrid.hpp
    fn_convert_epetra_operator_to_matrix.hpp
)

if(${MORIS_HAVE_PETSC})
    list(APPEND HEADERS cl_DLA_Linear_System_PETSc.hpp cl_DLA_Linear_Solver_PETSc.hpp
         cl_DLA_Preconditioner_PETSc.hpp
    )
endif()

if(${MORIS_HAVE_SLEPC})
    list(APPEND HEADERS cl_DLA_Eigen_Solver_SLEPc.hpp)
endif()

# distlinalg Library

# Library sources.
set(LIB_SOURCES
    cl_DLA_Eigen_Solver.cpp
    cl_DLA_Linear_Solver_Amesos.cpp
    cl_DLA_Linear_Solver_Aztec.cpp
    cl_DLA_Linear_Solver_ML.cpp
    cl_DLA_Linear_Solver_Belos.cpp
    cl_DLA_Linear_System_Trilinos.cpp
    cl_DLA_Linear_Solver_Algorithm_Trilinos.cpp
    cl_DLA_Linear_Solver.cpp
    cl_DLA_Linear_Problem.cpp
    cl_DLA_Geometric_Multigrid.cpp
    cl_DLA_Solver_Interface.cpp
    cl_DLA_Preconditioner_Trilinos.cpp
    cl_DLA_Solver_Factory.cpp
)

if(${MORIS_HAVE_PETSC})
    list(APPEND LIB_SOURCES cl_DLA_Linear_System_PETSc.cpp cl_DLA_Linear_Solver_PETSc.cpp
         cl_DLA_Preconditioner_PETSc.cpp
    )
endif()

if(${MORIS_HAVE_SLEPC})
    list(APPEND LIB_SOURCES cl_DLA_Eigen_Solver_SLEPc.cpp)
endif()

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${COM}-lib ${MTK}-lib ${SOL_CORE}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${DLA}-lib DLA_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${DLA_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${DLA}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${DLA_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${DLA_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${DLA})

foreach(TPL IN LISTS DLA_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${DLA_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${DLA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${DLA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    DLA
    TARGET
    ${DLA_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${DLA_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${DLA}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${DLA}-lib DESTINATION ${${DLA}_LIB_INSTALL_DIR} EXPORT ${DLA}Targets)

# Export library into *Targets file
install(EXPORT ${DLA}Targets DESTINATION ${${DLA}_CONFIG_INSTALL_DIR})
