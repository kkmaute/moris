# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Sol Core Source

# Package version.
set(${SOL_CORE}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_SOL_Amesos_GenOp.hpp
    cl_SOL_Warehouse.hpp
    cl_SOL_Matrix_Vector_Factory.hpp
    cl_SOL_Dist_Vector.hpp
    cl_SOL_Dist_Matrix.hpp
    cl_SOL_Dist_Map.hpp
    cl_Communicator_Epetra.hpp
    cl_Map_Epetra.hpp
    cl_Sparse_Matrix_EpetraFECrs.hpp
    cl_Vector_Epetra.hpp
)

if(${MORIS_HAVE_PETSC})
    list(APPEND HEADERS "cl_MatrixPETSc.hpp" "cl_Vector_PETSc_Multi.hpp"
         "cl_SOL_Dist_Map_Custom.hpp"
    )
endif()

# Sol Core Library

# Library sources.
set(LIB_SOURCES
    cl_SOL_Amesos_GenOp.cpp
    cl_SOL_Warehouse.cpp
    cl_SOL_Matrix_Vector_Factory.cpp
    cl_SOL_Dist_Vector.cpp
    cl_Communicator_Epetra.cpp
    cl_Map_Epetra.cpp
    cl_Sparse_Matrix_EpetraFECrs.cpp
    cl_Vector_Epetra.cpp
)

if(${MORIS_HAVE_PETSC})
    list(APPEND LIB_SOURCES "cl_MatrixPETSc.cpp" "cl_Vector_PETSc_Multi.cpp"
         "cl_SOL_Dist_Map_Custom.cpp"
    )
endif()

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${DLA}-lib ${ENM}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${SOL_CORE}-lib SOL_CORE_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${SOL_CORE_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                                    $<INSTALL_INTERFACE:${${SOL_CORE}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${SOL_CORE_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${SOL_CORE_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${SOL_CORE})

foreach(TPL IN LISTS SOL_CORE_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${SOL_CORE_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${SOL_CORE_BUILD_TARGET} INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS}
                                           ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${SOL_CORE_BUILD_TARGET} INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    SOL_CORE
    TARGET
    ${SOL_CORE_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${SOL_CORE_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${SOL_CORE}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${SOL_CORE}-lib DESTINATION ${${SOL_CORE}_LIB_INSTALL_DIR}
        EXPORT ${SOL_CORE}Targets
)

# Export library into *Targets file
install(EXPORT ${SOL_CORE}Targets DESTINATION ${${SOL_CORE}_CONFIG_INSTALL_DIR})
