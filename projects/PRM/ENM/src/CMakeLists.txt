# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Enums Source

# Package version.
set(${ENM}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_FEM_Enums.hpp
    GEN_Data_Types.hpp
    cl_MSI_Dof_Type_Enums.hpp
    cl_MTK_Enums.hpp
    cl_NLA_Nonlinear_Solver_Enums.hpp
    cl_SOL_Enums.hpp
    cl_TSA_Time_Solver_Enums.hpp
    cl_VIS_Output_Enums.hpp
    cl_XTK_Enums.hpp
    fn_enum_macros.hpp
    fn_enum_macros_generated.hpp
    enums.hpp
)

# Enums Library

# Library sources.
set(LIB_SOURCES GEN_Data_Types.cpp cl_MTK_Enums.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${ASR}-lib ${CNT}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Enums Target

# Define the library and its public build requirements.
moris_add_module(${ENM}-lib ENM_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${ENM_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${ENM}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${ENM_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${ENM_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${ENM})

foreach(TPL IN LISTS MAP_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${ENM_BUILD_TARGET} PUBLIC ${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${ENM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${ENM_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    ENM
    TARGET
    ${ENM_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${ENM_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${ENM}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${ENM}-lib DESTINATION ${${ENM}_LIB_INSTALL_DIR} EXPORT ${ENM}Targets)

# Install exported target definitions.
install(EXPORT ${ENM}Targets DESTINATION ${${ENM}_CONFIG_INSTALL_DIR})
