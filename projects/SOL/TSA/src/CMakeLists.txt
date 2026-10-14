# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Time Solver Algebra Source

# Package version.
set(${TSA}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_TSA_Time_Solver_Algorithm.hpp cl_TSA_Monolithic_Time_Solver.hpp
            cl_TSA_Staggered_Time_Solver.hpp cl_TSA_Time_Solver_Factory.hpp cl_TSA_Time_Solver.hpp
)

# Time Solver Library

# Library sources.
set(LIB_SOURCES
    cl_TSA_Time_Solver_Algorithm.cpp cl_TSA_Monolithic_Time_Solver.cpp
    cl_TSA_Staggered_Time_Solver.cpp cl_TSA_Time_Solver_Factory.cpp cl_TSA_Time_Solver.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${SOL_CORE}-lib ${DLA}-lib ${NLA}-lib ${MSI}-lib ${ENM}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${TSA}-lib TSA_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${TSA_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${TSA}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${TSA_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${TSA_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${TSA})

foreach(TPL IN LISTS TSA_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${TSA_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${TSA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${TSA_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    TSA
    TARGET
    ${TSA_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${TSA_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${TSA}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${TSA}-lib DESTINATION ${${TSA}_LIB_INSTALL_DIR} EXPORT ${TSA}Targets)

# Export library into *Targets file
install(EXPORT ${TSA}Targets DESTINATION ${${TSA}_CONFIG_INSTALL_DIR})
