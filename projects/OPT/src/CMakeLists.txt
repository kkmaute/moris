# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Optimization Source

# Package version.
set(${OPT}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_OPT_Algorithm.hpp
    cl_OPT_Algorithm_GCMMA.hpp
    cl_OPT_Algorithm_LBFGS.hpp
    cl_OPT_Algorithm_ADAM.hpp
    cl_OPT_Algorithm_SQP.hpp
    cl_OPT_Algorithm_Sweep.hpp
    cl_OPT_Manager.hpp
    cl_OPT_Problem.hpp
    cl_OPT_Criteria_Interface.hpp
    cl_OPT_Interface_Manager.hpp
    cl_OPT_Interface_User_Defined.hpp
    fn_OPT_create_problem.hpp
    fn_OPT_create_interface.hpp
    fn_OPT_create_algorithm.hpp
)

# Optimization Library

# Sources.
set(LIB_SOURCES
    cl_OPT_Algorithm.cpp
    cl_OPT_Algorithm_GCMMA.cpp
    cl_OPT_Algorithm_LBFGS.cpp
    cl_OPT_Algorithm_ADAM.cpp
    cl_OPT_Algorithm_SQP.cpp
    cl_OPT_Algorithm_Sweep.cpp
    cl_OPT_Manager.cpp
    cl_OPT_Problem.cpp
    cl_OPT_Criteria_Interface.cpp
    cl_OPT_Interface_Manager.cpp
    cl_OPT_Problem_User_Defined.cpp
    cl_OPT_Interface_User_Defined.cpp
    fn_OPT_create_problem.cpp
    fn_OPT_create_interface.cpp
    fn_OPT_create_algorithm.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${CNT}-lib ${LINALG}-lib ${COM}-lib ${MORIS_BASE_LIBS})
#message(INFO "base lib = ${MORIS_BASE_LIBS}\n\n\n\n\n")
# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${OPT}-lib OPT_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${OPT_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${OPT}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${OPT_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${OPT_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${OPT})
target_compile_definitions(${OPT_BUILD_TARGET} PUBLIC "-DF77ADD_")

# Third-party dependencies.
foreach(TPL IN LISTS OPT_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${OPT_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${OPT_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${OPT_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    OPT
    TARGET
    ${OPT_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${OPT_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${OPT}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${OPT}-lib DESTINATION ${${OPT}_LIB_INSTALL_DIR} EXPORT ${OPT}Targets)

# Export library into *Targets file
install(EXPORT ${OPT}Targets DESTINATION ${${OPT}_CONFIG_INSTALL_DIR})
