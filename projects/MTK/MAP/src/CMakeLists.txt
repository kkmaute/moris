# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# MAP Source

# Package version.
set(${MAP}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS cl_MTK_Mapper_Node.hpp cl_MTK_Mapper.hpp)

# MAP Library

# Library sources.
set(LIB_SOURCES cl_MTK_Mapper.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${CNT}-lib
    ${MTK}-lib
    ${NLA}-lib
    ${INT}-lib
    ${TSA}-lib
    ${MDL}-lib
)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${MAP}-lib MAP_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${MAP_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${MAP}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${MAP_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MAP_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${MAP})

foreach(TPL IN LISTS MAP_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${MAP_BUILD_TARGET} PUBLIC ${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${MAP_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${MAP_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    MAP
    TARGET
    ${MAP_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${MAP_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${MAP}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${MAP}-lib DESTINATION ${${MAP}_LIB_INSTALL_DIR} EXPORT ${MAP}Targets)

# Export library into *Targets file
install(EXPORT ${MAP}Targets DESTINATION ${${MAP}_CONFIG_INSTALL_DIR})
