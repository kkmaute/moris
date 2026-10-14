# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Algorithms Source

# Package version.
set(${ALG}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    algorithms.hpp
    fn_equal_to.hpp
    fn_equal_to.tpp
    fn_check_equal.hpp
    fn_hash_value.hpp
    fn_num_digits.hpp
    fn_clip_value.hpp
)

# Algorithms Library

# Sources.
set(LIB_SOURCES fn_equal_to.cpp)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library and set properties
moris_add_module(${ALG}-lib ALG_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${ALG_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${ALG}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${ALG_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${ALG_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${ALG})

# Third-party dependencies.
foreach(TPL IN LISTS ALG_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${ALG_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    # Later versions of CMake (3.11+) allow include directories and
    # compile definitions to be included in an imported library. If this
    # system ever gets updated to a newer version of CMake, shifting both
    # to the imported library would be good practice. Imported library
    # targets were added to add inter-tpl dependencies where needed
    # (PETSc needs SuperLU, for example).
    #     -- Matthew Ryan, June 18, 2019
    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${ALG_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${ALG_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# Generate and install the common component package configuration.
moris_install_component_package(
    ALG
    TARGET
    ${ALG_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${ALG_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${ALG}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${ALG}-lib DESTINATION ${${ALG}_LIB_INSTALL_DIR} EXPORT ${ALG}Targets)

# Export library into *Targets file
install(EXPORT ${ALG}Targets DESTINATION ${${ALG}_CONFIG_INSTALL_DIR})
