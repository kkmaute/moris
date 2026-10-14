# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Exceptions Source

# Public headers.
set(HEADERS exceptions.hpp cl_invalid_size.hpp cl_out_of_range.hpp cl_tic_toc_error.hpp)

# Exceptions Library

add_library(${EXC}-lib INTERFACE)
target_include_directories(
    ${EXC}-lib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${${EXC}_HEADER_INSTALL_DIR}>
)

# Exceptions Executable

# Generate and install the common component package configuration.
moris_install_component_package(EXC TARGET ${EXC}-lib DEPENDENCIES TPLS ${EXC_TPL_DEPENDENCIES})

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${EXC}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${EXC}-lib DESTINATION ${${EXC}_LIB_INSTALL_DIR} EXPORT ${EXC}Targets)

# Install exported target definitions.
install(EXPORT ${EXC}Targets DESTINATION ${${EXC}_CONFIG_INSTALL_DIR})
