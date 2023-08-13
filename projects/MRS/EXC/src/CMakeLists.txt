#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Exceptions Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    exceptions.hpp
    cl_invalid_size.hpp
    cl_out_of_range.hpp
    cl_tic_toc_error.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Library

add_library(${EXC}-lib INTERFACE)
target_include_directories(${EXC}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${EXC}_HEADER_INSTALL_DIR}> )

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)
	
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${EXC_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND EXC_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${EXC}Config.cmake.in
    ${${EXC}_CONFIG_BUILD_DIR}/${EXC}Config.cmake
    INSTALL_DESTINATION 
    ${${EXC}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${EXC}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${EXC}_CONFIG_BUILD_DIR}/${EXC}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Installation

# Install *Config and *ConfigVersion files
install(FILES ${${EXC}_CONFIG_BUILD_DIR}/${EXC}Config.cmake 
    ${${EXC}_CONFIG_BUILD_DIR}/${EXC}ConfigVersion.cmake
    DESTINATION ${${EXC}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${EXC}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${EXC}-lib
    DESTINATION ${${EXC}_LIB_INSTALL_DIR}
    EXPORT ${EXC}Targets )

# Export library in to *Targets file
install(EXPORT ${EXC}Targets
    DESTINATION ${${EXC}_CONFIG_INSTALL_DIR} )

