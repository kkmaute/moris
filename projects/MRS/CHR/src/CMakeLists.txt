#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Chronos Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    chronos.hpp
    cl_Stopwatch.hpp
    cl_Waitbar.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Library

add_library(${CHR}-lib INTERFACE)
target_include_directories(${CHR}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${CHR}_HEADER_INSTALL_DIR}> )

foreach(TPL ${CHR_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${CHR}-lib
    	INTERFACE ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${CHR}-lib 
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${CHR}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${CHR_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND CHR_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CHR}Config.cmake.in
    ${${CHR}_CONFIG_BUILD_DIR}/${CHR}Config.cmake
    INSTALL_DESTINATION 
    ${${CHR}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${CHR}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${CHR}_CONFIG_BUILD_DIR}/${CHR}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Installation

# Install *Config and *ConfigVersion files
install(FILES ${${CHR}_CONFIG_BUILD_DIR}/${CHR}Config.cmake 
    ${${CHR}_CONFIG_BUILD_DIR}/${CHR}ConfigVersion.cmake
    DESTINATION ${${CHR}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${CHR}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${CHR}-lib
    DESTINATION ${${CHR}_LIB_INSTALL_DIR}
    EXPORT ${CHR}Targets )

# Export library in to *Targets file
install(EXPORT ${CHR}Targets
    DESTINATION ${${CHR}_CONFIG_INSTALL_DIR} )
	
