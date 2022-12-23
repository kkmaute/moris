#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Containers Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    cl_Array.hpp
    cl_Bi_Map.hpp
    cl_Bitbool.hpp
    cl_Bitset.hpp
    cl_BoostBitset.hpp
    cl_Cell.hpp
    cl_Dist_Map.hpp
    cl_Map.hpp
    cl_Param_List.hpp
    cl_Tuple.hpp
    containers.hpp
    fn_zip.hpp
    cl_Fortran.hpp )
    
# List library source files
set(LIB_SOURCES
	cl_Cell.cpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Library

#add_library(${CON}-lib INTERFACE)
add_library(${CON}-lib STATIC ${LIB_SOURCES})
target_include_directories(${CON}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${CON}_HEADER_INSTALL_DIR}> )

foreach(TPL ${CON_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${CON}-lib
    	INTERFACE ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${CON}-lib 
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${CON}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${CON_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND CON_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CON}Config.cmake.in
    ${${CON}_CONFIG_BUILD_DIR}/${CON}Config.cmake
    INSTALL_DESTINATION 
    ${${CON}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${CON}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${CON}_CONFIG_BUILD_DIR}/${CON}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Installation

# Install *Config and *ConfigVersion files
install(FILES ${${CON}_CONFIG_BUILD_DIR}/${CON}Config.cmake 
    ${${CON}_CONFIG_BUILD_DIR}/${CON}ConfigVersion.cmake
    DESTINATION ${${CON}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${CON}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${CON}-lib
    DESTINATION ${${CON}_LIB_INSTALL_DIR}
    EXPORT ${CON}Targets )

# Export library in to *Targets file
install(EXPORT ${CON}Targets
    DESTINATION ${${CON}_CONFIG_INSTALL_DIR} )

