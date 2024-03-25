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
    cl_Vector.hpp
    cl_Dist_Map.hpp
    cl_Map.hpp
    cl_Mini_Maps.hpp
	cl_Validator.hpp
	cl_Parameter.hpp
	cl_Parameter_List.hpp
    cl_Tuple.hpp
    containers.hpp
    fn_zip.hpp
    cl_Fortran.hpp )
    
# List library source files
set(LIB_SOURCES
	cl_Vector.cpp
	cl_Validator.cpp
	cl_Parameter.cpp
	cl_Parameter_List.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Library

#add_library(${CNT}-lib INTERFACE)
add_library(${CNT}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${CNT}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${CNT}_HEADER_INSTALL_DIR}> )
target_link_libraries(${CNT}-lib PUBLIC ${LIB_DEPENDENCIES})

foreach(TPL ${CNT_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${CNT}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${CNT}-lib 
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${CNT}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${CNT_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND CNT_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CNT}Config.cmake.in
    ${${CNT}_CONFIG_BUILD_DIR}/${CNT}Config.cmake
    INSTALL_DESTINATION 
    ${${CNT}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${CNT}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${CNT}_CONFIG_BUILD_DIR}/${CNT}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Installation

# Install *Config and *ConfigVersion files
install(FILES ${${CNT}_CONFIG_BUILD_DIR}/${CNT}Config.cmake 
    ${${CNT}_CONFIG_BUILD_DIR}/${CNT}ConfigVersion.cmake
    DESTINATION ${${CNT}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${CNT}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${CNT}-lib
    DESTINATION ${${CNT}_LIB_INSTALL_DIR}
    EXPORT ${CNT}Targets )

# Export library in to *Targets file
install(EXPORT ${CNT}Targets
    DESTINATION ${${CNT}_CONFIG_INSTALL_DIR} )

