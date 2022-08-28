#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Tools Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${TOL}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Debug.hpp
    cl_Enums.hpp
    cl_Geometry.hpp
    cl_TOL_Memory_Map.hpp
    fn_TOL_Capacities.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Library

# List source files
set(LIB_SOURCES
    cl_Debug.cpp
    cl_Geometry.cpp
    cl_TOL_Memory_Map.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${ALG}-lib
	)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${TOL}-lib STATIC ${LIB_SOURCES})
target_include_directories(${TOL}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${TOL}_HEADER_INSTALL_DIR}> )
target_link_libraries(${TOL}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${TOL}-lib PROPERTIES OUTPUT_NAME ${TOL})

foreach(TPL ${TOL_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${TOL}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${TOL}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${TOL}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${TOL_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND TOL_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND TOL_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(TOL_COMPILE_DEFINITIONS ${TOL}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${TOL}Config.cmake.in
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}Config.cmake
    INSTALL_DESTINATION 
    ${${TOL}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${TOL}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tools Installation

# Install *Config and *ConfigVersion files
install(FILES ${${TOL}_CONFIG_BUILD_DIR}/${TOL}Config.cmake 
    ${${TOL}_CONFIG_BUILD_DIR}/${TOL}ConfigVersion.cmake
    DESTINATION ${${TOL}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${TOL}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${TOL}-lib
    DESTINATION ${${TOL}_LIB_INSTALL_DIR}
    EXPORT ${TOL}Targets )

# Export library into *Targets file
install(EXPORT ${TOL}Targets
    DESTINATION ${${TOL}_CONFIG_INSTALL_DIR} )

