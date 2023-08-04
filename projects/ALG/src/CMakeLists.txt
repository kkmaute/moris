#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Algorithms Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${ALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    algorithms.hpp
    fn_equal_to.hpp
    fn_equal_to.tpp
	fn_check_equal.hpp
    fn_hash_value.hpp
    fn_num_digits.hpp
    fn_clip_value.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Library

# List source files
set(LIB_SOURCES
    fn_equal_to.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library and set properties
add_library(${ALG}-lib STATIC ${LIB_SOURCES})
target_include_directories(${ALG}-lib PUBLIC 
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${ALG}_LIB_INSTALL_DIR}>
	$<INSTALL_INTERFACE:${${ALG}_HEADER_INSTALL_DIR}>
	)
target_link_libraries(${ALG}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${ALG}-lib PROPERTIES OUTPUT_NAME ${ALG})

# Link third party libraries
foreach(TPL ${ALG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${ALG}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    # Later versions of CMake (3.11+) allow include directories and
    # compile definitions to be included in an imported library. If this
    # system ever gets updated to a newer version of CMake, shifting both
    # to the imported library would be good practice. Imported library
    # targets were added to add inter-tpl dependencies where needed
    # (PETSc needs SuperLU, for example).
    # 	-- Matthew Ryan, June 18, 2019
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${ALG}-lib 
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${ALG}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${ALG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND ALG_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND ALG_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	#string(APPEND ALG_CONFIG_MORIS_DEPENDENCIES
	#	"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )

    string(APPEND ALG_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(ALG_COMPILE_DEFINITIONS ${ALG}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${ALG}Config.cmake.in
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}Config.cmake
    INSTALL_DESTINATION 
    ${${ALG}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${ALG}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Algorithms Installation

# Install *Config and *ConfigVersion files
install(FILES ${${ALG}_CONFIG_BUILD_DIR}/${ALG}Config.cmake 
    ${${ALG}_CONFIG_BUILD_DIR}/${ALG}ConfigVersion.cmake
    DESTINATION ${${ALG}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${ALG}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${ALG}-lib
    DESTINATION ${${ALG}_LIB_INSTALL_DIR}
    EXPORT ${ALG}Targets )

# Export library into *Targets file
install(EXPORT ${ALG}Targets
    DESTINATION ${${ALG}_CONFIG_INSTALL_DIR} )

