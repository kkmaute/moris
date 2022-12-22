
#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Enums Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${ENM}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_FEM_Enums.hpp
    cl_Mesh_Enums.hpp
    cl_MSI_Dof_Type_Enums.hpp
    cl_MTK_Enums.hpp
    cl_NLA_Nonlinear_Solver_Enums.hpp
    cl_SOL_Enums.hpp
    cl_TSA_Time_Solver_Enums.hpp
    cl_VIS_Output_Enums.hpp
    cl_XTK_Enums.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Enums Library

# List library source files
set(LIB_SOURCES
    cl_MTK_Enums.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${ASR}-lib
    ${CON}-lib 
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Enums Target

# Create library
add_library(${ENM}-lib STATIC ${LIB_SOURCES})
target_include_directories(${ENM}-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${${ENM}_HEADER_INSTALL_DIR}> )
target_link_libraries(${ENM}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${ENM}-lib PROPERTIES OUTPUT_NAME ${ENM})

foreach(TPL ${MAP_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${ENM}-lib
    	PUBLIC ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${ENM}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${ENM}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Enums Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Enums Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${ENM_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	list(APPEND ENM_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND ENM_CONFIG_MORIS_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be copied into config file
get_target_property(ENM_COMPILE_DEFINITIONS ${ENM}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${ENM}Config.cmake.in
    ${${ENM}_CONFIG_BUILD_DIR}/${ENM}Config.cmake
    INSTALL_DESTINATION 
    ${${ENM}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${ENM}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${ENM}_CONFIG_BUILD_DIR}/${ENM}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Enums Installation

# Install *Config and *ConfigVersion files
install(FILES ${${ENM}_CONFIG_BUILD_DIR}/${ENM}Config.cmake 
    ${${ENM}_CONFIG_BUILD_DIR}/${ENM}ConfigVersion.cmake
    DESTINATION ${${ENM}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${ENM}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${ENM}-lib
    DESTINATION ${${ENM}_LIB_INSTALL_DIR}
    EXPORT ${ENM}Targets )

# Export library in to *Targets file
install(EXPORT ${ENM}Targets
    DESTINATION ${${ENM}_CONFIG_INSTALL_DIR} )

