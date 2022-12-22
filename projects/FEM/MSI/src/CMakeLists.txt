#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Model Solver Interface Source -------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MSI}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_MSI_Adof.hpp
    cl_MSI_Dof_Manager.hpp
    cl_MSI_Model_Solver_Interface.hpp
    cl_MSI_Equation_Model.hpp
    cl_MSI_Equation_Object.hpp
    cl_MSI_Equation_Set.hpp
    cl_MSI_Pdof_Host.hpp
    cl_MSI_Solver_Interface.hpp
    cl_MSI_Node_Proxy.hpp
    cl_MSI_Element_Proxy.hpp
    cl_MSI_Multigrid.hpp
    cl_MSI_Design_Variable_Interface.hpp
    fn_MSI_get_mesh_index_for_dof_type.hpp
    )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Library

# List library source files
set(LIB_SOURCES
    cl_MSI_Dof_Manager.cpp
    cl_MSI_Pdof_Host.cpp
    cl_MSI_Model_Solver_Interface.cpp
    cl_MSI_Equation_Object.cpp
    cl_MSI_Equation_Set.cpp
    cl_MSI_Equation_Model.cpp
    cl_MSI_Multigrid.cpp
    cl_MSI_Solver_Interface.cpp
    cl_MSI_Design_Variable_Interface.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${COM}-lib
    ${CON}-lib
    ${INT}-lib
    ${DLA}-lib
    ${PRM}-lib
    ${ENM}-lib
    #${MDL}-lib
   )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MSI}-lib STATIC ${LIB_SOURCES})
target_include_directories(${MSI}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${MSI}_HEADER_INSTALL_DIR}> )
target_link_libraries(${MSI}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MSI}-lib PROPERTIES OUTPUT_NAME ${MSI})

foreach(TPL ${MSI_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MSI}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MSI}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MSI}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${MSI_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND MSI_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND MSI_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(MSI_COMPILE_DEFINITIONS ${MSI}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MSI}Config.cmake.in
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}Config.cmake
    INSTALL_DESTINATION 
    ${${MSI}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${MSI}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MSI}_CONFIG_BUILD_DIR}/${MSI}Config.cmake 
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}ConfigVersion.cmake
    DESTINATION ${${MSI}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${MSI}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${MSI}-lib
    DESTINATION ${${MSI}_LIB_INSTALL_DIR}
    EXPORT ${MSI}Targets )

# Export library into *Targets file
install(EXPORT ${MSI}Targets
    DESTINATION ${${MSI}_CONFIG_INSTALL_DIR} )

