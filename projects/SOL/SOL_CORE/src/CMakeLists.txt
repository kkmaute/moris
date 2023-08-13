#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Sol Core Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${SOL_CORE}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_SOL_Amesos_GenOp.hpp
    cl_SOL_Warehouse.hpp
	cl_SOL_Matrix_Vector_Factory.hpp
    cl_SOL_Dist_Vector.hpp
    cl_SOL_Dist_Matrix.hpp
    cl_SOL_Dist_Map.hpp
    cl_Communicator_Epetra.hpp
	cl_Map_Epetra.hpp
    cl_Sparse_Matrix_EpetraFECrs.hpp
    cl_Vector_Epetra.hpp )

if (${MORIS_HAVE_PETSC})
list( APPEND HEADERS
    "cl_Map_PETSc.hpp"
    "cl_MatrixPETSc.hpp"
    "cl_Vector_PETSc.hpp"
    )
endif()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Library

# List library source files
set(LIB_SOURCES
    cl_SOL_Amesos_GenOp.cpp
    cl_SOL_Warehouse.cpp
	cl_SOL_Matrix_Vector_Factory.cpp
	cl_SOL_Dist_Vector.cpp
    cl_Communicator_Epetra.cpp
    cl_Map_Epetra.cpp
    cl_Sparse_Matrix_EpetraFECrs.cpp
    cl_Vector_Epetra.cpp
)

if (${MORIS_HAVE_PETSC})
list( APPEND LIB_SOURCES
    "cl_Map_PETSc.cpp"
    "cl_MatrixPETSc.cpp"
    "cl_Vector_PETSc.cpp"
    )
endif()

# List library dependencies
set(LIB_DEPENDENCIES
    ${DLA}-lib
    ${ENM}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${SOL_CORE}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${SOL_CORE}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${SOL_CORE}_HEADER_INSTALL_DIR}> )
target_link_libraries(${SOL_CORE}-lib ${LIB_DEPENDENCIES})
set_target_properties(${SOL_CORE}-lib PROPERTIES OUTPUT_NAME ${SOL_CORE})

include(${CMAKE_SOURCE_DIR}/projects/tgt_includes.cmake)

foreach(TPL ${SOL_CORE_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${SOL_CORE}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${SOL_CORE}-lib
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${SOL_CORE}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${SOL_CORE_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND SOL_CORE_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND SOL_CORE_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(SOL_CORE_COMPILE_DEFINITIONS ${SOL_CORE}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${SOL_CORE}Config.cmake.in
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}Config.cmake
    INSTALL_DESTINATION 
    ${${SOL_CORE}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${SOL_CORE}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Installation

# Install *Config and *ConfigVersion files
install(FILES ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}Config.cmake 
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}ConfigVersion.cmake
    DESTINATION ${${SOL_CORE}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${SOL_CORE}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${SOL_CORE}-lib
    DESTINATION ${${SOL_CORE}_LIB_INSTALL_DIR}
    EXPORT ${SOL_CORE}Targets )

# Export library into *Targets file
install(EXPORT ${SOL_CORE}Targets
    DESTINATION ${${SOL_CORE}_CONFIG_INSTALL_DIR} )

