#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Distributed Linear Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${DLA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_DLA_Eigen_Solver.hpp
    cl_DLA_Linear_Solver_Amesos2.hpp
    cl_DLA_Linear_Solver_Amesos.hpp
    cl_DLA_Linear_Solver_Aztec.hpp
    cl_DLA_Linear_Solver_Belos.hpp
    cl_DLA_Linear_Solver.hpp
    cl_DLA_Linear_Problem.hpp
    cl_DLA_Linear_System_Trilinos.hpp
    cl_DLA_Solver_Factory.hpp
    cl_DLA_Solver_Interface.hpp
    cl_DLA_Linear_Solver_Algorithm.hpp
    cl_DLA_Preconditioner_Trilinos.hpp
    cl_DLA_Geometric_Multigrid.hpp)

if(${MORIS_HAVE_PETSC})
list( APPEND HEADERS
    cl_DLA_Linear_System_PETSc.hpp
    cl_DLA_Linear_Solver_PETSc.hpp
    cl_DLA_Preconditioner_PETSc.hpp
    )
endif()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Library

# List library source files
set(LIB_SOURCES
    cl_DLA_Eigen_Solver.cpp
    cl_DLA_Linear_Solver_Amesos2.cpp
    cl_DLA_Linear_Solver_Amesos.cpp
    cl_DLA_Linear_Solver_Aztec.cpp
    cl_DLA_Linear_Solver_Belos.cpp
    cl_DLA_Linear_System_Trilinos.cpp
    cl_DLA_Linear_Solver.cpp
    cl_DLA_Linear_Problem.cpp
    cl_DLA_Geometric_Multigrid.cpp
    cl_DLA_Solver_Interface.cpp
    cl_DLA_Preconditioner_Trilinos.cpp
    cl_DLA_Solver_Factory.cpp)

if(${MORIS_HAVE_PETSC})
list( APPEND LIB_SOURCES
    cl_DLA_Linear_System_PETSc.cpp
    cl_DLA_Linear_Solver_PETSc.cpp
    cl_DLA_Preconditioner_PETSc.cpp
    )
endif()

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${COM}-lib
    ${MTK}-lib
    ${SOL_CORE}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${DLA}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${DLA}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${DLA}_HEADER_INSTALL_DIR}> )
target_link_libraries(${DLA}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${DLA}-lib PROPERTIES OUTPUT_NAME ${DLA})

foreach(TPL ${DLA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${DLA}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${DLA}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${DLA}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${DLA_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND DLA_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND DLA_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(DLA_COMPILE_DEFINITIONS ${DLA}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${DLA}Config.cmake.in
    ${${DLA}_CONFIG_BUILD_DIR}/${DLA}Config.cmake
    INSTALL_DESTINATION 
    ${${DLA}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${DLA}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${DLA}_CONFIG_BUILD_DIR}/${DLA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Installation

# Install *Config and *ConfigVersion files
install(FILES ${${DLA}_CONFIG_BUILD_DIR}/${DLA}Config.cmake 
    ${${DLA}_CONFIG_BUILD_DIR}/${DLA}ConfigVersion.cmake
    DESTINATION ${${DLA}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${DLA}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${DLA}-lib
    DESTINATION ${${DLA}_LIB_INSTALL_DIR}
    EXPORT ${DLA}Targets )

# Export library into *Targets file
install(EXPORT ${DLA}Targets
    DESTINATION ${${DLA}_CONFIG_INSTALL_DIR} )

