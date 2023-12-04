#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# MIG Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MIG}_VERSION ${MORIS_VERSION})


# List header files
set(HEADERS
cl_MIG_Periodic_2D.hpp
cl_MIG_Periodic_3D.hpp
cl_MIG_Mesh_Editor.hpp
cl_MIG.cpp
)

include_directories(${MORIS_PACKAGE_DIR}/include)

set(GEN_INCLUDES
	"field"
	"field/geometry"
	"field/property"
	"pdv"
	)

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
endforeach()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MIG Library

# List library source files
set( LIB_SOURCES
cl_MIG_Periodic_2D.cpp
cl_MIG_Periodic_3D.cpp
cl_MIG_Mesh_Editor.cpp
cl_MIG.cpp
		cl_MIG_Spatial_Indexer.cpp
        cl_MIG_Spatial_Indexer.h
		cl_MIG_Spatial_Indexer_BruteForce.cpp
		cl_MIG_Spatial_Indexer_BruteForce.h
)



# List library dependencies
set(LIB_DEPENDENCIES
    ${GEN_MAIN}-lib
	${COM}-lib
	${LINALG}-lib
	${HMR}-lib
    ${TOL}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MIG}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${MIG}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:${${MIG}_HEADER_INSTALL_DIR}> )
target_link_libraries(${MIG}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MIG}-lib PROPERTIES OUTPUT_NAME ${MIG})

foreach(TPL ${MIG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${MIG}-lib
    	PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MIG}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${MIG}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MIG Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MIG Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${MIG_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND MIG_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND MIG_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})

    string(APPEND MIG_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} REQUIRED HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(MIG_COMPILE_DEFINITIONS ${MIG}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MIG}Config.cmake.in
    ${${MIG}_CONFIG_BUILD_DIR}/${MIG}Config.cmake
    INSTALL_DESTINATION 
    ${${MIG}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${MIG}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${MIG}_CONFIG_BUILD_DIR}/${MIG}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MIG Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MIG}_CONFIG_BUILD_DIR}/${MIG}Config.cmake 
    ${${MIG}_CONFIG_BUILD_DIR}/${MIG}ConfigVersion.cmake
    DESTINATION ${${MIG}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${MIG}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${MIG}-lib
    DESTINATION ${${MIG}_LIB_INSTALL_DIR}
    EXPORT ${MIG}Targets )

# Export library into *Targets file
install(EXPORT ${MIG}Targets
    DESTINATION ${${MIG}_CONFIG_INSTALL_DIR} )

