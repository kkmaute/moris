#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(GEN_INCLUDES
	"ADV"
	"ADV/fields"
	"PDV"
	)

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/GEN/src/${GEN_INCLUDES})
	list(APPEND CMAKE_CURRENT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${GEN_INCLUDES}")
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library
set(HEADERS

	# GEN
	cl_GEN_Geometry_Engine.hpp
	cl_GEN_Geometric_Query_Interface.hpp
	cl_GEN_Phase_Table.hpp
	cl_GEN_Node_Manager.hpp
	cl_GEN_Geometric_Proximity.hpp
    st_GEN_Geometry_Engine_Parameters.hpp

	# ADV
	ADV/cl_GEN_ADV.hpp
    ADV/cl_GEN_ADV_Manager.hpp
	ADV/cl_GEN_Design_Factory.hpp
	ADV/cl_GEN_Design_Field.hpp
	ADV/cl_GEN_Geometry.hpp
	ADV/cl_GEN_Level_Set_Geometry.hpp
	ADV/cl_GEN_Property.hpp

	# ADV/fields
	ADV/fields/cl_GEN_Field.hpp
	ADV/fields/cl_GEN_Field_Analytic.hpp
	ADV/fields/cl_GEN_Field_Discrete_Integration.hpp
	ADV/fields/cl_GEN_Field_Discrete_Interpolation.hpp
	ADV/fields/cl_GEN_Constant_Field.hpp
	ADV/fields/cl_GEN_Scaled_Field.hpp
	ADV/fields/cl_GEN_BSpline_Field.hpp
	ADV/fields/cl_GEN_Circle.hpp
	ADV/fields/cl_GEN_Superellipse.hpp
	ADV/fields/cl_GEN_Sphere.hpp
	ADV/fields/cl_GEN_Superellipsoid.hpp
	ADV/fields/cl_GEN_Plane.hpp
	ADV/fields/cl_GEN_Mesh_Field.hpp
	ADV/fields/cl_GEN_Voxel_Input.hpp
	ADV/fields/cl_GEN_Single_Grain.hpp
	ADV/fields/cl_GEN_Stored_Field.hpp
	ADV/fields/cl_GEN_Signed_Distance_Field.hpp
	ADV/fields/cl_GEN_Image_Signed_Distance_Field.hpp
	ADV/fields/cl_GEN_Combined_Fields.hpp
	ADV/fields/cl_GEN_User_Defined_Field.hpp
	ADV/fields/fn_GEN_create_field.hpp
	ADV/fields/cl_GEN_Field_Array_Factory.hpp

	# PDV
	PDV/cl_GEN_Pdv_Host_Manager.hpp
	PDV/cl_GEN_Child_Node.hpp
	PDV/cl_GEN_Interpolation_Pdv_Host.hpp
	PDV/cl_GEN_Pdv.hpp
	PDV/cl_GEN_Pdv_Value.hpp
	PDV/cl_GEN_Pdv_Property.hpp
	PDV/cl_GEN_Node.hpp
	PDV/cl_GEN_Base_Node.hpp
	PDV/cl_GEN_Derived_Node.hpp
	PDV/cl_GEN_Intersection_Node.hpp
	PDV/cl_GEN_Intersection_Node_Level_Set.hpp
	PDV/cl_GEN_Intersection_Node_Linear.hpp
	PDV/cl_GEN_Intersection_Node_Bilinear.hpp
    )

# List library source files
set(LIB_SOURCES

	# GEN
	cl_GEN_Geometry_Engine.cpp
	cl_GEN_Node_Manager.cpp
	cl_GEN_Phase_Table.cpp
	cl_GEN_Geometric_Proximity.cpp

	# ADV
	ADV/cl_GEN_ADV.cpp
	ADV/cl_GEN_ADV_Manager.cpp
	ADV/cl_GEN_Design_Factory.cpp
	ADV/cl_GEN_Design_Field.cpp
	ADV/cl_GEN_Geometry.cpp
	ADV/cl_GEN_Level_Set_Geometry.cpp
	ADV/cl_GEN_Property.cpp

	# ADV/fields
	ADV/fields/cl_GEN_Field.cpp
	ADV/fields/cl_GEN_Field_Discrete_Integration.cpp
	ADV/fields/cl_GEN_Field_Discrete_Interpolation.cpp
	ADV/fields/cl_GEN_Constant_Field.cpp
	ADV/fields/cl_GEN_Scaled_Field.cpp
	ADV/fields/cl_GEN_BSpline_Field.cpp
	ADV/fields/cl_GEN_Circle.cpp
	ADV/fields/cl_GEN_Superellipse.cpp
	ADV/fields/cl_GEN_Sphere.cpp
	ADV/fields/cl_GEN_Superellipsoid.cpp
	ADV/fields/cl_GEN_Plane.cpp
	ADV/fields/cl_GEN_Mesh_Field.cpp
	ADV/fields/cl_GEN_Voxel_Input.cpp
	ADV/fields/cl_GEN_Single_Grain.cpp
	ADV/fields/cl_GEN_Stored_Field.cpp
	ADV/fields/cl_GEN_Signed_Distance_Field.cpp
	ADV/fields/cl_GEN_Image_Signed_Distance_Field.cpp
	ADV/fields/cl_GEN_Combined_Fields.cpp
	ADV/fields/cl_GEN_User_Defined_Field.cpp
	ADV/fields/fn_GEN_create_field.cpp
	ADV/fields/cl_GEN_Field_Array_Factory.cpp

	# PDV
	PDV/cl_GEN_Pdv_Host_Manager.cpp
	PDV/cl_GEN_Interpolation_Pdv_Host.cpp
	PDV/cl_GEN_Pdv.cpp
	PDV/cl_GEN_Pdv_Value.cpp
	PDV/cl_GEN_Pdv_Property.cpp
	PDV/cl_GEN_Node.cpp
	PDV/cl_GEN_Base_Node.cpp
	PDV/cl_GEN_Derived_Node.cpp
	PDV/cl_GEN_Child_Node.cpp
	PDV/cl_GEN_Intersection_Node.cpp
	PDV/cl_GEN_Intersection_Node_Level_Set.cpp
	PDV/cl_GEN_Intersection_Node_Linear.cpp
	PDV/cl_GEN_Intersection_Node_Bilinear.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${COM}-lib
	${LINALG}-lib
	${MAP}-lib
	${SDF}-lib
    ${ENM}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${GEN}-lib ${LIB_LINK_MODE} ${LIB_SOURCES})
target_include_directories(${GEN}-lib PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/geometry;>"
	$<INSTALL_INTERFACE:${${GEN}_HEADER_INSTALL_DIR}>
	$<INSTALL_INTERFACE:${HEADER_INSTALL_DIR}/${MORIS}> )
target_link_libraries(${GEN}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${GEN}-lib PROPERTIES OUTPUT_NAME ${GEN})

foreach(TPL ${GEN_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${GEN}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${GEN}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${GEN}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${GEN_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
	string(APPEND GEN_CONFIG_TPL_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	#string(APPEND XTK_CONFIG_MORIS_DEPENDENCIES
    #		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )

    string(APPEND GEN_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	#string(APPEND GEN_CONFIG_MORIS_DEPENDENCIES
	#	"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )

        string(APPEND GEN_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(GEN_COMPILE_DEFINITIONS ${GEN}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)


configure_package_config_file(
    ${GEN}Config.cmake.in
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}Config.cmake
    INSTALL_DESTINATION 
    ${${GEN}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${GEN}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Installation

# Install *Config and *ConfigVersion files
install(FILES ${${GEN}_CONFIG_BUILD_DIR}/${GEN}Config.cmake
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}ConfigVersion.cmake
    DESTINATION ${${GEN}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${GEN}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${GEN}-lib
    DESTINATION ${${GEN}_LIB_INSTALL_DIR}
    EXPORT ${GEN}Targets )

# Export library into *Targets file
install(EXPORT ${GEN}Targets
    DESTINATION ${${GEN}_CONFIG_INSTALL_DIR} )

