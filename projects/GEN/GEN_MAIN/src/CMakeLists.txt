# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN_MAIN}_VERSION ${MORIS_VERSION})

# List header dependencies                                                                                                                                                                                                                                                     
set(GEN_INCLUDES
	"field"
	"field/geometry"
	"field/property"
	"pdv"
	"topology"
	)

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})
	include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
	list(APPEND CMAKE_CURRENT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${GEN_INCLUDES}")
endforeach()

# Include WRK directory
include_directories(${MORIS_PACKAGE_DIR}/WRK/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library
set(HEADERS

	cl_GEN_Geometry_Engine.hpp
	cl_GEN_Phase_Table.hpp
	cl_GEN_Geometric_Proximity.hpp
    st_GEN_Geometry_Engine_Parameters.hpp

	field/cl_GEN_Field.hpp
	field/cl_GEN_Field_Analytic.hpp
	field/cl_GEN_Field_Discrete_Integration.hpp
	field/cl_GEN_Field_Discrete_Interpolation.hpp
	field/cl_GEN_User_Defined_Field.hpp
	field/cl_GEN_BSpline_Field.hpp

	field/geometry/cl_GEN_Geometry.hpp
	field/geometry/cl_GEN_Circle.hpp
	field/geometry/cl_GEN_Superellipse.hpp
	field/geometry/cl_GEN_Sphere.hpp
	field/geometry/cl_GEN_Superellipsoid.hpp
	field/geometry/cl_GEN_Plane.hpp
	field/geometry/cl_GEN_Mesh_Field_Geometry.hpp
	field/geometry/cl_GEN_Voxel_Input.hpp
	field/geometry/cl_GEN_Single_Grain.hpp
	field/geometry/cl_GEN_User_Defined_Geometry.hpp
	field/geometry/cl_GEN_BSpline_Geometry.hpp
	field/geometry/cl_GEN_Stored_Geometry.hpp
	field/geometry/cl_GEN_Geometry_SDF.hpp
	field/geometry/cl_GEN_Multigeometry.hpp
	field/geometry/cl_GEN_Swiss_Cheese_Slice.hpp
	field/geometry/fn_GEN_create_geometries.hpp
	field/geometry/cl_GEN_Child_Node.hpp
    field/geometry/st_GEN_Geometry_Parameters.hpp
    

	field/property/cl_GEN_Property.hpp
	field/property/cl_GEN_Constant_Property.hpp
	field/property/cl_GEN_Scaled_Field.hpp
	field/property/cl_GEN_User_Defined_Property.hpp
	field/property/cl_GEN_BSpline_Property.hpp
	field/property/fn_GEN_create_properties.hpp
    field/property/st_GEN_Property_Parameters.hpp
    field/st_GEN_Field_Parameters.hpp

	pdv/cl_GEN_Pdv_Host_Manager.hpp
	pdv/cl_GEN_Interpolation_Pdv_Host.hpp
	pdv/cl_GEN_Pdv.hpp
	pdv/cl_GEN_Pdv_Value.hpp
	pdv/cl_GEN_Pdv_Property.hpp
	pdv/cl_GEN_Intersection_Node.hpp
	pdv/cl_GEN_Intersection_Node_Linear.hpp
	pdv/cl_GEN_Intersection_Node_Bilinear.hpp
    )
# List library source files
set(LIB_SOURCES

	cl_GEN_Geometry_Engine.cpp
	cl_GEN_Phase_Table.cpp
	cl_GEN_Geometric_Proximity.cpp

	field/cl_GEN_Field.cpp
	field/cl_GEN_Field_Analytic.cpp
	field/cl_GEN_Field_Discrete_Integration.cpp
	field/cl_GEN_Field_Discrete_Interpolation.cpp
	field/cl_GEN_User_Defined_Field.cpp
	field/cl_GEN_BSpline_Field.cpp

	field/geometry/cl_GEN_Geometry.cpp
	field/geometry/cl_GEN_Circle.cpp
	field/geometry/cl_GEN_Superellipse.cpp
	field/geometry/cl_GEN_Sphere.cpp
	field/geometry/cl_GEN_Superellipsoid.cpp
	field/geometry/cl_GEN_Plane.cpp
	field/geometry/cl_GEN_Mesh_Field_Geometry.cpp
	field/geometry/cl_GEN_Voxel_Input.cpp
	field/geometry/cl_GEN_Single_Grain.cpp
	field/geometry/cl_GEN_User_Defined_Geometry.cpp
	field/geometry/cl_GEN_BSpline_Geometry.cpp
	field/geometry/cl_GEN_Stored_Geometry.cpp
	field/geometry/cl_GEN_Geometry_SDF.cpp
	field/geometry/cl_GEN_Multigeometry.cpp
	field/geometry/cl_GEN_Swiss_Cheese_Slice.cpp
	field/geometry/fn_GEN_create_geometries.cpp

	field/property/cl_GEN_Property.cpp
	field/property/cl_GEN_Constant_Property.cpp
	field/property/cl_GEN_Scaled_Field.cpp
	field/property/cl_GEN_User_Defined_Property.cpp
	field/property/cl_GEN_BSpline_Property.cpp
	field/property/fn_GEN_create_properties.cpp

	pdv/cl_GEN_Pdv_Host_Manager.cpp
	pdv/cl_GEN_Interpolation_Pdv_Host.cpp
	pdv/cl_GEN_Pdv.cpp
	pdv/cl_GEN_Pdv_Value.cpp
	pdv/cl_GEN_Pdv_Property.cpp
	pdv/cl_GEN_Child_Node.cpp
	pdv/cl_GEN_Intersection_Node.cpp
	pdv/cl_GEN_Intersection_Node_Linear.cpp
	pdv/cl_GEN_Intersection_Node_Bilinear.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${COM}-lib
	${LINALG}-lib
	${MAP}-lib
	${SDF}-lib
 	${XTK}-lib
    #${WRK}-lib
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${GEN_MAIN}-lib STATIC ${LIB_SOURCES})
target_include_directories(${GEN_MAIN}-lib PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_SOURCE_DIR}/geometry;>"
	$<INSTALL_INTERFACE:${${GEN_MAIN}_HEADER_INSTALL_DIR}>
	$<INSTALL_INTERFACE:${HEADER_INSTALL_DIR}/${MORIS}> )
target_link_libraries(${GEN_MAIN}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${GEN_MAIN}-lib PROPERTIES OUTPUT_NAME ${GEN_MAIN})

foreach(TPL ${GEN_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${GEN_MAIN}-lib
    	PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${GEN_MAIN}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${GEN_MAIN}-lib
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

    string(APPEND GEN_MAIN_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} QUIET HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	#string(APPEND GEN_CONFIG_MORIS_DEPENDENCIES
	#	"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )

        string(APPEND GEN_MAIN_CONFIG_MORIS_DEPENDENCIES
		"find_package(${LIB} HINTS \"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/\")\n" )
endforeach()

# Properties to be coppied into config file
get_target_property(GEN_COMPILE_DEFINITIONS ${GEN_MAIN}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)


configure_package_config_file(
    ${GEN_MAIN}Config.cmake.in
    ${${GEN_MAIN}_CONFIG_BUILD_DIR}/${GEN_MAIN}Config.cmake
    INSTALL_DESTINATION 
    ${${GEN_MAIN}_CONFIG_INSTALL_DIR}
    PATH_VARS 
    ${GEN_MAIN}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${GEN_MAIN}_CONFIG_BUILD_DIR}/${GEN_MAIN}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Installation

# Install *Config and *ConfigVersion files
install(FILES ${${GEN_MAIN}_CONFIG_BUILD_DIR}/${GEN_MAIN}Config.cmake 
    ${${GEN_MAIN}_CONFIG_BUILD_DIR}/${GEN_MAIN}ConfigVersion.cmake
    DESTINATION ${${GEN_MAIN}_CONFIG_INSTALL_DIR} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${${GEN_MAIN}_HEADER_INSTALL_DIR} )

# Install library
install(TARGETS ${GEN_MAIN}-lib
    DESTINATION ${${GEN_MAIN}_LIB_INSTALL_DIR}
    EXPORT ${GEN_MAIN}Targets )

# Export library into *Targets file
install(EXPORT ${GEN_MAIN}Targets
    DESTINATION ${${GEN_MAIN}_CONFIG_INSTALL_DIR} )
