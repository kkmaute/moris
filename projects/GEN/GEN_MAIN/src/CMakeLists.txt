# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN_MAIN}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	additional/cl_GEN_Basis_Function.hpp
	additional/cl_GEN_Edge_Topology.hpp
	additional/cl_GEN_Enums.hpp
	additional/cl_GEN_Hexahedron_8_Topology.hpp
	additional/cl_GEN_Hexahedron_8_Basis_Function.hpp	
	additional/cl_GEN_Interpolation.hpp
	additional/cl_GEN_Linear_Basis_Functions.hpp
	additional/cl_GEN_Matrix_Base_Utilities.hpp
	additional/cl_GEN_Pending_Node.hpp
	additional/cl_GEN_Phase_Table.hpp
	additional/cl_GEN_Quad_4_Basis_Function.hpp
	additional/cl_GEN_Quad_4_Topology.hpp
	additional/cl_GEN_Topology.hpp
	additional/fn_GEN_approximatate.hpp
	additional/fn_GEN_bubble_sort.hpp
    
	geomeng/cl_GEN_Design_Variable_Interface.hpp
	gemoeng/cl_GEN_Field.hpp
	geomeng/cl_GEN_Geometry_Engine.hpp
	geomeng/cl_GEN_Geometry_Onject_Manager.hpp
	geomeng/cl_GEN_Geometry_Object.hpp
	geomeng/cl_GEN_Pdv.hpp
	geomeng/cl_GEN_Pdv_Host.hpp
	geomeng/cl_GEN_Pdv_Host_Manager.hpp
	geomeng/fn_GEN_triangle_geometry.hpp

    geometry/cl_GEN_Geometry_Analytic.hpp
	geometry/cl_GEN_Geometry_Discrete.hpp
	geometry/cl_GEN_Circle.hpp
#	geometry/cl_GEN_Cylinder_With_End_Caps.hpp
#	geometry/cl_GEN_Discrete_Level_Set.hpp
#	geometry/cl_GEN_Geom_Data.hpp
#	geometry/cl_GEN_Geom_Field_HMR.hpp
#	geometry/cl_GEN_Geom_Field.hpp
#	geometry/cl_GEN_Sphere.hpp
#	geometry/fn_bounding_box.hpp
	geometry/fn_GEN_create_geometry.hpp
	
	property/cl_GEN_Property.hpp
	
#	cl_GE_Analytic.hpp
#	cl_GE_Core.hpp
#	cl_GE_Discrete.hpp
#	cl_GE_Element.hpp
#	cl_GE_Enums.hpp
#	cl_GE_Factory.hpp
#	cl_GE_Geometry_Library.hpp
#	cl_GE_Geometry.hpp
#	cl_GE_Intersection_Object_Line.hpp
#	cl_GE_Intersection_Object.hpp
#	cl_GE_Main.hpp
#	cl_GE_Nodal_Info.hpp
#	cl_GE_Node.hpp
#	cl_GE_SDF.hpp
	
    )
    
# List header dependencies                                                                                                                                                                                                                                                     
set(GEN_INCLUDES                                                                                                                                                                                                                                                               
    "additional"
    "geomeng"
    "geometry"
    "property"
    )

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
    list(APPEND CMAKE_CURRENT_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${GEN_INCLUDES}")
endforeach()


# include_directories(${MORIS_PACKAGE_DIR}/GEN/src/include)

#FIXME: Remove this once xtk has removed subdirectories from its package
set(XTK_INCLUDES
    "assert"
    "containers"
    "core"
    "geomeng"
    "geometry"
    "ios"
    "linalg"
    "mesh"
    "tools"
    "topology"
    "xtk"
    )

# Include depedency directories
foreach(XTK_INCLUDE ${XTK_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})
endforeach()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library

# List library source files
set(LIB_SOURCES
    cl_GE_Factory.cpp
    cl_GE_PDV_Info.cpp

	additional/cl_GEN_Phase_Table.cpp
    
	geomeng/cl_GEN_Geometry_Engine.cpp
	geomeng/cl_GEN_Geometry_Object.cpp

    geometry/cl_GEN_Geometry_Analytic.cpp
	geometry/cl_GEN_Geometry_Discrete.cpp
	geometry/fn_GEN_create_geometry.cpp
	geometry/cl_GEN_Circle.cpp
	
	property/cl_GEN_Property.cpp
	
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${COM}-lib
	${INT}-lib
	${LINALG}-lib
	${MAP}-lib
	${SDF}-lib
	${XTK}-lib
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

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
	string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
	string(APPEND GEN_CONFIG_MORIS_DEPENDENCIES
		"include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
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
