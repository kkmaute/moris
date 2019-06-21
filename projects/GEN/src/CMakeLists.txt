# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_GE_Analytic.hpp
    cl_GE_Core.hpp
    cl_GE_Discrete.hpp
    cl_GE_Element.hpp
    cl_GE_Enums.hpp
    cl_GE_Factory.hpp
    cl_GE_Geometry_Library.hpp
    cl_GE_Geometry.hpp
  	cl_GE_Node.hpp
  	cl_GE_Nodal_Info.hpp
  	cl_GE_Primitive.hpp
 	cl_GE_SDF.hpp 
    )

# List header dependencies
set(GEN_INCLUDES
    #${GEN}
    #${COM}
    #${FEM}
    #${FEM}/${INT}  
    #${HMR}
    #${IOS}
    #${LINALG}
    #${MTK}
    #${MTK}/${MAP} 
    #${SDF}
    #${TOL}
    )

foreach(GEN_INCLUDE ${GEN_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${GEN_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library

# List library source files
set(LIB_SOURCES
    cl_GE_Factory.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${LINALG}-lib
	${COM}-lib
	${INT}-lib
	${MAP}-lib
	${SDF}-lib
	#${CHR}-lib
    #${MTK}-lib
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    #${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${GEN}-lib STATIC ${LIB_SOURCES})
target_include_directories(${GEN}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries(${GEN}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${GEN}-lib PROPERTIES OUTPUT_NAME ${GEN})

foreach(TPL ${GEN_LIB_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${GEN}-lib
    	PUBLIC ${tpl})
    
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

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${GEN}Config.cmake.in
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${GEN}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Installation

# Install *Config and *ConfigVersion files
install(FILES ${${GEN}_CONFIG_BUILD_DIR}/${GEN}Config.cmake 
    ${${GEN}_CONFIG_BUILD_DIR}/${GEN}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${GEN} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${GEN} )

# Install library and executable
install(TARGETS ${GEN}-lib
    EXPORT ${GEN}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${GEN})

# Export library in to *Targets file
install(EXPORT ${GEN}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${GEN} )
