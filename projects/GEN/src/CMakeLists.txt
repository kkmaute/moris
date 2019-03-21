# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_GEN_Geometry_Engine.hpp
    cl_GE_Geometry.hpp
 	cl_GE_Enums.hpp
  	cl_GE_Factory.hpp
  	cl_GE_LS.hpp
 	cl_GE_SDF.hpp
 	cl_GE_Node.hpp
    cl_GE_Element.hpp
    cl_GE_Analytical.hpp
    cl_GE_Main.hpp
)

# List header dependencies
set(GEN_INCLUDES
    ${LINALG}
    ${TOL}
    ${IOS}
    ${MTK}
    ${COM}
    ${FEM}/${INT} # for L2
    )

foreach(GEN_INCLUDE ${GEN_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${GEN_INCLUDE}/src)
endforeach()

# Additional includes
include_directories(${MORIS_PACKAGE_DIR}/${STK}/src/Hierarchical)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Geometry Engine Library

# List library source files
set(LIB_SOURCES
    cl_GEN_Dummy.cpp
    cl_GE_Factory.cpp
    cl_GE_Main.cpp
    cl_GE_SDF.cpp
    cl_GEN_Geometry_Engine.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MTK}-lib
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_BASE_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${GEN}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${GEN}-lib ${LIB_DEPENDENCIES})
set_target_properties(${GEN}-lib PROPERTIES OUTPUT_NAME ${GEN})

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
