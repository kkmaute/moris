# Geometry Engine Source --------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${GEN}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    SDF/cl_Ge_SDF_Core.hpp
    SDF/cl_Ge_SDF_Data.hpp
    SDF/cl_Ge_SDF_Gen.hpp
    SDF/cl_Ge_SDF_Mesh_Data.hpp
    SDF/cl_Ge_SDF_Triangle_File.hpp
    SDF/cl_Ge_SDF_Triangle.hpp
    SDF/GeUtilities.hpp)

# List header dependencies
set(GEN_INCLUDES
    ${XTK}
    ${TOL}
    ${LINALG}
    ${XTK}
    ${MTK}
    ${ALG}
    ${COM}
    ${STK}
    ${DLA}
#    ${MOD} // < -- this library is obsolete
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
    #SDF/cl_Ge_SDF_Core.cpp
    #SDF/cl_Ge_SDF_Data.cpp
    #SDF/cl_Ge_SDF_Gen.cpp
    #SDF/cl_Ge_SDF_Mesh_Data.cpp
    #SDF/cl_Ge_SDF_Triangle_File.cpp
    #SDF/cl_Ge_SDF_Triangle.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    "" )

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
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${GEN}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${GEN}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${GEN} )

# Export library in to *Targets file
install(EXPORT ${GEN}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${GEN} )
