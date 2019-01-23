# MAP Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MAP}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    )

# List header dependencies
set(MAP_INCLUDES
    ${LINALG}
    ${COM}
    ${SOL}/${DLA}
    ${SOL}/${NLA}
    ${MTK}
    ${HMR}
    ${MTK}/${MAP}
    ${FEM}/${INT}
    ${FEM}/${MDL}
    ${FEM}/${MSI}
    )
    
# Include depedency directories
foreach(MAP_INCLUDE ${MAP_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MAP_INCLUDE}/src)
endforeach()
    
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MAP Library

# List library source files
set( LIB_SOURCES
    cl_MTK_Mapper.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MTK}-lib
    ${COM}-lib
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_BOOST_LIBS}
    )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MAP}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${MAP}-lib ${LIB_DEPENDENCIES})
set_target_properties(${MAP}-lib PROPERTIES OUTPUT_NAME ${MAP})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# map Config Files

# Build *Config and *ConfigVersion files
#include(CMakePackageConfigHelpers)

#configure_package_config_file(
#    ${MAP}Config.cmake.in
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}Config.cmake
#    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MAP}
#    PATH_VARS LIB_INSTALL_DIR )

#write_basic_package_version_file(
#    ${${MAP}_CONFIG_BUILD_DIR}/${MAP}ConfigVersion.cmake
#    COMPATIBILITY AnyNewerVersion )
    
    