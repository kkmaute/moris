# Sol Core Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${SOL_CORE}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_SOL_Warehouse.hpp
 )
    
# List header dependencies
set(SOL_CORE_INCLUDES
    ${LINALG}
    ${SOL}/${DLA}
    ${IOS}
    ${COM}
    ${SOL}/${NLA}
    ${FEM}/${MSI} )

# Include depedency directories
foreach(SOL_CORE_INCLUDE ${SOL_CORE_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${SOL_CORE_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Library

# List library source files
set(LIB_SOURCES
    cl_SOL_Warehouse.cpp
 )

# List library dependencies
set(LIB_DEPENDENCIES
    ${DLA}-lib
    ${COM}-lib
    ${MSI}-lib
    ${NLA}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS} )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${SOL_CORE}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${SOL_CORE}-lib ${LIB_DEPENDENCIES})
set_target_properties(${SOL_CORE}-lib PROPERTIES OUTPUT_NAME ${SOL_CORE})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${SOL_CORE}Config.cmake.in
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${SOL_CORE}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sol Core Installation

# Install *Config and *ConfigVersion files
install(FILES ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}Config.cmake 
    ${${SOL_CORE}_CONFIG_BUILD_DIR}/${SOL_CORE}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${SOL_CORE} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${SOL_CORE} )

# Install library and executable
install(TARGETS ${SOL_CORE}-lib
    EXPORT ${SOL_CORE}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${SOL_CORE}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${SOL_CORE}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${SOL_CORE} )

# Export library in to *Targets file
install(EXPORT ${SOL_CORE}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${SOL_CORE} )
