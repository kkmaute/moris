# NonLinear Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${NLA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_NLA_Newton_Solver.hpp
    cl_NLA_Solver_Input_Test.hpp )
    

# List header dependencies
set(NLA_INCLUDES
    ${LINALG}
    ${DLA}
    ${IOS}
    ${COM}
)

# Include depedency directories
foreach(NLA_INCLUDE ${NLA_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${NLA_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Library

# List library source files
set(LIB_SOURCES
    cl_NLA_Newton_Solver.cpp
    cl_NLA_Solver_Input_Test.cpp)

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_PETSC_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS} )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${NLA}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${NLA}-lib ${LIB_DEPENDENCIES})
set_target_properties(${NLA}-lib PROPERTIES OUTPUT_NAME ${NLA})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${NLA}Config.cmake.in
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${NLA}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Installation

# Install *Config and *ConfigVersion files
install(FILES ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake 
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${NLA} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${NLA} )

# Install library and executable
install(TARGETS ${NLA}-lib
    EXPORT ${NLA}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${NLA}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${NLA}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${NLA} )

# Export library in to *Targets file
install(EXPORT ${NLA}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${NLA} )
