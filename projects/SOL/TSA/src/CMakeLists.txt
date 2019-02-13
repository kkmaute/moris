# Time Solver Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${TSA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_TSA_Time_Solver.hpp
    cl_TSA_Monolithic_Time_Solver.hpp
    cl_TSA_Staggered_Time_Solver.hpp
    cl_TSA_Solver_Interface_Proxy.hpp
    cl_TSA_Solver_Interface_Proxy2.hpp
    cl_TSA_Time_Solver_Factory.hpp
    cl_TSA_Time_Solver_Enums.hpp
 )
    
# List header dependencies
set(TSA_INCLUDES
    ${LINALG}
    ${SOL}/${DLA}
    ${SOL}/${NLA}
    ${SOL}/${TSA}
    ${IOS}
    ${COM}
    ${FEM}/${MSI} )

# Include depedency directories
foreach(TSA_INCLUDE ${TSA_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${TSA_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Library

# List library source files
set(LIB_SOURCES
    cl_TSA_Time_Solver.cpp
    cl_TSA_Solver_Interface_Proxy.cpp
    cl_TSA_Solver_Interface_Proxy2.cpp
    cl_TSA_Monolithic_Time_Solver.cpp
    cl_TSA_Staggered_Time_Solver.cpp
    cl_TSA_Time_Solver_Factory.cpp
 )

# List library dependencies
set(LIB_DEPENDENCIES
    ${DLA}-lib
    ${NLA}-lib
    ${COM}-lib
    ${MSI}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    ${MORIS_ARPACK_LIBS} )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${TSA}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${TSA}-lib ${LIB_DEPENDENCIES})
set_target_properties(${TSA}-lib PROPERTIES OUTPUT_NAME ${TSA})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${TSA}Config.cmake.in
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${TSA}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Time Solver Installation

# Install *Config and *ConfigVersion files
install(FILES ${${TSA}_CONFIG_BUILD_DIR}/${TSA}Config.cmake 
    ${${TSA}_CONFIG_BUILD_DIR}/${TSA}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${TSA} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${TSA} )

# Install library and executable
install(TARGETS ${TSA}-lib
    EXPORT ${TSA}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${TSA}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${TSA}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${TSA} )

# Export library in to *Targets file
install(EXPORT ${TSA}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${TSA} )
