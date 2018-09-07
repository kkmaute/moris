# Model Solver Interface Source -------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MSI}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_MSI_Adof.hpp
    cl_MSI_Dof_Manager.hpp
    cl_MSI_Dof_Type_Enums.hpp
    cl_MSI_Model_Solver_Interface.hpp
    cl_MSI_Equation_Object.hpp
    cl_MSI_Node_Obj.hpp
    cl_MSI_Pdof_Host.hpp
    cl_MSI_Solver_Interface.hpp
    cl_MSI_Node.hpp
    )

# List header dependencies
set(MSI_INCLUDES
#     ${LNA}
      ${COM}
#     ${LNA}
#     ${DLA}
      ${FEM}/${INT}
    )

foreach(MSI_INLCUDE ${MSI_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MSI_INCLUDE}/src)
endforeach()

include_directories(${MORIS_PACKAGE_DIR}/LNA/src)
include_directories(${MORIS_PACKAGE_DIR}/IOS/src)
include_directories(${MORIS_PACKAGE_DIR}/DLA/src)
include_directories(${MORIS_PACKAGE_DIR}/COM/src)
include_directories(${MORIS_PACKAGE_DIR}/MTK/src)
include_directories(${MORIS_PACKAGE_DIR}/FEM/INT/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Library

# List library source files
set(LIB_SOURCES
    cl_MSI_Dof_Manager.cpp
    cl_MSI_Pdof_Host.cpp
    cl_MSI_Model_Solver_Interface.cpp
    cl_MSI_Equation_Object.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
#     ${COM}-lib
      ${DLA}-lib
#     ${MORIS_BASE_LIBS}
#     ${MORIS_PETSC_LIBS}
#     ${MORIS_T_LIBS}
      ${MORIS_T_LIBS}
   )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MSI}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${MSI}-lib ${LIB_DEPENDENCIES})
set_target_properties(${MSI}-lib PROPERTIES OUTPUT_NAME ${MSI})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${MSI}Config.cmake.in
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MSI}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Model Solver Interface Installation

# Install *Config and *ConfigVersion files
install(FILES ${${MSI}_CONFIG_BUILD_DIR}/${MSI}Config.cmake 
    ${${MSI}_CONFIG_BUILD_DIR}/${MSI}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${MSI} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MSI} )

# Install library and executable
install(TARGETS ${MSI}-lib
    EXPORT ${MSI}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MSI}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${MSI}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${MSI} )

# Export library in to *Targets file
install(EXPORT ${MSI}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MSI} )
