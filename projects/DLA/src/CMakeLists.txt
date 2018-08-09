# Distributed Linear Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${DLA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Communicator_Epetra.hpp
    cl_DistLinAlg_Enums.hpp
    cl_Linear_Solver_Amesos2.hpp
    cl_Linear_Solver_Amesos.hpp
    cl_Linear_Solver_Aztec.hpp
    cl_Linear_Solver.hpp
    cl_Linear_Solver_PETSc.hpp
    cl_Linear_Solver_Trilinos.hpp
    cl_Map_Class.hpp
    cl_Map_Epetra.hpp
    cl_Map_PETSc.hpp
    cl_MatrixPETSc.hpp
    cl_Matrix_Vector_Factory.hpp
    cl_Model_Solver_Interface_Solver.hpp
    cl_Solver_Factory.hpp
    cl_Solver_Input.hpp
    cl_Solver_Input_Test.hpp
    cl_Sparse_Matrix_EpetraFECrs.hpp
    cl_Sparse_Matrix.hpp
    cl_Vector_Epetra.hpp
    cl_Vector.hpp
    cl_VectorPETSc.hpp )

# List header dependencies
set(DLA_INCLUDES
    ${LNA}
    ${IOS}
    ${COM})

# Include depedency directories
foreach(DLA_INCLUDE ${DLA_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${DLA_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Library

# List library source files
set(LIB_SOURCES
    cl_Communicator_Epetra.cpp
    cl_Linear_Solver_Amesos2.cpp
    cl_Linear_Solver_Amesos.cpp
    cl_Linear_Solver_Aztec.cpp
    cl_Linear_Solver_PETSc.cpp
    cl_Linear_Solver_Trilinos.cpp
    cl_Map_Epetra.cpp
    cl_Map_PETSc.cpp
    cl_MatrixPETSc.cpp
    cl_Matrix_Vector_Factory.cpp
    cl_Model_Solver_Interface_Solver.cpp
    cl_Solver_Factory.cpp
    cl_Solver_Input_Test.cpp
    cl_Sparse_Matrix_EpetraFECrs.cpp
    cl_Vector_Epetra.cpp
    cl_VectorPETSc.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${MORIS_PETSC_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS} )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${DLA}-lib STATIC ${LIB_SOURCES})
target_link_libraries(${DLA}-lib ${LIB_DEPENDENCIES})
set_target_properties(${DLA}-lib PROPERTIES OUTPUT_NAME ${DLA})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${DLA}Config.cmake.in
#     ${${DLA}_CONFIG_BUILD_DIR}/${DLA}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${DLA}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${DLA}_CONFIG_BUILD_DIR}/${DLA}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# distlinalg Installation

# # Install *Config and *ConfigVersion files
# install(FILES ${${DLA}_CONFIG_BUILD_DIR}/${DLA}Config.cmake 
#     ${${DLA}_CONFIG_BUILD_DIR}/${DLA}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${DLA} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${DLA} )

# Install library and executable
install(TARGETS ${DLA}-lib
    EXPORT ${DLA}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${DLA}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${DLA}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${DLA} )

# Export library in to *Targets file
install(EXPORT ${DLA}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${DLA} )
