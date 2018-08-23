# src-----------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MSI}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_Adof.hpp
    cl_Dof_Manager.hpp
    cl_Dof_Type_Enums.hpp
    cl_Model_Solver_Interface.hpp
    cl_Equation_Object.hpp
    cl_Node_Obj.hpp
    cl_Pdof_Host.hpp
    cl_MSI_Solver_Interface.hpp
    )

# List header dependencies
set(MSI_INCLUDES
#     ${LNA}
      ${COM}
#     ${IOS}
#     ${LNA}
#     ${IOS}
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
# project Library

# List library source files
set(LIB_SOURCES
    cl_Dof_Manager.cpp
<<<<<<< HEAD
=======
    cl_Pdof_Host.cpp
>>>>>>> 3f85f560ee1d084502bc9c8d0eeba462b14145d0
    cl_Model_Solver_Interface.cpp
    cl_Equation_Object.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
#     ${COM}-lib
      ${DLA}-lib
#     ${MORIS_LDLIBS}
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
