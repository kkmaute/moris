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
    ${LNA}
    ${COM}
    ${IOS}
    ${LNA}
    ${IOS}
    ${DLA})

foreach(MSI_INLCUDE ${MSI_INCLUDES})
    include_directories(../../../${MSI_INCLUDE}/src)
endforeach()

include_directories(../../../LNA/src)
include_directories(../../../IOS/src)
include_directories(../../../DLA/src)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# project Library

# List library source files
set(LIB_SOURCES
    cl_Model_Solver_Interface.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${COM}-lib
    ${DLA}-lib
    ${MORIS_PETSC_LIBS}
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
