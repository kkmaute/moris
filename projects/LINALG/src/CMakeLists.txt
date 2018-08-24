# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${LINALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	cl_LINALG_Enums.hpp
	cl_LINALG_Parameters.hpp
	cl_Dense_Backend_Matrix_Factory.hpp
	cl_Matrix_Base.hpp
	cl_Matrix.hpp	)
	
if(MORIS_USE_EIGEN)
list(APPEND HEADERS
	Eigen_Impl/cl_Matrix_Base_Eigen.hpp)
endif()

if(MORIS_USE_ARMA)
list(APPEND HEADERS)
endif()


# List header dependencies
set(LINALG_INCLUDES
    ${COM}
    ${LNA})

foreach(LINALG_INCLUDE ${LINALG_INCLUDES})
    include_directories(../../${LINALG_INCLUDE}/src)
endforeach()

# List dependencies

include_directories(${LINALG_IMPLEMENTATION_INCLUDES})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Executable

# LINALG Executable
# List source files
set(EXE_SOURCES
        main.cpp)

# List library dependencies
set(EXE_DEPENDENCIES 
	${MORIS_ARMADILLO_EIGEN_LIBS} 
	${COM}-lib)

# Set the output path for library
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${bin})

# Create library
add_executable(${LINALG} ${EXE_SOURCES})
target_link_libraries(${LINALG} ${EXE_DEPENDENCIES})




include_directories(${MORIS_PACKAGE_DIR}/XTK/include)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${LINALG} )
