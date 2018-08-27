# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${LINALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	cl_Matrix.hpp
	linalg_typedefs.hpp
	fn_trans.hpp
	op_times.hpp
	fn_comp_abs.hpp
	fn_det.hpp
	fn_print.hpp	)
	
if(MORIS_USE_EIGEN)
list(APPEND HEADERS
	Eigen_Impl/cl_Matrix_Eigen_3x3.hpp
	Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp
	Eigen_Impl/fn_comp_abs_Eigen.hpp
	Eigen_Impl/fn_det_Eigen.hpp
	Eigen_Impl/op_times_Eigen.hpp
	Eigen_Impl/fn_trans_Eigen.hpp)
endif()

if(MORIS_USE_ARMA)
list(APPEND HEADERS
	Arma_Impl/cl_Matrix_Arma_Dynamic.hpp
	Arma_Impl/op_times_Arma.hpp
	Arma_Impl/fn_comp_abs_Arma.hpp
	Arma_Impl/fn_det_Arma.hpp
	Arma_Impl/fn_trans_Arma.hpp)
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
	${MORIS_LDLIBS} 
	${COM}-lib
	${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_SUPERLU_LIBS})
	
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
