# LINALG --------------------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${LINALG}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
	cl_Matrix_Base.hpp
	cl_LINALG_Enums.hpp
	cl_Dense_Backend_Matrix_Factory.hpp)

# List header dependencies
set(LINALG_INCLUDES
    ${ARMADILLO_EIGEN} )

foreach(PROJ_INLCUDE ${LINALG_INCLUDES})
    include_directories(../../${LINALG_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Linear Algebra Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${LINALG} )
