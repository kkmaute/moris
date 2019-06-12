# Exceptions Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    exceptions.hpp
    cl_invalid_size.hpp
    cl_out_of_range.hpp
    cl_tic_toc_error.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Target
add_library(${EXC}-target INTERFACE)
target_include_directories(${EXC}-target INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Exceptions Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${EXC} )

