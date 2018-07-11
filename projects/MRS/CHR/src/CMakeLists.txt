# Chronos Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    chronos.hpp
    cl_Stopwatch.hpp
    cl_Waitbar.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${CHR} )

