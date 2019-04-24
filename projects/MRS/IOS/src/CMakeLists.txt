# IOS --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    ios.hpp
    cl_Logger.hpp
    fn_Logger.hpp
    fn_to_stdio.hpp
    IO_Tools.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Installation
# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${IOS} )

