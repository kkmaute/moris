# Core Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    core.hpp
    banner.hpp
    common.hpp
    typedefs.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Library
add_library(${COR}-lib INTERFACE)
target_include_directories(${COR}-lib INTERFACE .)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Core Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${COR} )

