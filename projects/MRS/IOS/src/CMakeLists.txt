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

# Install *Config and *ConfigVersion files
install(FILES ${${IOS}_CONFIG_BUILD_DIR}/${IOS}Config.cmake 
    ${${IOS}_CONFIG_BUILD_DIR}/${IOS}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${IOS} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${IOS} )

