# MTK Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    cl_Database.hpp
    cl_Mesh_Enums.hpp
    cl_Mesh.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MTK Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MTK Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MTK Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# MTK Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MTK} )

