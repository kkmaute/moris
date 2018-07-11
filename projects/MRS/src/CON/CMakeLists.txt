# Containers Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    cl_Array.hpp
    cl_Bi_Map.hpp
    cl_Bitbool.hpp
    cl_Bitset.hpp
    cl_BoostBitset.hpp
    cl_Cell.hpp
    cl_Dist_Map.hpp
    cl_Map.hpp
    cl_Param_List.hpp
    cl_Tuple.hpp
    containers.hpp
    fn_zip.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${CON} )

