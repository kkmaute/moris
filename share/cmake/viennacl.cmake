# -------------------------------------------------------------------------
# ViennaCL libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(ViennaCL)

message(STATUS "VIENNACL_LIBRARIES: ${VIENNACL_LIBRARIES}")

list(APPEND MORIS_INCDIRS ${VIENNACL_INCLUDE_DIRS})

mark_as_advanced(ViennaCL_DIR)
