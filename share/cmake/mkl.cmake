# -------------------------------------------------------------------------
# MKL libraries -----------------------------------------------------------
# -------------------------------------------------------------------------

find_package(MKL)

message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_MKL")
list(APPEND MORIS_INCDIRS ${MKL_INCLUDE_DIRS})
list(APPEND MORIS_LDLIBS ${MKL_LIBRARIES})
