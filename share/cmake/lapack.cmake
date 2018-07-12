# -------------------------------------------------------------------------
# LAPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(LAPACK)

message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_LAPACK")
list(APPEND MORIS_LDLIBS ${LAPACK_LIBRARIES})
