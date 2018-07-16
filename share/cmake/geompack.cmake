# -------------------------------------------------------------------------
# GEOMPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(GEOMPACK)

message(STATUS "GEOMPACK_LIBRARIES: ${GEOMPACK_LIBRARIES}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_GEOMPACK")
list(APPEND MORIS_LDLIBS ${GEOMPACK_LIBRARIES})
