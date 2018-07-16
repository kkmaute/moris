# -------------------------------------------------------------------------
# ARPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(ARPACK)

message(STATUS "ARPACK_LIBRARIES: ${ARPACK_LIBRARIES}")

list(APPEND MORIS_LDLIBS ${ARPACK_LIBRARIES})
