# -----------------------------------------------------------------------------
# ACML libraries --------------------------------------------------------------
# -----------------------------------------------------------------------------

find_package(ACML)

message(STATUS "ACML_LIBRARIES: ${ACML_LIBRARIES}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_ACML")
list(APPEND MORIS_LDLIBS ${ACML_LIBRARIES})
