# -------------------------------------------------------------------------
# LBFGSB libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(LBFGSB)

message(STATUS "LBFGSB_LIBRARIES: ${LBFGSB_LIBRARIES}")

list(APPEND MORIS_LDFLAGS "${LBFGSB_LIBRARY_DIRS}")
list(APPEND MORIS_LDLIBS "${LBFGSB_LIBRARIES}")
