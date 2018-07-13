# -------------------------------------------------------------------------
# GCMMA libraries --------------------------------------------------------
# -------------------------------------------------------------------------

find_package(GCMMA)

message(STATUS "GCMMA_LIBRARIES: ${GCMMA_LIBRARIES}")

list(APPEND MORIS_INCDIRS "${GCMMA_INCLUDE_DIRS}")
list(APPEND MORIS_LDFLAGS "${GCMMA_LIBRARY_DIRS}")
list(APPEND MORIS_LDLIBS "${GCMMA_LIBRARIES}")
