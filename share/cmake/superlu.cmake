# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

find_package(SuperLU)

message(STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
message(STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")

# list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_SUPERLU")
list(APPEND MORIS_INCDIRS ${SUPERLU_INCLUDES})
list(APPEND MORIS_LDLIBS ${SUPERLU_LIBRARIES})

# -------------------------------------------------------------------------
# SuperLU_DIST

find_package(SuperLU_DIST)

message(STATUS "SUPERLU_DIST_LIBRARIES: ${SUPERLU_DIST_LIBRARIES}")

list(APPEND MORIS_LDLIBS ${SUPERLU_DIST_LIBRARIES})
