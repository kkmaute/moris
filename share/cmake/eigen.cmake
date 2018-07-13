# -------------------------------------------------------------------------
# Eigen libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

find_package(Eigen)

message(STATUS "EIGEN_LIBRARIES: ${EIGEN_LIBRARIES}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_EIGEN")
list(APPEND MORIS_INCDIRS ${EIGEN_DIRS})
