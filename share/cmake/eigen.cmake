# -------------------------------------------------------------------------
# Eigen libraries ---------------------------------------------------------
# -------------------------------------------------------------------------
set(EIGEN_ENV_VARS
    $ENV{EIGENDIR}
    $ENV{EIGEN_DIR}
    $ENV{Eigen_DIR}
    $ENV{Eigen3_DIR}
    $ENV{EIGEN_ROOT}
    $ENV{Eigen_ROOT}
    $ENV{Eigen3_ROOT}
    $ENV{EIGEN_PATH}
    $ENV{Eigen_PATH}
    $ENV{Eigen3_PATH})

find_package(Eigen
    REQUIRED
    NAMES Eigen Eigen3
    HINTS ${EIGEN_ENV_VARS})

message(STATUS "EIGEN_LIBRARIES: ${EIGEN_LIBRARIES}")

# list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_EIGEN")
# list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_EIGEN")
# list(APPEND MORIS_INCDIRS ${EIGEN_DIRS})
add_definitions("-DMORIS_USE_EIGEN")
include_directories(/home/maute/tpls/eigen-3.3.3/include)

set(MORIS_ARMADILLO_EIGEN_LIBS "Eigen3::Eigen")
