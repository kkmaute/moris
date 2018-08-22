# -------------------------------------------------------------------------
# Eigen libraries ---------------------------------------------------------
# -------------------------------------------------------------------------
if(NOT EIGEN_FOUND_ONCE)
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
        NAMES Eigen3
        HINTS ${EIGEN_ENV_VARS})
message(${EIGEN3_INCLUDE_DIRS})
    message(STATUS "EIGEN_LIBRARIES: ${EIGEN_LIBRARIES}")
endif()

add_definitions("-DMORIS_USE_EIGEN")
include_directories(${EIGEN3_INCLUDE_DIRS})
message(${EIGEN3_INCLUDE_DIRS})
set(MORIS_ARMADILLO_EIGEN_LIBS "Eigen3::Eigen")
