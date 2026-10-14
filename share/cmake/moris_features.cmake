# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

if(MORIS_USE_ACML)
    if(MORIS_USE_MKL OR MORIS_USE_OPENBLAS OR MORIS_USE_LAPACK)
        message(FATAL_ERROR "Use either MKL or ACML, OPENBLAS, LAPACK")
    endif()
endif()

if(MORIS_USE_OPENBLAS)
    if(MORIS_USE_MKL OR MORIS_USE_ACML OR MORIS_USE_LAPACK)
        message(FATAL_ERROR "Use either MKL or ACML, OPENBLAS, LAPACK")
    endif()
endif()

if(MORIS_USE_LAPACK)
    if(MORIS_USE_MKL OR MORIS_USE_ACML OR MORIS_USE_OPENBLAS)
        message(FATAL_ERROR "Use either MKL or ACML, OPENBLAS, LAPACK")
    endif()
endif()

if(MORIS_USE_ACML)
    set(ACML_LAPACK_MKL_OPENBLAS "acml")
elseif(MORIS_USE_MKL)
    set(ACML_LAPACK_MKL_OPENBLAS "mkl")
elseif(MORIS_USE_LAPACK)
    set(ACML_LAPACK_MKL_OPENBLAS "lapack")
elseif(MORIS_USE_OPENBLAS)
    set(ACML_LAPACK_MKL_OPENBLAS "openblas")
else()
    message(
        FATAL_ERROR
            "MORIS only supports the AMDBLAS (amdblis/amdlibflame), OPENBLAS, LAPACK, and MKL linear algebra libraries."
    )
endif()

# Armadillo / Eigen

if(MORIS_USE_ARMA AND MORIS_USE_EIGEN)
    message(FATAL_ERROR "Can only set MORIS_USE_EIGEN or MORIS_USE_ARMA; not both")
endif()

set(ARMADILLO_EIGEN)

if(MORIS_USE_ARMA)
    list(APPEND ARMADILLO_EIGEN "armadillo")
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_ARMA")

    if(NOT MORIS_HAVE_DEBUG)
        list(APPEND MORIS_DEFINITIONS "-DARMA_NO_DEBUG")
    endif()

endif()

# Eigen3

if(MORIS_USE_EIGEN)
    list(APPEND ARMADILLO_EIGEN "eigen")
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_EIGEN")
endif()

# These switches affect public headers and must also reach installed consumers.
foreach(
    feature IN
    ITEMS MORIS_USE_32BIT
          MORIS_HAVE_PETSC
          MORIS_HAVE_SLEPC
          MORIS_USE_PARDISO
          MORIS_USE_MUMPS
          MORIS_HAVE_GCMMA
          MORIS_HAVE_SNOPT
          MORIS_HAVE_LBFGS
          MORIS_HAVE_DEBUG
          MORIS_PERFORM_CHECK
          MORIS_HAVE_PARALLEL
)
    if(${feature})
        list(APPEND MORIS_DEFINITIONS "${feature}")
    endif()
endforeach()

if(MORIS_USE_XTK)
    list(APPEND MORIS_DEFINITIONS USE_XTK XTK_USE_MORIS)
endif()
if(MORIS_USE_MATRIX_FILL)
    list(APPEND MORIS_DEFINITIONS MATRIX_FILL)
endif()
if(MORIS_USE_CHECK_MEMORY)
    list(APPEND MORIS_DEFINITIONS CHECK_MEMORY)
endif()
list(APPEND MORIS_DEFINITIONS F77ADD_)
