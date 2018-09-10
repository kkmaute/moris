# -------------------------------------------------------------------------
# Linear algebra library fixing -------------------------------------------
# -------------------------------------------------------------------------

# Replace LAPACK libraries with ACML; MKL needs to stay for Paradiso (Trilinos)
if(MORIS_USE_ACML)
    string(REGEX REPLACE
        "[^;]*/lib64/liblapack\\.so;[^;]*/lib64/libblas\\.so"
        "$ENV{ACML_DIR}/lib/libacml.so"
        MORIS_T_LIBS
        "${MORIS_T_LIBS}" )
    string(REGEX REPLACE
        "[^;]*/lib/intel64(_lin)?/libmkl_intel_lp64\\.so;[^;]*/lib/intel64(_lin)?/libmkl_sequential\\.so;[^;]*/lib/intel64(_lin)?/libmkl_core\\.so"
        "$ENV{ACML_DIR}/lib/libacml.so"
        MORIS_T_LIBS
        "${MORIS_T_LIBS}" )
#     string(REGEX REPLACE 
#         "[^;]*/acml/gfortran64/lib/libacml\.so"
#         ""
#         MORIS_T_LIBS 
#         "${MORIS_T_LIBS}" )
endif()

# Replace ACML libraries with LAPACK; MKL needs to stay for Paradiso (Trilinos)
if(MORIS_USE_LAPACK)
    string(REGEX REPLACE 
        "[^;]*/acml/gfortran64/lib/libacml\\.so"
        "$ENV{LAPACK_DIR}/lib64/liblapack.so;$ENV{LAPACK_DIR}/lib64/libblas.so"
        MORIS_T_LIBS 
        "${MORIS_T_LIBS}" )
    string(REGEX REPLACE
        "[^;]*/lib/intel64(_lin)?/libmkl_intel_lp64\\.so;[^;]*/lib/intel64(_lin)?/libmkl_sequential\\.so;[^;]*/lib/intel64(_lin)?/libmkl_core\\.so"
        "$ENV{LAPACK_DIR}/lib64/liblapack.so;$ENV{LAPACK_DIR}/lib64/libblas.so"
        MORIS_T_LIBS
        "${MORIS_T_LIBS}" )
endif()

# Replace ACML and LAPACK libraries with MKL
if(MORIS_USE_MKL)
    string(REGEX REPLACE 
        "[^;]*/acml/gfortran64/lib/libacml\\.so"
        "$ENV{MKL_DIR}/lib/intel64_lin/libmkl_intel_lp64.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_core.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_sequential.so;/usr/lib64/libpthread.so"
        MORIS_T_LIBS 
        "${MORIS_T_LIBS}" )
    string(REGEX REPLACE 
        "[^;]*/lib64/liblapack\\.so;[^;]*/lib64/libblas\\.so"
        "$ENV{MKL_DIR}/lib/intel64_lin/libmkl_intel_lp64.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_core.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_sequential.so;/usr/lib64/libpthread.so"
        MORIS_T_LIBS
        "${MORIS_T_LIBS}" )
endif()
