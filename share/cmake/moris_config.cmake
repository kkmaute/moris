list( APPEND MORIS_INCDIRS "${SRC}" )
list( APPEND MORIS_INCDIRS "snippets" )
list( APPEND MORIS_INCDIRS "${INC}" )
list( APPEND MORIS_INCDIRS "/usr/include" )

# -----------------------------------------------------------------------------
# Options ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

# MPI support.
option( MORIS_HAVE_PARALLEL
    "Build MORIS with support for MPI communication." ON )

# MPI package.
set(MORIS_MPI_LIBS "MPICH" "OPENMPI")
set(MORIS_USE_MPI "OPENMPI" CACHE STRING "Set MPI package: ${MORIS_MPI_LIBS}")

option( MORIS_USE_ARMA
    "Use the Armadillo library as the dense linear algebra package." ON )

option( MORIS_USE_EIGEN
    "Use the Eigen3 library as the dense linear algebra package." OFF )

option( HAVE_GEOMPACK
    "Build MORIS with support for Geompack." OFF )

option( MORIS_HAVE_GEOMPACK
    "Use Geompack3 library to perform geometry computations." OFF )

option( MORIS_USE_32BIT
    "Use 32bit numbers for identifiers and sizes." ON )

option( MORIS_HAVE_DEBUG
    "Use the debug version of the code." OFF )
    
option( USE_XTK "Have XTK Library" ON)	

option( MORIS_HAVE_SYMBOLIC
    "Use symbolic information in executable." OFF )

option( MORIS_USE_ACML
    "Use the AMD ACML library as the linear algebra package." ON )

option( MORIS_USE_MKL
    "Use the AMD ACML library as the linear algebra package." OFF )

option( MORIS_USE_BLAS
    "Use the BLAS and LAPACK as the linear algebra packages." OFF )

option( HAVE_GCMMA
    "Build MORIS with support for GCMMA." ON )

option( HAVE_LBFGS
    "Build MORIS with support for LBFGS." ON )

option( HAVE_ARPACK
    "Build MORIS with support for LBFGS." ON )

option( HAVE_SNOPT
    "Build MORIS with support for SNOPT." ON )
    
option( HAVE_PERF_MAT
    "Run matrix performance tests." OFF )
    
option( HAVE_PERF_SP_MAT
    "Run sparse matrix performance tests." OFF )
    
option( HAVE_PERF_LIN_SOLVE
    "Run linear solver performance tests." OFF )

option( MORIS_USE_TESTS
    "Compile unit tests." ON )
    
option( MORIS_USE_REGRESSION
    "Compile regression tests." ON )

option( MORIS_PERFORM_CHECK
    "Use performance check." ON )
    
option( USE_TACC "USE TACC" OFF)
    
# -----------------------------------------------------------------------------

# Set default C compiler.
if ( NOT MORIS_C_COMPILER )
    set( MORIS_C_COMPILER "mpicc" )
endif()
set( MORIS_C_LINK_EXECUTABLE ${MORIS_C_COMPILER} )

# Set default C++ compiler
if ( NOT MORIS_CXX_COMPILER )
    if ( USE_TACC ) 
        set( MORIS_CXX_COMPILER "mpicxx" )
    else()
        set( MORIS_CXX_COMPILER "mpic++" )
    endif()
endif()
set( MORIS_CXX_LINK_EXECUTABLE ${MORIS_CXX_COMPILER} )

message( STATUS "MORIS recognized the C compiler type MORIS_C_COMPILER=${MORIS_C_COMPILER}." )
message( STATUS "MORIS recognized the C++ compiler type MORIS_CXX_COMPILER=${MORIS_CXX_COMPILER}." )

# -----------------------------------------------------------------------------
# XTK -------------------------------------------------------------------------
# -----------------------------------------------------------------------------    
if ( USE_XTK )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/include" )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/src" )
    list( APPEND MORIS_DEFINITIONS "-DUSE_XTK" )
    list( APPEND MORIS_DEFINITIONS "-DXTK_USE_MORIS")
endif()

    
# -----------------------------------------------------------------------------
# PETSc -----------------------------------------------------------------------
# -----------------------------------------------------------------------------

include(${SHARE}/cmake/PETSc.cmake)

if ( USE_TACC )
   list( APPEND MORIS_INCDIRS "$ENV{$PETSC_DIR}/include" )
   list( APPEND MORIS_INCDIRS "$ENV{$PETSC_DIR}/$ENV{$PETSC_ARCH}/include" )
else()
   list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/petsc/gcc-openmpi/include" )
endif()

# -----------------------------------------------------------------------------
# ViennaCL  -------------------------------------------------------------------
# -----------------------------------------------------------------------------

list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/ViennaCL/include" )

# -----------------------------------------------------------------------------
# Trilinos---------------------------------------------------------------------
# -----------------------------------------------------------------------------

include(${SHARE}/cmake/trilinos.cmake)

# -----------------------------------------------------------------------------
# Boost -----------------------------------------------------------------------
# -----------------------------------------------------------------------------

include(${SHARE}/cmake/boost.cmake)

# -----------------------------------------------------------------------------
# ACML/MKL/LAPACK libraries ---------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( MORIS_USE_ACML AND MORIS_USE_MKL )
    message( FATAL_ERROR "Use either MKL or ACML" )
endif()

if ( MORIS_USE_ACML AND MORIS_USE_BLAS )
    message( FATAL_ERROR "Use either BLAS or ACML" )
endif()

if ( MORIS_USE_BLAS AND MORIS_USE_MKL )
    message( FATAL_ERROR "Use either MKL or BLAS" )
endif()


if ( MORIS_USE_ACML )
    set( ACML_LIBRARIES "$ENV{HOME}/tpls/acml/gfortran64/lib/libacml.so" )   
    include(${SHARE}/cmake/acml.cmake)
elseif( MORIS_USE_MKL )
    list( APPEND MORIS_DEFINITIONS "-DMORIS_USE_MKL" )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/mkl/include" )
    list( APPEND MORIS_LDLIBS "-L$ENV{HOME}/tpls/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread" )
elseif( MORIS_USE_BLAS )
    list( APPEND MORIS_DEFINITIONS "-DMORIS_USE_BLAS" )
    list( APPEND MORIS_LDLIBS "-L$ENV{HOME}/tpls/lapack/lib64 -llapack -lblas" )
else()
    message( FATAL_ERROR "MORIS only supports the BLAS and ACML linear algebra libraries." )
endif()

if ( MORIS_USE_32BIT )
    list( APPEND MORIS_DEFINITIONS "-DMORIS_USE_32BIT" )
endif()

# Check for compiler flags.
include( CheckCXXCompilerFlag )

# Check if include file exists.
include( CheckIncludeFile )

# -----------------------------------------------------------------------------
# Compiler --------------------------------------------------------------------
# -----------------------------------------------------------------------------

set( CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition" )

if ( USE_TACC )
    list( APPEND MORIS_LDLIBS "-lgfortran -lrt -ldl" )
else()
    list( APPEND MORIS_LDLIBS "-lmpi_mpifh -lgfortran -lrt -ldl -lssl -lcrypto" )
endif()

list( APPEND MORIS_DEFINITIONS "-DF77ADD_")

# -g Request that the compiler and linker generate 
# and retain symbol information in the executable itself,
# for debugging purposes.
if ( MORIS_HAVE_SYMBOLIC )
    check_cxx_compiler_flag( -g HAVE_DEBUG )
    if ( HAVE_DEBUG )
        set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g" )
    endif()
endif()

# Debug flags.
if ( MORIS_HAVE_DEBUG )
    list( APPEND MORIS_DEFINITIONS "-DDEBUG" )
else()
    list( APPEND MORIS_DEFINITIONS "-DNDEBUG" )
endif()

# Performace test flags
if ( HAVE_PERF_MAT )
    list( APPEND MORIS_DEFINITIONS "-DPERF_MAT" )
endif()

if ( HAVE_PERF_SP_MAT )
    list( APPEND MORIS_DEFINITIONS "-DPERF_SP_MAT" )
endif()

if ( HAVE_PERF_LIN_SOLVE )
    list( APPEND MORIS_DEFINITIONS "-DPERF_LIN_SOLVE" )
endif()

# Performace check flag for logger
if(MORIS_PERFORM_CHECK)
        list(APPEND MORIS_DEFINITIONS "-DMORIS_PERFORM_CHECK")
endif()

# Set optimization type.
# Compile with -O compiler optimizations.
if ( NOT MORIS_HAVE_DEBUG )
    check_cxx_compiler_flag( -O3 HAVE_O3_OPTIMISATION )
    if ( HAVE_O3_OPTIMISATION )
        set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -O3 -m64" )
    endif()
else()
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g -ggdb -m64" )
endif()

# Add some strict compiler checks.
# -pedantic-errors     Make all pedantic warnings into errors. -pedantic issues 
#                      all the warnings demanded by strict ISO C and ISO C++.
# -Wall                This enables all the warnings about constructions 
#                      that some users consider questionable, and
#                      that are easy to avoid (or modify to prevent the warning),
#                      even in conjunction with macros.
#                      However, it does not enable all warnings available.
# -Werror              Make all warnings into errors.
# -Wconversion         Give warning if conversion between data types occurs.
# -Wno-long-long       Do not issue a warning if a long long variable type is used.
# -fno-strict-aliasing Do not enforce strict aliasing.
#                      Strict aliasing means that pointer arguments in a function are assumed to not alias.
#                      For example, the following code would not compile: foo * a; bar * b; b = ( foo * ) a;
#                      Because the pointers point to fundamentally different types.
check_cxx_compiler_flag( "-Wall -Werror -Wno-long-long -pedantic-errors" HAVE_PEDANTIC )
if ( HAVE_PEDANTIC )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wall -Werror -Wno-long-long -pedantic-errors")
endif()

# Check for C++11 or C++0x support.
check_cxx_compiler_flag( -std=c++11 HAVE_STD_CPP11 )
if ( HAVE_STD_CPP11 )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++11" )
else()
    check_cxx_compiler_flag( -std=c++0x HAVE_STD_CPP0x )
    if ( HAVE_STD_CPP0x )
        set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++0x" )
    endif()
endif()

# Build 64-bit binaries.
check_cxx_compiler_flag( -m64 HAVE_M64 )
if ( HAVE_M64 )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -m64" )
endif()

# Add UTF-8 support for identifier names.
check_cxx_compiler_flag( -fextended-identifiers HAVE_FEXT_IDENTIFIERS )
if ( HAVE_FEXT_IDENTIFIERS )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -fextended-identifiers" )
endif()

# Add support to use shared libraries as input files.
# -rdynamic Pass the flag -export-dynamic to the ELF linker, 
#           on targets that support it. This instructs the linker
#           to add all symbols, not only used ones,
#           to the dynamic symbol table. This option is needed for
#           some uses of dlopen or to allow obtaining backtraces
#           from within a program.
# -fPIC     Generate position-independent code (PIC) suitable
#           for use in a shared library, if supported for the target machine.
check_cxx_compiler_flag( -rdynamic HAVE_RDYNAMIC )
if ( HAVE_RDYNAMIC )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -rdynamic" )
endif()

check_cxx_compiler_flag( -fPIC HAVE_FPIC )
if ( HAVE_FPIC )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -fPIC" )
endif()

# -----------------------------------------------------------------------------
# MPI libraries and includes --------------------------------------------------
# -----------------------------------------------------------------------------

include(${SHARE}/cmake/mpi.cmake)

# -----------------------------------------------------------------------------
# Geompack3 -------------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( HAVE_GEOMPACK )
    list( APPEND MORIS_DEFINITIONS "-DHAVE_GEOMPACK" )
    list( APPEND MORIS_LDLIBS "-lgeompack3" )
else()
#   message( FATAL_ERROR "Geompack3 is required to run project ${PROJECT_NAME}" )
endif()


# -----------------------------------------------------------------------------
# Armadillo -------------------------------------------------------------------
# -----------------------------------------------------------------------------
if ( MORIS_USE_ARMA AND MORIS_USE_EIGEN )
	    message( FATAL_ERROR "Can only set MORIS_USE_EIGEN or MORIS_USE_ARMA; but not both" )
endif()



if ( MORIS_USE_ARMA )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/armadillo/${INC}" )
    list( APPEND MORIS_LDFLAGS "$ENV{HOME}/tpls/armadillo/lib" )
    list( APPEND MORIS_LDLIBS "-larmadillo" )
    list( APPEND MORIS_DEFINITIONS "-DMORIS_USE_ARMA" )
endif()

# -----------------------------------------------------------------------------
# Eigen3 ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( MORIS_USE_EIGEN )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/eigen/include/Eigen" )
    list( APPEND MORIS_DEFINITIONS "-DMORIS_USE_EIGEN" )
endif()

# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

include(${SHARE}/cmake/superlu.cmake)

if ( USE_TACC )
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/SuperLU-DIST/intel-parallel/include" )
    list( APPEND MORIS_LDLIBS "$ENV{HOME}/tpls/SuperLU-DIST/intel-parallel/lib64/libsuperlu_dist.a" )
else()
    list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/SuperLU-DIST-5.1/gcc-openmpi/include" )
    list( APPEND MORIS_LDLIBS "$ENV{HOME}/tpls/SuperLU-DIST-5.1/gcc-openmpi/lib/libsuperlu_dist.a" )
endif()

# -----------------------------------------------------------------------------
# GCMMA ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

list( APPEND MORIS_LDFLAGS "$ENV{HOME}/tpls/gcmma/lib" )
list( APPEND MORIS_LDLIBS "-lgcmma" )
list( APPEND MORIS_INCDIRS "$ENV{HOME}/tpls/gcmma/include" )

# -----------------------------------------------------------------------------
# SQP -------------------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( HAVE_SNOPT )
    list( APPEND MORIS_LDFLAGS "$ENV{HOME}/tpls/Snopt/lib" )
    list( APPEND MORIS_LDLIBS "-lsnopt" )
endif()

# -----------------------------------------------------------------------------
# LBFGS -----------------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( HAVE_LBFGS )
    list( APPEND MORIS_LDFLAGS "$ENV{HOME}/tpls/Lbfgsb/lib" )
    list( APPEND MORIS_LDLIBS "-llbfgsb" )
endif()


# -----------------------------------------------------------------------------
# ARPACL ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

if ( HAVE_ARPACK )
    list( APPEND MORIS_LDFLAGS "$ENV{HOME}/tpls/ARPACK/lib" )
    list( APPEND MORIS_LDLIBS "-larpack" )
endif()

# -----------------------------------------------------------------------------
# Replace all acml library calls with calls to lapack and blas ----------------
# This should be the last entry in this file ----------------------------------
# -----------------------------------------------------------------------------

if ( MORIS_USE_BLAS )
    STRING( REGEX REPLACE "$ENV{HOME}/tpls/acml/gfortran64/lib/libacml.so" "$ENV{HOME}/tpls/lapack/lib64/liblapack.so;$ENV{HOME}/tpls/lapack/lib64/libblas.so" MORIS_LDLIBS "${MORIS_LDLIBS}" )
    STRING( REGEX REPLACE "/home/maute/tpls/acml/gfortran64/lib/libacml.so" "$ENV{HOME}/tpls/lapack/lib64/liblapack.so;$ENV{HOME}/tpls/lapack/lib64/libblas.so" MORIS_LDLIBS "${MORIS_LDLIBS}" )
endif()

if ( MORIS_USE_MKL )
    STRING( REGEX REPLACE "$ENV{HOME}/tpls/acml/gfortran64/lib/libacml.so" "$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_intel_lp64.so;$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_core.so;$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_sequential.so;/usr/lib64/libpthread.so" MORIS_LDLIBS "${MORIS_LDLIBS}" )
    STRING( REGEX REPLACE "/home/maute/tpls/acml/gfortran64/lib/libacml.so" "$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_intel_lp64.so;$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_core.so;$ENV{HOME}/tpls/mkl/lib/intel64/libmkl_sequential.so;/usr/lib64/libpthread.so" MORIS_LDLIBS "${MORIS_LDLIBS}" )
endif()
