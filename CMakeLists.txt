# -------------------------------------------------------------------------
# MORIS CMakeLists.txt ----------------------------------------------------
# -------------------------------------------------------------------------

# NOTES:
# The command 'add_library(${MORIS} STATIC ${SOURCES})' automatically
# adds the 'lib' prefix to the library and the 'a' file extension.
# The extension 'a' is added for a static library.
# Result is 'Linking CXX static library libmoris.a'.

# Require CMake 3.1
cmake_minimum_required(VERSION 3.0.2)
message("${CMAKE_BINARY_DIR}")
# Guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt from this directory.")
endif()

# Project name and version
project(MORIS VERSION 1.0.0)

set(MORIS_VERSION_RELEASE 0)

if(NOT MORIS_VERSION_RELEASE)
    set(MORIS_VERSION "${MORIS_VERSION}.dbg")
endif()

# -------------------------------------------------------------------------
## MORIS Directories

# Top directory
set(MORIS "moris")

# Project directories
set(BIN "bin")
set(BUILD "build")
set(CMAKE "cmake")
set(INC "include")
set(LIB "lib")
set(SHARE "share")
set(SRC "projects")
set(TMP "tmp")
set(TEST "test")
set(REGRESSION "regression")
set(UNIT "unit")
set(MAIN "mains")

# Source directories
set(ALG "ALG")
set(ASR "ASR")
set(CHR "CHR")
set(COM "COM")
set(CON "CON")
set(COR "COR")
set(DLA "DLA")
set(EXC "EXC")
set(FEM "FEM")
set(GEN "GEN")
set(HIA "Hierarchical")
set(HMR "HMR")
set(INT "INT")
set(IOS "IOS")
#set(LNA "LNA")
set(MDL "MDL")
set(LINALG "LINALG")
set(MOD "MOD")
set(MRS "MRS")
set(MSI "MSI")
set(MTK "MTK")
set(NLA "NLA")
set(OPT "OPT")
set(SDF "SDF")
set(STK "STK")
set(TIN "TIN")
set(TOL "TOL")
set(XTK "XTK")

# Common paths
set(MORIS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MORIS_CMAKE_DIR ${MORIS_DIR}/${SHARE}/${CMAKE})
set(MORIS_DEPENDS_DIR ${MORIS_CMAKE_DIR}/depends)
set(MORIS_PACKAGE_DIR ${MORIS_DIR}/projects)

# Build directories
set(CONFIG_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE})
set(LIB_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(EXE_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
set(HEADER_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${INC})
set(TEST_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${TEST})
set(REGRESSION_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/${REGRESSION})

# Output path for all static and shared libraries.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Installation directories
set(CONFIG_INSTALL_DIR ${SHARE})
set(LIB_INSTALL_DIR ${LIB})
set(EXE_INSTALL_DIR ${BIN})
set(HEADER_INSTALL_DIR ${INC})

# -------------------------------------------------------------------------
## CMake Paths

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/${SHARE}/${CMAKE}/modules")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")

# -------------------------------------------------------------------------
## MORIS Important Include Directories

list(APPEND MORIS_INCDIRS "${SRC}" "snippets" "${INC}" "/usr/include")

# -------------------------------------------------------------------------
## Options

# MPI support.
option(MORIS_HAVE_PARALLEL
    "Build MORIS with support for MPI communication." ON )

# MPI package.
set(MORIS_MPI_LIBS "MPICH" "OPENMPI")
set(MORIS_USE_MPI "OPENMPI" CACHE STRING "Set MPI package: ${MORIS_MPI_LIBS}")

option(MORIS_USE_ARMA
    "Use the Armadillo library as the dense linear algebra package." ON )

option(MORIS_USE_EIGEN
    "Use the Eigen3 library as the dense linear algebra package." OFF )

# option(HAVE_GEOMPACK
#     "Build MORIS with support for Geompack." OFF )

# option(MORIS_HAVE_GEOMPACK
#     "Use Geompack3 library to perform geometry computations." OFF )

option(MORIS_USE_32BIT
    "Use 32bit numbers for identifiers and sizes." ON )

option(MORIS_HAVE_DEBUG
    "Use the debug version of the code." OFF )
    
option(USE_XTK "Have XTK Library" ON)	

option(MORIS_HAVE_SYMBOLIC
    "Use symbolic information in executable." OFF )

option(MORIS_USE_ACML
    "Use the AMD ACML library as the linear algebra package." OFF )

option(MORIS_USE_MKL
    "Use the Inel Math kernel library as the linear algebra package." OFF )

option(MORIS_USE_LAPACK
    "Use the LAPACK and LAPACK as the linear algebra packages." ON )

option(HAVE_GCMMA
    "Build MORIS with support for GCMMA." ON )

# option(HAVE_LBFGSB
#     "Build MORIS with support for LBFGSB." ON )

# option(HAVE_ARPACK
#     "Build MORIS with support for LBFGSB." ON )

option(HAVE_SNOPT
    "Build MORIS with support for SNOPT." ON )
    
option(HAVE_PERF_MAT
    "Run matrix performance tests." OFF )
    
option(HAVE_PERF_SP_MAT
    "Run sparse matrix performance tests." OFF )
    
option(HAVE_PERF_LIN_SOLVE
    "Run linear solver performance tests." OFF )

option(MORIS_USE_TESTS
    "Compile unit tests." ON )
    
# option(MORIS_USE_REGRESSION
#     "Compile regression tests." ON )

option(MORIS_PERFORM_CHECK
    "Use performance check." ON )
    
option(USE_TACC "USE TACC" OFF)

# -------------------------------------------------------------------------
# Package Options

option(BUILD_MAIN "Build main executable." OFF)

option(BUILD_ALL "Build all executables." OFF)

option(BUILD_ALG "Build the algorithms executable." OFF)

option(BUILD_ASR "Build the assert executable." OFF)

option(BUILD_CHR "Build the chronos executable." OFF)

option(BUILD_COM "Build the communication executable." OFF)

option(BUILD_CON "Build the containers executable." OFF)

option(BUILD_DLA "Build the distributed linear algebra executable." OFF)

option(BUILD_EXC "Build the exceptions executable." OFF)

option(BUILD_FEM "Build the FEM executable." OFF)

option(BUILD_GEN "Build the geometry engine executable." OFF)

option(BUILD_HMR "Build the HMR executable." OFF)

option(BUILD_INT "Build the integration and interpolation executable." OFF)

option(BUILD_IOS "Build the IOS executable." OFF)

#option(BUILD_LNA "Build the linear algebra executable." OFF)

option(BUILD_LINALG "Build the linear algebra executable." OFF)

option(BUILD_MOD "Build the model executable." OFF)

option(BUILD_MSI "Build the model solver interface executable." OFF)

option(BUILD_MDL "Build the model executable." OFF)

option(BUILD_MTK "Build the MTK executable." OFF)

option(BUILD_NLA "Build the non-linear algebra executable." OFF)

option(BUILD_OPT "Build the optimization executable." OFF)

option(BUILD_STK "Build the STK executable." OFF)

option(BUILD_TIN "Build the TIN executable." OFF)

option(BUILD_TOL "Build the tools executable." OFF)

option(BUILD_TUTORIALS "Build the tutorial executables." ON)

option(BUILD_XTK "Build the XTK executable." OFF)


# -------------------------------------------------------------------------
# Compiler ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Set default C compiler.
if (NOT MORIS_C_COMPILER )
    set(MORIS_C_COMPILER "mpicc" )
endif()
set(MORIS_C_LINK_EXECUTABLE ${MORIS_C_COMPILER} )

# Set default C++ compiler
if (NOT MORIS_CXX_COMPILER )
    if (USE_TACC ) 
        set(MORIS_CXX_COMPILER "mpicxx" )
    else()
        set(MORIS_CXX_COMPILER "mpic++")
    endif()
endif()
set(MORIS_CXX_LINK_EXECUTABLE ${MORIS_CXX_COMPILER})

message(STATUS "MORIS recognized the C compiler type MORIS_C_COMPILER=${MORIS_C_COMPILER}.")
message(STATUS "MORIS recognized the C++ compiler type MORIS_CXX_COMPILER=${MORIS_CXX_COMPILER}.")

# Check for compiler flags.
include(CheckCXXCompilerFlag)

set(CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition")

if (USE_TACC)
    list(APPEND MORIS_BASE_LIBS "-lgfortran -lrt -ldl")
else()
#     list(APPEND MORIS_BASE_LIBS "-lmpi_mpifh -lgfortran -lrt -ldl -lssl -lcrypto")
   list(APPEND MORIS_BASE_LIBS "-lmpi_mpifh -lgfortran -ldl -lssl -lcrypto")
endif()

list(APPEND MORIS_DEFINITIONS "-DF77ADD_")

# -g Request that the compiler and linker generate 
# and retain symbol information in the executable itself,
# for debugging purposes.
if (MORIS_HAVE_SYMBOLIC)
    check_cxx_compiler_flag(-g HAVE_DEBUG)
    if (HAVE_DEBUG)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g")
    endif()
endif()

# Debug flags.
if (MORIS_HAVE_DEBUG)
    list(APPEND MORIS_DEFINITIONS "-DDEBUG")
else()
    list(APPEND MORIS_DEFINITIONS "-DNDEBUG")
endif()

# Performace test flags
if (HAVE_PERF_MAT)
    list(APPEND MORIS_DEFINITIONS "-DPERF_MAT")
endif()

if (HAVE_PERF_SP_MAT)
    list(APPEND MORIS_DEFINITIONS "-DPERF_SP_MAT")
endif()

if (HAVE_PERF_LIN_SOLVE)
    list(APPEND MORIS_DEFINITIONS "-DPERF_LIN_SOLVE")
endif()

# Performace check flag for logger
if(MORIS_PERFORM_CHECK)
        list(APPEND MORIS_DEFINITIONS "-DMORIS_PERFORM_CHECK")
endif()

# Set optimization type.
# Compile with -O compiler optimizations.
if (NOT MORIS_HAVE_DEBUG)
    check_cxx_compiler_flag(-O3 HAVE_O3_OPTIMISATION)
    if (HAVE_O3_OPTIMISATION)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -O3 -m64")
    endif()
else()
#> redundant? see stuff starting at line 210
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g -ggdb -m64")
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
# -Wno-error=maybe-uninitialized
#                      Do not issue an error if a maybe-uninitialized warning is thrown
# -fno-strict-aliasing Do not enforce strict aliasing.
#                      Strict aliasing means that pointer arguments in a function are assumed to not alias.
#                      For example, the following code would not compile: foo * a; bar * b; b = (foo *) a;
#                      Because the pointers point to fundamentally different types.
check_cxx_compiler_flag("-Wall -Werror -Wno-long-long -pedantic-errors -Wno-error=maybe-uninitialized" HAVE_PEDANTIC )
if (HAVE_PEDANTIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wall -Werror -Wno-long-long -pedantic-errors -Wno-error=maybe-uninitialized" )
endif()

# Check for C++11 or C++0x support.
check_cxx_compiler_flag(-std=c++11 HAVE_STD_CPP11)
if (HAVE_STD_CPP11)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++11")
else()
    check_cxx_compiler_flag(-std=c++0x HAVE_STD_CPP0x)
    if (HAVE_STD_CPP0x)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++0x")
    endif()
endif()

# Build 64-bit binaries.
check_cxx_compiler_flag(-m64 HAVE_M64)
if (HAVE_M64)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -m64")
endif()

# Add UTF-8 support for identifier names.
check_cxx_compiler_flag(-fextended-identifiers HAVE_FEXT_IDENTIFIERS)
if (HAVE_FEXT_IDENTIFIERS)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -fextended-identifiers")
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
check_cxx_compiler_flag(-rdynamic HAVE_RDYNAMIC)
if (HAVE_RDYNAMIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -rdynamic")
endif()

check_cxx_compiler_flag(-fPIC HAVE_FPIC)
if (HAVE_FPIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -fPIC")
endif()


# -------------------------------------------------------------------------
# GPERFTOOLS : Make sure that paths are set correctly.
# -------------------------------------------------------------------------
option(USE_GPERFTOOLS "Use Google Profiling tools (check paths in CMakeLists.txt)." OFF)
if (USE_GPERFTOOLS)
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -I/home/messe/opt/gperftools/include -DWITHGPERFTOOLS -g" )
    list(APPEND MORIS_BASE_LIBS "-L/home/messe/opt/gperftools/lib64/ -lprofiler -ltcmalloc" )
endif()

# -------------------------------------------------------------------------
# XTK --------------------------------------------------------------------------
# -------------------------------------------------------------------------   
if (USE_XTK)
    list(APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/include") #> ask if this is separate from the files in MORIS/src/xtk
    list(APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/src")
    list(APPEND MORIS_DEFINITIONS "-DUSE_XTK")
    list(APPEND MORIS_DEFINITIONS "-DXTK_USE_MORIS")
endif()


# -------------------------------------------------------------------------
# PETSc -------------------------------------------------------------------
# -------------------------------------------------------------------------

# include(${MORIS_CMAKE_DIR}/PETSc.cmake)

# if (USE_TACC)
# #    list(APPEND MORIS_INCDIRS "$ENV{$PETSC_DIR}/include") #> this is adding /include
# #    list(APPEND MORIS_INCDIRS "$ENV{$PETSC_DIR}/$ENV{$PETSC_ARCH}/include") #> same as above
# else()
#    list(APPEND MORIS_INCDIRS ${PETSC_INCLUDE_DIR})
# endif()


# -------------------------------------------------------------------------
# ViennaCL  --------------------------------------------------------------------
# -------------------------------------------------------------------------

# find_package(ViennaCL)
# list(APPEND MORIS_INCDIRS ${VIENNACL_INCLUDE_DIRS})
# include(${MORIS_CMAKE_DIR}/viennacl.cmake)


# -------------------------------------------------------------------------
# Trilinos----------------------------------------------------------------------
# -------------------------------------------------------------------------

# include(${MORIS_CMAKE_DIR}/trilinos.cmake)


# -------------------------------------------------------------------------
# Boost ------------------------------------------------------------------------
# -------------------------------------------------------------------------

# include(${MORIS_CMAKE_DIR}/boost.cmake)


# -------------------------------------------------------------------------
# ACML/MKL/LAPACK libraries -----------------------------------------------
# -------------------------------------------------------------------------

if (MORIS_USE_ACML AND MORIS_USE_MKL)
    message(FATAL_ERROR "Use either MKL or ACML")
endif()

if (MORIS_USE_ACML AND MORIS_USE_LAPACK)
    message(FATAL_ERROR "Use either LAPACK or ACML")
endif()

if (MORIS_USE_LAPACK AND MORIS_USE_MKL)
    message(FATAL_ERROR "Use either MKL or LAPACK")
endif()

if (MORIS_USE_ACML)
    #include(${SHARE}/${CMAKE}/acml.cmake)
    set(ACML_LAPACK_MKL "acml")
elseif(MORIS_USE_MKL)
    #include(${SHARE}/${CMAKE}/mkl.cmake)
    set(ACML_LAPACK_MKL "mkl")
elseif(MORIS_USE_LAPACK)
    #include(${SHARE}/${CMAKE}/lapack.cmake)
    set(ACML_LAPACK_MKL "lapack")
else()
    message(FATAL_ERROR "MORIS only supports the ACML and LAPACK linear algebra libraries." )
endif()

if (MORIS_USE_32BIT)
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_32BIT")
endif()


# -------------------------------------------------------------------------
# MPI libraries and includes ----------------------------------------------
# -------------------------------------------------------------------------

#include(${MORIS_CMAKE_DIR}/mpi.cmake)


# -------------------------------------------------------------------------
# Geompack3 ---------------------------------------------------------------
# -------------------------------------------------------------------------

# if (HAVE_GEOMPACK)
#     include(${MORIS_CMAKE_DIR}/geompack.cmake)
# #     list(APPEND MORIS_DEFINITIONS "-DHAVE_GEOMPACK")
# #     list(APPEND MORIS_BASE_LIBS "-lgeompack3")
# endif()


# -------------------------------------------------------------------------
# Armadillo ---------------------------------------------------------------
# -------------------------------------------------------------------------

if (MORIS_USE_ARMA AND MORIS_USE_EIGEN)
    message(FATAL_ERROR "Can only set MORIS_USE_EIGEN or MORIS_USE_ARMA; not both" )
endif()
set(ARMADILLO_EIGEN)
if (MORIS_USE_ARMA)
    list(APPEND ARMADILLO_EIGEN "armadillo")
#     include(${MORIS_CMAKE_DIR}/armadillo.cmake)
#     find_package(Armadillo REQUIRED HINTS "$ENV{Armadillo_DIR}")
#     list(APPEND MORIS_INCDIRS "${ARMADILLO_INCLUDE_DIRS}")
#     list(APPEND MORIS_LDFLAGS "${ARMADILLO_LIBRARY_DIRS}")
#     list(APPEND MORIS_BASE_LIBS "-larmadillo")
#     list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_ARMA")
endif()


# -------------------------------------------------------------------------
# Eigen3 ------------------------------------------------------------------
# -------------------------------------------------------------------------

if (MORIS_USE_EIGEN)
    list(APPEND ARMADILLO_EIGEN "eigen")
#     list(APPEND MORIS_INCDIRS "$ENV{Eigen3_DIR}/include/Eigen") #> can't just do find_package becuase .../include/Eigen isn't in a variable and MORIS source seems written to need it
#     list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_EIGEN")
endif()


# -------------------------------------------------------------------------
# SuperLU libraries and includes ------------------------------------------
# -------------------------------------------------------------------------

#include(${MORIS_CMAKE_DIR}/superlu.cmake)

# list(APPEND MORIS_INCDIRS "$ENV{SuperLU_DIR}/include")
# 
# if(USE_TACC)
#     list(APPEND MORIS_BASE_LIBS "$ENV{SuperLU_DIST_DIR}/lib64/libsuperlu_dist.a")
# else()
#     list(APPEND MORIS_BASE_LIBS "$ENV{SuperLU_DIST_DIR}/lib/libsuperlu_dist.a")
# endif()


# -------------------------------------------------------------------------
# GCMMA -------------------------------------------------------------------
# -------------------------------------------------------------------------

# include(${MORIS_CMAKE_DIR}/gcmma.cmake)
# list(APPEND MORIS_LDFLAGS "$ENV{GCMMA_DIR}/lib")
# list(APPEND MORIS_BASE_LIBS "-lgcmma")
# list(APPEND MORIS_INCDIRS "$ENV{GCMMA_DIR}/include")


# -------------------------------------------------------------------------
# SQP ---------------------------------------------------------------------
# -------------------------------------------------------------------------

# Moris uses SNOPT for SQP

# if (HAVE_SNOPT)
#     include(${MORIS_CMAKE_DIR}/snopt.cmake)
# #     list(APPEND MORIS_LDFLAGS "$ENV{SNOPT_DIR}/lib") 
# #     list(APPEND MORIS_BASE_LIBS "-lsnopt")
# endif()


# -------------------------------------------------------------------------
# LBFGSB ------------------------------------------------------------------
# -------------------------------------------------------------------------
#> old was just called LBFGS, but we link the -B version so this must be specified

# if (HAVE_LBFGSB)
#     include(${MORIS_CMAKE_DIR}/lbfgsb.cmake)
# #     list(APPEND MORIS_LDFLAGS "$ENV{LBFGSB_DIR}/lib")
# #     list(APPEND MORIS_BASE_LIBS "-llbfgsb")
# endif()


# -------------------------------------------------------------------------
# ARPACK ------------------------------------------------------------------
# -------------------------------------------------------------------------

# if (HAVE_ARPACK)
    include(${MORIS_CMAKE_DIR}/arpack.cmake)
# #     list(APPEND MORIS_LDFLAGS "$ENV{ARPACK_DIR}/lib")
# #     list(APPEND MORIS_BASE_LIBS "-larpack")
# endif()


# # -------------------------------------------------------------------------
# # Linear algebra library fixing -------------------------------------------
# # -------------------------------------------------------------------------
# 
# # Replace LAPACK libraries with ACML; MKL needs to stay for Paradiso (Trilinos)
# if(MORIS_USE_ACML)
#     string(REGEX REPLACE
#         "[^;]*/lib64/liblapack\\.so;[^;]*/lib64/libblas\\.so"
#         "$ENV{ACML_DIR}/lib/libacml.so"
#         MORIS_BASE_LIBS
#         "${MORIS_BASE_LIBS}" )
#     string(REGEX REPLACE
#         "[^;]*/lib/intel64(_lin)?/libmkl_intel_lp64\\.so;[^;]*/lib/intel64(_lin)?/libmkl_sequential\\.so;[^;]*/lib/intel64(_lin)?/libmkl_core\\.so"
#         "$ENV{ACML_DIR}/lib/libacml.so"
#         MORIS_BASE_LIBS
#         "${MORIS_BASE_LIBS}" )
# #     string(REGEX REPLACE 
# #         "[^;]*/acml/gfortran64/lib/libacml\.so"
# #         ""
# #         MORIS_BASE_LIBS 
# #         "${MORIS_BASE_LIBS}" )
# endif()
# 
# # Replace ACML libraries with LAPACK; MKL needs to stay for Paradiso (Trilinos)
# if(MORIS_USE_LAPACK)
#     string(REGEX REPLACE 
#         "[^;]*/acml/gfortran64/lib/libacml\\.so"
#         "$ENV{LAPACK_DIR}/lib64/liblapack.so;$ENV{LAPACK_DIR}/lib64/libblas.so"
#         MORIS_BASE_LIBS 
#         "${MORIS_BASE_LIBS}" )
#     string(REGEX REPLACE
#         "[^;]*/lib/intel64(_lin)?/libmkl_intel_lp64\\.so;[^;]*/lib/intel64(_lin)?/libmkl_sequential\\.so;[^;]*/lib/intel64(_lin)?/libmkl_core\\.so"
#         "$ENV{LAPACK_DIR}/lib64/liblapack.so;$ENV{LAPACK_DIR}/lib64/libblas.so"
#         MORIS_BASE_LIBS
#         "${MORIS_BASE_LIBS}" )
# endif()
# 
# # Replace ACML and LAPACK libraries with MKL
# if(MORIS_USE_MKL)
#     string(REGEX REPLACE 
#         "[^;]*/acml/gfortran64/lib/libacml\\.so"
#         "$ENV{MKL_DIR}/lib/intel64_lin/libmkl_intel_lp64.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_core.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_sequential.so;/usr/lib64/libpthread.so"
#         MORIS_BASE_LIBS 
#         "${MORIS_BASE_LIBS}" )
#     string(REGEX REPLACE 
#         "[^;]*/lib64/liblapack\\.so;[^;]*/lib64/libblas\\.so"
#         "$ENV{MKL_DIR}/lib/intel64_lin/libmkl_intel_lp64.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_core.so;$ENV{MKL_DIR}/lib/intel64_lin/libmkl_sequential.so;/usr/lib64/libpthread.so"
#         MORIS_BASE_LIBS
#         "${MORIS_BASE_LIBS}" )
# endif()
# 
# file(WRITE LDLIBS.txt "${MORIS_BASE_LIBS}")
# -------------------------------------------------------------------------
# MORIS variable handling ------------------------------------------------------
# -------------------------------------------------------------------------

# C++ compiler flags.
set(CMAKE_C_COMPILER ${MORIS_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MORIS_CXX_COMPILER})
set(CMAKE_CXX_FLAGS ${MORIS_CXX_FLAGS})

# Check compiler version.
# Check for GCC version - earlier versions have insufficient C++11
# support, or bugs. Version info is only available in CMake 2.8.0 and
# later.
# if (NOT CMAKE_VERSION VERSION_LESS 2.8.10)
#     if (CMAKE_COMPILER_IS_GNUCXX)
#         if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
#             message(FATAL_ERROR "GCC version must be at least 4.8 (for sufficient C++11  support. You have version ${CMAKE_CXX_COMPILER_VERSION}" )
#         endif()
#     endif()
# endif()

# Linker directories.
link_directories(${MORIS_LDFLAGS})

# Preprocessor definitions.
add_definitions(${MORIS_DEFINITIONS})

# Include directories.
include_directories(${MORIS_INCDIRS})

# Executable.
set(EXECUTABLE "${MORIS}.exe")

message(STATUS "MORIS recognized the C++ flags MORIS_CXX_FLAGS=${MORIS_CXX_FLAGS}.")
message(STATUS "MORIS recognized the C++ linker libraries MORIS_BASE_LIBS=${MORIS_BASE_LIBS}.")
message(STATUS "MORIS recognized the C++ definitions MORIS_DEFINITIONS=${MORIS_DEFINITIONS}.")
message(STATUS "MORIS recognized the C++ linker directories MORIS_LDFLAGS=${MORIS_LDFLAGS}.")
message(STATUS "MORIS recognized the C++ include directories MORIS_INCDIRS=${MORIS_INCDIRS}.")

# Enable verbose output from Makefile builds.
set(CMAKE_VERBOSE_MAKEFILE ON)

# EXECUTABLE_OUTPUT_PATH and LIBRARY_OUTPUT_PATH are deprecated since CMake 2.6.
# http://www.cmake.org/pipermail/cmake-commits/2008-February/003333.html
# set(EXECUTABLE_OUTPUT_PATH "${BIN}" CACHE PATH "Single directory for all executables." )
# set(LIBRARY_OUTPUT_PATH "${LIB}" CACHE PATH "Single directory for all libraries." )
# mark_as_advanced(EXECUTABLE_OUTPUT_PATH)
# mark_as_advanced(LIBRARY_OUTPUT_PATH)

# # Libraries mains must link to; some must be listed more than once
# set(MAIN_LIBS
#     ${ALGORITHMS}
#     ${DISTLINALG}
#     ${GEOMENG}
#     ${OPTIMIZATION}
#     ${MESH}
#     ${MODEL}
#     ${DISTLINALG}
#     ${MESH}
#     ${MORIS_BASE_LIBS}
#     ${COMM}
#     ${TOOLS}
#     ${XTK} )
# 

# -------------------------------------------------------------------------
# Depends files

#include_directories(${MORIS_DEPENDS_DIR})
include(${MORIS_DEPENDS_DIR}/MRS_includes.cmake)

if(MORIS_USE_TESTS)
    include(${MORIS_DEPENDS_DIR}/test_includes.cmake)
endif()

if(BUILD_MAIN)
    include(${MORIS_DEPENDS_DIR}/main_includes.cmake)
endif()

if(BUILD_ALL)
    include(${MORIS_CMAKE_DIR}/utilities/build_all_exe.cmake)
endif()

if(BUILD_ALG)
    include(${MORIS_DEPENDS_DIR}/ALG_Depends.cmake)
endif()

if(BUILD_COM)
    include(${MORIS_DEPENDS_DIR}/COM_Depends.cmake)
endif()

if(BUILD_DLA)
    include(${MORIS_DEPENDS_DIR}/DLA_Depends.cmake)
endif()

if(BUILD_FEM)
    include(${MORIS_DEPENDS_DIR}/FEM_Depends.cmake)
endif()

if(BUILD_GEN)
    include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
endif()

if(BUILD_HMR)
   include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
endif()

if(BUILD_INT)
    include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
endif()

if(BUILD_LNA)
    include(${MORIS_DEPENDS_DIR}/LNA_Depends.cmake)
endif()

if(BUILD_LINALG)
    include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
endif()

if(BUILD_MDL)
    include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
endif()

if(BUILD_MOD)
    include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
endif()

if(BUILD_MSI)
    include(${MORIS_DEPENDS_DIR}/MSI_Depends.cmake)
endif()

if(BUILD_MTK)
    include(${MORIS_DEPENDS_DIR}/MTK_Depends.cmake)
endif()

 if(BUILD_NLA)
    include(${MORIS_DEPENDS_DIR}/NLA_Depends.cmake)
 endif()

if(BUILD_OPT)
    include(${MORIS_DEPENDS_DIR}/OPT_Depends.cmake)
endif()

if(BUILD_STK)
    include(${MORIS_DEPENDS_DIR}/STK_Depends.cmake)
endif()

if(BUILD_SDF)
    include(${MORIS_DEPENDS_DIR}/SDF_Depends.cmake)
endif()

# if(BUILD_TIN)
#    include(${MORIS_DEPENDS_DIR}/TIN_Depends.cmake)
# endif()

if(BUILD_TOL)
    include(${MORIS_DEPENDS_DIR}/TOL_Depends.cmake)
endif()

if(BUILD_XTK)
    include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
endif()


# # -------------------------------------------------------------------------
# # MORIS Config files -----------------------------------------------------------
# # -------------------------------------------------------------------------
# 
# # Needed for the following macros
# include(CMakePackageConfigHelpers)
# 
# # Configure *Config file
# configure_package_config_file(
#     ${MORIS}Config.cmake.in
#     ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS}
#     PATH_VARS MORIS_CXX_FLAGS MORIS_LDFLAGS MORIS_DEFINITIONS MORIS_INCDIRS )
# 
# # Create *ConfigVersion file
# write_basic_package_version_file(
#     ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )
# 
# # Install *Config and *ConfigVersion files
# install(FILES ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake 
#     ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS} )
# 

# -------------------------------------------------------------------------
# MORIS subprojects -------------------------------------------------------
# -------------------------------------------------------------------------
set(EXE_EXT ".exe")
enable_testing()

add_subdirectory(${SRC})


# # -------------------------------------------------------------------------
# # MORIS unit tests -------------------------------------------------------------
# # -------------------------------------------------------------------------
# 
# if(MORIS_USE_TESTS)
#     add_subdirectory(${TEST})
# endif()
# 
# 
# # -------------------------------------------------------------------------
# # MORIS Regression tests -------------------------------------------------------
# # -------------------------------------------------------------------------
# 
# if(MORIS_USE_REGRESSION)
#     add_subdirectory(${REGRESSION})
# endif()

# -------------------------------------------------------------------------
# MORIS Documentation -----------------------------------------------------
# -------------------------------------------------------------------------

add_subdirectory(share/doc)


