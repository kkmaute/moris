#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE file in the project root for details.
#
# -------------------------------------------------------------------------
# MORIS CMakeLists.txt ----------------------------------------------------
# -------------------------------------------------------------------------

# NOTES:
# The command 'add_library(${MORIS} STATIC ${SOURCES})' automatically
# adds the 'lib' prefix to the library and the 'a' file extension.
# The extension 'a' is added for a static library.
# Result is 'Linking CXX static library libmoris.a'.

# Require CMake 3.1
cmake_minimum_required(VERSION 3.3.2)

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
set(EXA "EXA")
set(EXC "EXC")
set(FEM "FEM")
#set(GE "GE")
set(GEN "GEN")
set(GEN_CORE "GEN_CORE")
set(GEN_MAIN "GEN_MAIN")
set(HIA "Hierarchical")
set(HMR "HMR")
set(INT "INT")
set(IOS "IOS")
set(MAP "MAP")
set(MDL "MDL")
set(LINALG "LINALG")
set(MOD "MOD")
set(MRS "MRS")
set(MSI "MSI")
set(MTK "MTK")
set(NLA "NLA")
set(OPT "OPT")
set(PRM "PRM")
set(SDF "SDF")
set(SOL "SOL")
set(SOL_CORE "SOL_CORE")
# set(TIN "TIN")
set(TOL "TOL")
set(TSA "TSA")
set(VIS "VIS")
set(WRK "WRK")
set(XTK "XTK")
set(MIG "MIG")

# Common paths
set(MORIS_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MORIS_CMAKE_DIR ${MORIS_DIR}/${SHARE}/${CMAKE})
set(MORIS_DEPENDS_DIR ${MORIS_CMAKE_DIR}/dependencies)
set(MORIS_TPL_DIR ${MORIS_CMAKE_DIR}/find_wrappers)
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

list(APPEND CMAKE_MODULE_PATH 
    "${CMAKE_SOURCE_DIR}/${SHARE}/${CMAKE}/find_modules"
    "${CMAKE_SOURCE_DIR}/${SHARE}/${CMAKE}/find_wrappers" )
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

# MPI execution command
set(MORIS_EXECUTE_COMMAND "mpirun" CACHE STRING "The command used to execute MORIS tests.")

option(MORIS_USE_ARMA
    "Use the Armadillo library as the dense linear algebra package." ON )

option(MORIS_USE_EIGEN
    "Use the Eigen3 library as the dense linear algebra package." OFF )

option(MORIS_USE_32BIT
    "Use 32bit numbers for identifiers and sizes." ON )

option(MORIS_HAVE_DEBUG
    "Use the debug version of the code." OFF )

option(MORIS_USE_VALGRIND
    "Use valgrind to run tests." OFF )
    
option(MORIS_USE_MATRIX_FILL
    "Fill matrices by default with NANs or MaxValues." OFF )

option(MORIS_USE_CHECK_MEMORY
    "Check memory usage for matrices and cells." OFF )
    
option(MORIS_USE_XTK "Have XTK Library" ON)    

option(MORIS_HAVE_SYMBOLIC
    "Use symbolic information in executable." OFF )

option(MORIS_HAVE_SYMBOLIC_STRONG
    "Use symbolic information (strong) in executable." OFF )

option(MORIS_USE_ACML
    "Use the AMD ACML library as the linear algebra package." OFF )

option(MORIS_USE_OPENBLAS
    "Use the open blas library as the linear algebra package." ON )

option(MORIS_USE_MKL
    "Use the Inel Math kernel library as the linear algebra package." OFF )

option(MORIS_USE_LAPACK
    "Use the LAPACK and LAPACK as the linear algebra packages." OFF )

option(MORIS_HAVE_GCMMA
    "Build MORIS with support for GCMMA." ON )

option(MORIS_HAVE_SNOPT
    "Build MORIS with support for SNOPT." ON )
    
option(MORIS_HAVE_LBFGS
    "Build MORIS with support for LBFGS." ON )
    
# option(HAVE_PERF_MAT
#     "Run matrix performance tests." OFF )
#     
# option(HAVE_PERF_SP_MAT
#     "Run sparse matrix performance tests." OFF )
#     
# option(HAVE_PERF_LIN_SOLVE
#     "Run linear solver performance tests." OFF )

option(MORIS_USE_TESTS
    "Compile unit tests." ON )
    
option(MORIS_USE_PERF_TESTS 
    " Build the performance tests" OFF)
    
    
option(MORIS_USE_EXAMPLES
     "Build the examples." ON)
    

option(MORIS_HAVE_PARALLEL_TESTS
    "Run unit tests in parallel. MORIS_USE_TESTS must be ON." ON)

option(MORIS_PERFORM_CHECK
    "Use performance check." ON )
    
option(MORIS_USE_TACC "USE TACC" OFF)

option(MORIS_USE_OPENMP "Use OpenMP" OFF)

# -------------------------------------------------------------------------
# Package Options

option(BUILD_MAIN "Build main executable." OFF)

option(BUILD_ALL "Build all executables." OFF)

option(BUILD_NONE "Build no executables." OFF)

option(BUILD_ALG "Build the algorithms executable." OFF)

option(BUILD_ASR "Build the assert executable." OFF)

option(BUILD_CHR "Build the chronos executable." OFF)

option(BUILD_COM "Build the communication executable." OFF)

option(BUILD_CON "Build the containers executable." OFF)

option(BUILD_DLA "Build the distributed linear algebra executable." OFF)

option(BUILD_EXC "Build the exceptions executable." OFF)

option(BUILD_FEM "Build the FEM executable." OFF)

#option(BUILD_GE "Build the geometry engine executable." OFF)

option(BUILD_GEN "Build the geometry engine executable." OFF)

option(BUILD_HMR "Build the HMR executable." OFF)

option(BUILD_INT "Build the integration and interpolation executable." OFF)

option(BUILD_IOS "Build the IOS executable." OFF)

option(BUILD_LINALG "Build the linear algebra executable." OFF)

option(BUILD_MAP "Build the MTK mapper." OFF)

option(BUILD_MOD "Build the model executable." OFF)

option(BUILD_MSI "Build the model solver interface executable." OFF)

option(BUILD_MDL "Build the model executable." OFF)

option(BUILD_MTK "Build the MTK executable." OFF)

option(BUILD_NLA "Build the non-linear algebra executable." OFF)

option(BUILD_OPT "Build the optimization executable." OFF)

option(BUILD_PRM "Build the PRM executable." OFF)

option(BUILD_SDF "Build the SDF executable." OFF)

option(BUILD_STK "Build the STK executable." OFF)

# option(BUILD_TIN "Build the TIN executable." OFF)

option(BUILD_TOL "Build the tools executable." OFF)

option(BUILD_TSA "Build the time solver algebra executable." OFF)

option(BUILD_TUTORIALS "Build the tutorial executables." OFF)

option(BUILD_VIS "Build the visualization executable." OFF)

option(BUILD_WRK "Build the workflow executable." OFF)

option(BUILD_XTK "Build the XTK executable." OFF)

option(BUILD_MIG "Build the MIG executable." ON)

# -------------------------------------------------------------------------
# Macros
include(${MORIS_CMAKE_DIR}/utilities/moris_macros.cmake)

# ------------------------------------------------------------------------
## Dynamic link test function
include(${MORIS_CMAKE_DIR}/utilities/dynamic_link_input.cmake)

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
    if (MORIS_USE_TACC ) 
        set(MORIS_CXX_COMPILER "mpicxx" )
        if (HAVE_PEDANTIC)
            set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wunused-variable -Wsign-compare" )
        endif()
    else()
        set(MORIS_CXX_COMPILER "mpic++")
    endif()
endif()
set(MORIS_CXX_LINK_EXECUTABLE ${MORIS_CXX_COMPILER})

message(STATUS "MORIS recognized the C compiler type MORIS_C_COMPILER=${MORIS_C_COMPILER}.")
message(STATUS "MORIS recognized the C++ compiler type MORIS_CXX_COMPILER=${MORIS_CXX_COMPILER}.")

# Check for compiler flags.
include(CheckCXXCompilerFlag)

# Allow for mutliple defintions for the same routine
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--allow-multiple-definition")

# set libaries to be linked for all targets
# note: should be done similar to find_modules find_library

# generate search path from LD_LIBRARY_PATH
string(REPLACE ":" ";" LIB_SEARCH_PATH "$ENV{LD_LIBRARY_PATH}")
list(REMOVE_DUPLICATES LIB_SEARCH_PATH)
message(STATUS "Library search path: ${LIB_SEARCH_PATH}")

if (MORIS_USE_TACC)
   find_library(MORIS_LIB_Ifcore NAMES ifcore HINTS ${LIB_SEARCH_PATH})
   message(STATUS "MORIS_LIB_Ifcore = ${MORIS_LIB_Ifcore}")
   link_libraries(${MORIS_LIB_Ifcore})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Ifcore})
   
   find_library(MORIS_LIB_Ifcoremt NAMES ifcoremt HINTS ${LIB_SEARCH_PATH})
   message(STATUS "MORIS_LIB_Ifcoremt = ${MORIS_LIB_Ifcoremt}")
   link_libraries(${MORIS_LIB_Ifcoremt})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Ifcoremt})
   
   find_library(MORIS_LIB_GFortran NAMES libgfortran.so.3 HINTS ${LIB_SEARCH_PATH})
   message(STATUS "MORIS_LIB_GFortran = ${MORIS_LIB_GFortran}")
   link_libraries(${MORIS_LIB_GFortran})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_GFortran})

   find_library(MORIS_LIB_Rt NAMES rt HINTS ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Rt = ${MORIS_LIB_Rt}")
   link_libraries(${MORIS_LIB_Rt})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Rt})
    
   find_library(MORIS_LIB_Dl NAMES dl HINTS ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Dl = ${MORIS_LIB_Dl}")
   link_libraries(${MORIS_LIB_Dl})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Dl})

   find_library(MORIS_LIB_Ssl NAMES ssl HINTS ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Ssl = ${MORIS_LIB_Ssl}")
   link_libraries(${MORIS_LIB_Ssl})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Ssl})

   find_library(MORIS_LIB_Crypto NAMES crypto HINTS ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Crypto = ${MORIS_LIB_Crypto}")
   link_libraries(${MORIS_LIB_Crypto})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Crypto})
   
else(MORIS_USE_TACC)
   
   find_library(MORIS_LIB_Mpi NAMES mpi_mpifh HINTS ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Mpi = ${MORIS_LIB_Mpi}")
   link_libraries(${MORIS_LIB_Mpi})
   #list(APPEND MORIS_BASE_LIBS "${MORIS_LIB_Mpi}")
   
   find_library(MORIS_LIB_GFortran NAMES libgfortran.so.4 HINTS${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_GFortran = ${MORIS_LIB_GFortran}")
   link_libraries(${MORIS_LIB_GFortran})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_GFortran})

   find_library(MORIS_LIB_Dl NAMES dl HINTS  ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Dl = ${MORIS_LIB_Dl}")
   link_libraries(${MORIS_LIB_Dl})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Dl})

   find_library(MORIS_LIB_Ssl NAMES ssl HINTS  ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Ssl = ${MORIS_LIB_Ssl}")
   link_libraries(${MORIS_LIB_Ssl})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Ssl})

   find_library(MORIS_Lib_Crypto NAMES crypto HINTS${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_Lib_Crypto = ${MORIS_Lib_Crypto}")
   link_libraries(${MORIS_Lib_Crypto})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_Lib_Crypto})

   find_library(MORIS_LIB_Z NAMES z HINTS  ${LIB_SEARCH_PATH} )
   message(STATUS "MORIS_LIB_Z = ${MORIS_LIB_Z}")
   link_libraries(${MORIS_LIB_Z})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Z})

   find_library(MORIS_LIB_Sz NAMES sz HINTS  ${LIB_SEARCH_PATH})
   message(STATUS "MORIS_LIB_Sz = ${MORIS_LIB_Sz}")
   link_libraries(${MORIS_LIB_Sz})
   #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_Sz})
endif(MORIS_USE_TACC)

list(APPEND MORIS_DEFINITIONS "-DF77ADD_")

if(MORIS_USE_OPENMP)
    check_cxx_compiler_flag(-fopenmp HAVE_OPENMP)
    message("CHECK FOPENMP")
    if(HAVE_OPENMP)
    message("SETTING FOPENMP")
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -fopenmp")
    endif()
endif()

# -g Request that the compiler and linker generate 
# and retain symbol information in the executable itself,
# for debugging purposes.


# Set optimization type.
# Compile with -O compiler optimizations.
# Debug flags.
if (MORIS_HAVE_DEBUG)
    list(APPEND MORIS_DEFINITIONS "-DDEBUG")
endif()

if (MORIS_HAVE_SYMBOLIC)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type Debug" FORCE)
    
    if (MORIS_HAVE_SYMBOLIC_STRONG)
        # add here option to ease debugging
    else()
        check_cxx_compiler_flag(-Og HAVE_DEBUG)
        if (HAVE_DEBUG)
            set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Og")
         endif()
    endif()

    check_cxx_compiler_flag(-g3 HAVE_DEBUG_3)
    if (HAVE_DEBUG_3)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g3")
    else()
        check_cxx_compiler_flag(-g HAVE_DEBUG)
        if (HAVE_DEBUG)
            set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -g")
        endif()
    endif()
                
    check_cxx_compiler_flag(-ggdb HAVE_DEBUG_DBG)
    if (HAVE_DEBUG_DBG)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -ggdb")
    endif()
else()    
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type Release" FORCE)
    check_cxx_compiler_flag(-O3 HAVE_O3_OPTIMISATION)
    if (HAVE_O3_OPTIMISATION)
        if (MORIS_USE_TACC)
	        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -O3 -xCORE-AVX2 -mkl:sequential -ip")
        else()
            set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -O3 -m64")
        endif()
    endif()
endif()

if (MORIS_USE_MATRIX_FILL)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -DMATRIX_FILL ")
endif()

if (MORIS_USE_CHECK_MEMORY)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -DCHECK_MEMORY ")
endif()

if (MORIS_HAVE_PARALLEL)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -DMORIS_HAVE_PARALLEL ")
endif()

# # Performace test flags
# if (HAVE_PERF_MAT)
#     list(APPEND MORIS_DEFINITIONS "-DPERF_MAT")
# endif()
# 
# if (HAVE_PERF_SP_MAT)
#     list(APPEND MORIS_DEFINITIONS "-DPERF_SP_MAT")
# endif()
# 
# if (HAVE_PERF_LIN_SOLVE)
#     list(APPEND MORIS_DEFINITIONS "-DPERF_LIN_SOLVE")
# endif()

# Performace check flag for logger
if(MORIS_PERFORM_CHECK)
        list(APPEND MORIS_DEFINITIONS "-DMORIS_PERFORM_CHECK")
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
check_cxx_compiler_flag("-Wall -Werror -Wno-unknown-pragmas" HAVE_PEDANTIC )
if (HAVE_PEDANTIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wall -Werror -Wno-unknown-pragmas" )
endif()
check_cxx_compiler_flag("-Wno-long-long" HAVE_PEDANTIC )
if (HAVE_PEDANTIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wno-long-long -Wno-unknown-pragmas" )
endif()
check_cxx_compiler_flag("-pedantic-errors" HAVE_PEDANTIC )
if (HAVE_PEDANTIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -pedantic-errors -Wno-unknown-pragmas" )
endif()
check_cxx_compiler_flag("-Wno-error=maybe-uninitialized" HAVE_PEDANTIC )
if (HAVE_PEDANTIC)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wno-error=maybe-uninitialized -Wno-unknown-pragmas" )
endif()

# Check for C++17 and lower version support.
check_cxx_compiler_flag(-std=c++17 HAVE_STD_CPP17)
if (HAVE_STD_CPP17)
    set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++17")
else()
    check_cxx_compiler_flag(-std=c++14 HAVE_STD_CPP14)
    if (HAVE_STD_CPP14)
        set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++14")
    else()
        check_cxx_compiler_flag(-std=c++11 HAVE_STD_CPP11)
        if (HAVE_STD_CPP11)
            set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++11")
        else()
            check_cxx_compiler_flag(-std=c++0x HAVE_STD_CPP0x)
            if (HAVE_STD_CPP0x)
                set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -std=c++0x")
            endif()
        endif()
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

# Attempt
find_package(Threads)

# -------------------------------------------------------------------------
# GPERFTOOLS : Make sure that paths are set correctly.
# -------------------------------------------------------------------------
option(MORIS_USE_GPERFTOOLS "Use Google Profiling tools (check paths in CMakeLists.txt)." OFF)

if (MORIS_USE_GPERFTOOLS)
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -I/usr/include/gperftools/ -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free" )

    link_libraries("-Wl,--no-as-needed -lprofiler -Wl,--as-needed -ltcmalloc")
endif()

# -------------------------------------------------------------------------
# MEMKIND : for memory profiling.
# -------------------------------------------------------------------------
option(MORIS_USE_MEMORY_DEBUG "Use memkind library for memory debugging)." OFF)

if (MORIS_USE_MEMORY_DEBUG)
    find_library(MORIS_LIB_MemKind NAMES libmemkind.so.0  HINTS /usr/lib64)
    if ( NOT MORIS_LIB_MemKind )
        message( FATAL_ERROR "memkind lib not found - do not use MORIS_USE_MEMORY_DEBUG}")
    endif()

    message(STATUS "MORIS_LIB_MemKind = ${MORIS_LIB_MemKind}")
    link_libraries(${MORIS_LIB_MemKind})
    #list(APPEND MORIS_BASE_LIBS ${MORIS_LIB_MemKind})
endif()

# -------------------------------------------------------------------------
# SACADO : Make sure that paths are set correctly.
# -------------------------------------------------------------------------
option(MORIS_USE_SACADO "UseSacado (check paths in CMakeLists.txt)." OFF)
if (USE_SACADO)
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -I$ENV{HOME}/codes/Sacado/Trilinos-master/packages/sacado/src" )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -I$ENV{HOME}/codes/Sacado/Trilinos-master/packages/sacado/src/mpl" )
    set( MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -I$ENV{HOME}/codes/Sacado/Trilinos-master/packages/sacado/src/parameter" )
endif() #> clean this up eventually

# -------------------------------------------------------------------------
# XTK --------------------------------------------------------------------------
# -------------------------------------------------------------------------   
if (MORIS_USE_XTK)
#    list(APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/include") #> ask if this is separate from the files in MORIS/src/xtk#
#    list(APPEND MORIS_INCDIRS "$ENV{HOME}/codes/xtk/src")
    list(APPEND MORIS_DEFINITIONS "-DUSE_XTK")
    list(APPEND MORIS_DEFINITIONS "-DXTK_USE_MORIS")
endif()


# -------------------------------------------------------------------------
# ACML/MKL/LAPACK libraries -----------------------------------------------
# -------------------------------------------------------------------------

if (MORIS_USE_ACML)
    IF (MORIS_USE_MKL OR MORIS_USE_OPENBLAS OR MORIS_USE_LAPACK)
        message(FATAL_ERROR "Use either MKL or ACML, OPENBLAS, LAPACK")
    endif()
endif()

if (MORIS_USE_OPENBLAS)
    IF (MORIS_USE_MKL OR MORIS_USE_ACML OR MORIS_USE_LAPACK)
        message(FATAL_ERROR "Use either MKL or ACML")
    endif()
endif()

if (MORIS_USE_LAPACK)
    IF (MORIS_USE_MKL OR MORIS_USE_ACML OR MORIS_USE_OPENBLAS)
        message(FATAL_ERROR "Use either LAPACK or ACML")
    endif()
endif()

if (MORIS_USE_ACML)
    #include(${SHARE}/${CMAKE}/acml.cmake)
    set(ACML_LAPACK_MKL_OPENBLAS "acml")
elseif(MORIS_USE_MKL)
    #include(${SHARE}/${CMAKE}/mkl.cmake)
    set(ACML_LAPACK_MKL_OPENBLAS "mkl")
elseif(MORIS_USE_LAPACK)
    #include(${SHARE}/${CMAKE}/lapack.cmake)
    set(ACML_LAPACK_MKL_OPENBLAS "lapack")
elseif(MORIS_USE_OPENBLAS)
    #include(${SHARE}/${CMAKE}/openblas.cmake)
    set(ACML_LAPACK_MKL_OPENBLAS "openblas")
else()
    message(FATAL_ERROR "MORIS only supports the ACML, OPENBLAS, LAPACK, and MKL linear algebra libraries." )
endif()

if (MORIS_USE_32BIT)
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_32BIT")
endif()


# -------------------------------------------------------------------------
# Armadillo ---------------------------------------------------------------
# -------------------------------------------------------------------------

if (MORIS_USE_ARMA AND MORIS_USE_EIGEN)
    message(FATAL_ERROR "Can only set MORIS_USE_EIGEN or MORIS_USE_ARMA; not both" )
endif()
set(ARMADILLO_EIGEN)
if (MORIS_USE_ARMA)
    list(APPEND ARMADILLO_EIGEN "armadillo")
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_ARMA")
    
    if (NOT MORIS_HAVE_DEBUG)
        list(APPEND MORIS_DEFINITIONS "-DARMA_NO_DEBUG")
    endif()
    
#     include(${MORIS_TPL_DIR}/armadillo.cmake)
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
    list(APPEND MORIS_DEFINITIONS "-DMORIS_USE_EIGEN")
endif()

# -------------------------------------------------------------------------
# ARPACK ------------------------------------------------------------------
# -------------------------------------------------------------------------

#include(${MORIS_TPL_DIR}/arpack.cmake)

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
#link_directories(${MORIS_LDFLAGS})

# Preprocessor definitions.
add_definitions(${MORIS_DEFINITIONS})

# Include directories.
include_directories(${MORIS_INCDIRS})

# Executable.
set(EXECUTABLE "${MORIS}.exe")

message(STATUS "MORIS recognized the C++ flags MORIS_CXX_FLAGS=${MORIS_CXX_FLAGS}.")
message(STATUS "MORIS recognized the C++ linker libraries MORIS_BASE_LIBS=${MORIS_BASE_LIBS}.")
message(STATUS "MORIS recognized the C++ definitions MORIS_DEFINITIONS=${MORIS_DEFINITIONS}.")
#message(STATUS "MORIS recognized the C++ linker directories MORIS_LDFLAGS=${MORIS_LDFLAGS}.")
message(STATUS "MORIS recognized the C++ include directories MORIS_INCDIRS=${MORIS_INCDIRS}.")

# Enable verbose output from Makefile builds.
set(CMAKE_VERBOSE_MAKEFILE ON)

# -------------------------------------------------------------------------
# Depends files

#include_directories(${MORIS_DEPENDS_DIR})
include(${MORIS_DEPENDS_DIR}/MRS_includes.cmake)

if(MORIS_USE_TESTS)
    include(${MORIS_DEPENDS_DIR}/test_includes.cmake)
endif()

if(BUILD_MAIN)
    include(${MORIS_DEPENDS_DIR}/main_includes.cmake)
else()
    set(MORIS_USE_EXAMPLES OFF CACHE BOOL "Build the examples." FORCE)
endif()

if(BUILD_ALL)
    include(${MORIS_CMAKE_DIR}/utilities/build_all_exe.cmake)
endif()

if(BUILD_NONE)
    include(${MORIS_CMAKE_DIR}/utilities/build_none_exe.cmake)
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

#if(BUILD_GE)
#    include(${MORIS_DEPENDS_DIR}/GE_Depends.cmake)
#endif()

if(BUILD_GEN)
    include(${MORIS_DEPENDS_DIR}/GEN_Depends.cmake)
endif()

if(BUILD_HMR)
    include(${MORIS_DEPENDS_DIR}/HMR_Depends.cmake)
endif()

if(BUILD_INT)
    include(${MORIS_DEPENDS_DIR}/INT_Depends.cmake)
endif()

if(BUILD_LINALG)
    include(${MORIS_DEPENDS_DIR}/LINALG_Depends.cmake)
endif()

if(BUILD_MAP)
    include(${MORIS_DEPENDS_DIR}/MAP_Depends.cmake)
endif()

if(BUILD_MDL)
    include(${MORIS_DEPENDS_DIR}/MDL_Depends.cmake)
endif()

#if(BUILD_MOD)
#    include(${MORIS_DEPENDS_DIR}/MOD_Depends.cmake)
#endif()

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

if(BUILD_PRM)
    include(${MORIS_DEPENDS_DIR}/PRM_Depends.cmake)
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

if(BUILD_TSA)
    include(${MORIS_DEPENDS_DIR}/TSA_Depends.cmake)
endif()

if(BUILD_VIS)
    include(${MORIS_DEPENDS_DIR}/VIS_Depends.cmake)
endif()

if(BUILD_WRK)
    include(${MORIS_DEPENDS_DIR}/WRK_Depends.cmake)
endif()

if(BUILD_XTK)
    include(${MORIS_DEPENDS_DIR}/XTK_Depends.cmake)
endif()

if(BUILD_MIG)
    include(${MORIS_DEPENDS_DIR}/MIG_Depends.cmake)
endif()

##### Include LINALG Implementation

if(MORIS_USE_ARMA)
    #include(${MORIS_TPL_DIR}/armadillo_new.cmake)
endif()

if(MORIS_USE_EIGEN)
    #include(${MORIS_TPL_DIR}/eigen.cmake)
endif()

# -------------------------------------------------------------------------
# MORIS unit tests --------------------------------------------------------
# -------------------------------------------------------------------------
if(MORIS_HAVE_PARALLEL_TESTS)
    set(MORIS_TEST_PROCS 1 2 4)
    set(MORIS_TEST_PROCS_HMR 1 2 4 12)
endif()

if (MORIS_USE_VALGRIND)
    set(VALGRIND valgrind)
    set(VALGRIND_OPTIONS "")
    set(VALGRIND_OPTIONS_EXA --leak-check=full --track-origins=yes)
else()
    set(VALGRIND "")
    set(VALGRIND_OPTIONS "")
    set(VALGRIND_OPTIONS_EXA "")
endif()

# -------------------------------------------------------------------------
# MORIS examples  --------------------------------------------------------
# -------------------------------------------------------------------------
if(MORIS_USE_EXAMPLES)
    message(STATUS "Build examples")
    include(${MORIS_DEPENDS_DIR}/EXA_Depends.cmake)
endif()

# -------------------------------------------------------------------------
# MORIS subpackages -------------------------------------------------------
# -------------------------------------------------------------------------
set(EXE_EXT ".exe")

enable_testing()

add_subdirectory(${SRC})

# -------------------------------------------------------------------------
# MORIS Find modules and wrappers -----------------------------------------
# -------------------------------------------------------------------------
add_subdirectory(${SHARE}/${CMAKE}/find_modules)
add_subdirectory(${SHARE}/${CMAKE}/find_wrappers)

# -------------------------------------------------------------------------
# MORIS utility scripts ---------------------------------------------------
# -------------------------------------------------------------------------
add_subdirectory(${SHARE}/${CMAKE}/utilities)

# -------------------------------------------------------------------------
# MORIS Documentation -----------------------------------------------------
# -------------------------------------------------------------------------

add_subdirectory(share/doc)

# -------------------------------------------------------------------------
# MORIS Config and Installation -------------------------------------------
# -------------------------------------------------------------------------

# Target to catch general moris properties
add_library(${MORIS}_general INTERFACE)
#target_link_libraries(${MORIS}_general INTERFACE ${MORIS_BASE_LIBS})
#target_include_directories(${MORIS}_general INTERFACE ${MORIS_INCDIRS})
target_compile_definitions(${MORIS}_general INTERFACE ${MORIS_DEFINITIONS})

# Needed for the following macros
include(CMakePackageConfigHelpers)

# Configure *Config file
configure_package_config_file(
    ${MORIS_CMAKE_DIR}/${MORIS}Config.cmake.in
    ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake
    INSTALL_DESTINATION
    ${CONFIG_INSTALL_DIR}/${MORIS}
    PATH_VARS
    CONFIG_INSTALL_DIR )

# Create *ConfigVersion file
write_basic_package_version_file(
    ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# Install *Config and *ConfigVersion files
install(FILES ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}Config.cmake 
    ${CONFIG_BUILD_DIR}/${MORIS}/${MORIS}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS} )

# Install general moris target
install(TARGETS ${MORIS}_general
    DESTINATION ${LIB_INSTALL_DIR}/${MORIS}
    EXPORT ${MORIS}Targets )

install(EXPORT ${MORIS}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MORIS} )

install(FILES "include/catch.hpp"
    DESTINATION ${HEADER_INSTALL_DIR}/${MORIS} )
