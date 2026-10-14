# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Directory variables shared by module builds and installed package configurations.
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
set(CNT "CNT")
set(COR "COR")
set(DLA "DLA")
set(EXA "EXA")
set(EXC "EXC")
set(FEM "FEM")
set(GEN "GEN")
set(GEN_CORE "GEN_CORE")
set(GEN_MAIN "GEN_MAIN")
set(GUI "GUI")
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
set(ENM "ENM")
set(SDF "SDF")
set(SOL "SOL")
set(SOL_CORE "SOL_CORE")
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

## CMake Paths

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/${SHARE}/${CMAKE}/find_modules"
     "${CMAKE_SOURCE_DIR}/${SHARE}/${CMAKE}/find_wrappers"
)

## MORIS Important Include Directories

list(APPEND MORIS_INCDIRS "${PROJECT_SOURCE_DIR}/${SRC}" "${PROJECT_SOURCE_DIR}/snippets"
     "${PROJECT_SOURCE_DIR}/${INC}"
)
