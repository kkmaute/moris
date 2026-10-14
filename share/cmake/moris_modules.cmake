# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

include_guard(GLOBAL)

# Registry for BUILD_<module> options, ALL/NONE, and dependency selection.
set(MORIS_BUILD_MODULES
    ALG
    ASR
    CHR
    COM
    CNT
    DLA
    EXC
    FEM
    GEN
    HMR
    INT
    IOS
    LINALG
    MAP
    MOD
    MSI
    MDL
    MTK
    NLA
    OPT
    PRM
    GUI
    SDF
    STK
    TOL
    TSA
    VIS
    WRK
    XTK
    MIG
)

# GUI requires Qt6 and remains opt-in. MOD is disabled in its CMakeLists.txt;
# this checkout has no standalone STK source directory (Trilinos STK is separate).
set(MORIS_ALL_EXCLUSIONS GUI MOD STK)

# These modules are included unconditionally by MRS_includes.cmake.
set(MORIS_CORE_MODULES ASR CHR CNT EXC IOS)
