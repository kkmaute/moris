#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# BUILD_NONE Utility ------------------------------------------------------
# -------------------------------------------------------------------------

set(BUILD_ALG OFF CACHE BOOL "Build the algorithms executable." FORCE)

set(BUILD_ASR OFF CACHE BOOL "Build the assert executable." FORCE)

set(BUILD_CHR OFF CACHE BOOL "Build the chronos executable." FORCE)

set(BUILD_COM OFF CACHE BOOL "Build the communication executable." FORCE)

set(BUILD_CON OFF CACHE BOOL "Build the containers executable." FORCE)

set(BUILD_DLA OFF CACHE BOOL "Build the distributed linear algebra executable." FORCE)

set(BUILD_EXC OFF CACHE BOOL "Build the exceptions executable." FORCE)

set(BUILD_FEM OFF CACHE BOOL "Build the FEM executable." FORCE)

set(BUILD_GEN OFF CACHE BOOL "Build the geometry engine executable." FORCE)

set(BUILD_HMR OFF CACHE BOOL "Build the HMR executable." FORCE)

set(BUILD_INT OFF CACHE BOOL "Build the integration and interpolation executable." FORCE)

set(BUILD_IOS OFF CACHE BOOL "Build the IOS executable." FORCE)

set(BUILD_LINALG OFF CACHE BOOL "Build the linear algebra executable." FORCE)

set(BUILD_MDL OFF CACHE BOOL "Build the model executable." FORCE)

#set(BUILD_MOD OFF CACHE BOOL "Build the model executable." FORCE)

set(BUILD_MSI OFF CACHE BOOL "Build the model solver interface executable." FORCE)

set(BUILD_MTK OFF CACHE BOOL "Build the MTK executable." FORCE)

set(BUILD_MAP OFF CACHE BOOL "Build the MTK mapper." FORCE)

set(BUILD_NLA OFF CACHE BOOL "Build the non-linear algebra executable." FORCE)

set(BUILD_OPT OFF CACHE BOOL "Build the optimization executable." FORCE)

set(BUILD_SDF OFF CACHE BOOL "Build the SDF executable." FORCE)

set(BUILD_STK OFF CACHE BOOL "Build the STK executable." FORCE)

# set(BUILD_TIN OFF CACHE BOOL "Build the TIN executable." FORCE)

set(BUILD_TSA OFF CACHE BOOL "Build the TSA executable." FORCE)

set(BUILD_TOL OFF CACHE BOOL "Build the tools executable." FORCE)

set(BUILD_VIS OFF CACHE BOOL "Build the visualization executable." FORCE)

set(BUILD_XTK OFF CACHE BOOL "Build the XTK executable." FORCE)

set(BUILD_WRK OFF CACHE BOOL "Build the workflow executable." FORCE)

set(BUILD_MAIN OFF CACHE BOOL "Build main executable." FORCE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Turn off BUILD_NONE

set(BUILD_NONE OFF CACHE BOOL "Build no executables." FORCE)
