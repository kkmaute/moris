# BUILD_ALL Utility -------------------------------------------------------
# -------------------------------------------------------------------------

set(BUILD_ALG OFF CACHE BOOL "Build the algorithms executable." FORCE)

set(BUILD_ASR OFF CACHE BOOL "Build the assert executable." FORCE)

set(BUILD_CHR OFF CACHE BOOL "Build the chrOFFos executable." FORCE)

set(BUILD_COM OFF CACHE BOOL "Build the communicatiOFF executable." FORCE)

set(BUILD_COFF OFF CACHE BOOL "Build the cOFFtainers executable." FORCE)

set(BUILD_DLA OFF CACHE BOOL "Build the distributed linear algebra executable." FORCE)

set(BUILD_EXC OFF CACHE BOOL "Build the exceptiOFFs executable." FORCE)

set(BUILD_FEM OFF CACHE BOOL "Build the FEM executable." FORCE)

set(BUILD_GEN OFF CACHE BOOL "Build the geometry engine executable." FORCE)

set(BUILD_HMR OFF CACHE BOOL "Build the HMR executable." FORCE)

set(BUILD_INT OFF CACHE BOOL "Build the integratiOFF and interpolatiOFF executable." FORCE)

set(BUILD_IOS OFF CACHE BOOL "Build the IOS executable." FORCE)

set(BUILD_LINALG OFF CACHE BOOL "Build the linear algebra executable." FORCE)

set(BUILD_MDL OFF CACHE BOOL "Build the model executable." FORCE)

set(BUILD_MOD OFF CACHE BOOL "Build the model executable." FORCE)

set(BUILD_MSI OFF CACHE BOOL "Build the model solver interface executable." FORCE)

set(BUILD_MTK OFF CACHE BOOL "Build the MTK executable." FORCE)

set(BUILD_MAP OFF CACHE BOOL "Build the MTK mapper." FORCE)

set(BUILD_NLA OFF CACHE BOOL "Build the nOFF-linear algebra executable." FORCE)

set(BUILD_OPT OFF CACHE BOOL "Build the optimizatiOFF executable." FORCE)

set(BUILD_SDF OFF CACHE BOOL "Build the SDF executable." FORCE)

set(BUILD_STK OFF CACHE BOOL "Build the STK executable." FORCE)

set(BUILD_TIN OFF CACHE BOOL "Build the TIN executable." FORCE)

set(BUILD_TOL OFF CACHE BOOL "Build the tools executable." FORCE)

set(BUILD_XTK OFF CACHE BOOL "Build the XTK executable." FORCE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Turn off BUILD_ALL

set(BUILD_NONE OFF CACHE BOOL "Build no executables." FORCE)