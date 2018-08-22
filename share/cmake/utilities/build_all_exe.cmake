# BUILD_ALL Utility -------------------------------------------------------
# -------------------------------------------------------------------------

set(BUILD_ALG ON CACHE BOOL "Build the algorithms executable." FORCE)

set(BUILD_ASR ON CACHE BOOL "Build the assert executable." FORCE)

set(BUILD_CHR ON CACHE BOOL "Build the chronos executable." FORCE)

set(BUILD_COM ON CACHE BOOL "Build the communication executable." FORCE)

set(BUILD_CON ON CACHE BOOL "Build the containers executable." FORCE)

set(BUILD_DLA ON CACHE BOOL "Build the distributed linear algebra executable." FORCE)

set(BUILD_EXC ON CACHE BOOL "Build the exceptions executable." FORCE)

set(BUILD_FEM ON CACHE BOOL "Build the FEM executable." FORCE)

set(BUILD_GEN ON CACHE BOOL "Build the geometry engine executable." FORCE)

set(BUILD_HMR ON CACHE BOOL "Build the HMR executable." FORCE)

set(BUILD_INT ON CACHE BOOL "Build the integration and interpolation executable." FORCE)

set(BUILD_IOS ON CACHE BOOL "Build the IOS executable." FORCE)

set(BUILD_LNA ON CACHE BOOL "Build the linear algebra executable." FORCE)

set(BUILD_MOD ON CACHE BOOL "Build the model executable." FORCE)

set(BUILD_MSI ON CACHE BOOL "Build the model solver interface executable." FORCE)

set(BUILD_MTK ON CACHE BOOL "Build the MTK executable." FORCE)

set(BUILD_NLA ON CACHE BOOL "Build the non-linear algebra executable." FORCE)

set(BUILD_OPT ON CACHE BOOL "Build the optimization executable." FORCE)

set(BUILD_STK ON CACHE BOOL "Build the STK executable." FORCE)

set(BUILD_TIN ON CACHE BOOL "Build the TIN executable." FORCE)

set(BUILD_TOL ON CACHE BOOL "Build the tools executable." FORCE)

set(BUILD_XTK ON CACHE BOOL "Build the XTK executable." FORCE)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Turn off BUILD_ALL

set(BUILD_ALL OFF CACHE BOOL "Build all executables." FORCE)
