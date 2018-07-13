# -----------------------------------------------------------------------------
# PETSc libraries and includes ------------------------------------------------
# -----------------------------------------------------------------------------

#> does not seem to do anything
# if ( USE_TACC )
#     SET(PETSC_PATH "$ENV{PETSC_DIR}")
#     SET(PETSC_DIR  "$ENV{PETSC_DIR}")
#     SET(PETSC_LIB  "$ENV{PETSC_LIB}")
# else()
#     if ( NOT MORIS_HAVE_DEBUG )
#         #SET(PETSC_PATH "$ENV{PETSC_DIR}")
#         #SET(PETSC_DIR  "$ENV{PETSC_DIR}")
#     else()
#         SET(PETSC_PATH "$ENV{HOME}/tpls/petsc-dbg/gcc-openmpi") #> Overwritten by find_package
#         SET(PETSC_DIR  "$ENV{HOME}/tpls/petsc-dbg/gcc-openmpi")
#     endif()
# endif()
# 
# MESSAGE("\nLooking for ${PETSC_PATH}\n\n")

find_package(PETSc)

message(STATUS "PETSC_INCLUDE_DIR: ${PETSC_INCLUDE_DIR}")
message(STATUS "PETSC_LIBRARY_DIR: ${PETSC_LIBRARY_DIR}")
message(STATUS "PETSC_LIBRARY_RELEASE: ${PETSC_LIBRARY_RELEASE}")

list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_PETSC")
list(APPEND MORIS_INCDIRS ${PETSC_INCLUDE_DIR})
list(APPEND MORIS_INCDIRS ${PETSC_ARCH_INCLUDE_DIR})
list(APPEND MORIS_LDLIBS ${PETSC_LIBRARY_RELEASE})
