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

if(NOT PETSC_FOUND_ONCE)
    find_package(PETSc)
    if(PETSC_FOUND)
        set(PETSC_FOUND_ONCE TRUE CACHE INTERNAL "PETSc was found." FORCE)
    endif()
    message(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
	message(STATUS "PETSC_LIBRARY_DIR: ${PETSC_LIBRARY_DIR}")
	message(STATUS "PETSC_LIBRARY_RELEASE: ${PETSC_LIBRARY_RELEASE}")
endif()

add_definitions("-DMORIS_HAVE_PETSC")
include_directories(${PETSC_INCLUDE_DIR} ${PETSC_ARCH_INCLUDE_DIR})
set(MORIS_PETSC_LIBS ${PETSC_LIBRARY_RELEASE})
