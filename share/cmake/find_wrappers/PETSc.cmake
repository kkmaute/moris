# -----------------------------------------------------------------------------
# PETSc libraries and includes ------------------------------------------------
# -----------------------------------------------------------------------------

#> does not seem to do anything
# if ( USE_INTEL )
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
        
        set(MORIS_PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} ${PETSC_ARCH_INCLUDE_DIR}
        	CACHE PATH "PETSc include directories." )
        set(MORIS_PETSC_LIBRARY_DIR ${PETSC_LIBRARY_DIR}
        	CACHE PATH "PETSc library directories." )
        set(MORIS_PETSC_LIBRARIES ${PETSC_LIBRARY_RELEASE}
        	CACHE INTERNAL "PETSc libraries.")
        set(MORIS_PETSC_DEFINITIONS "-DMORIS_HAVE_PETSC"
        	CACHE INTERNAL "Moris preprocessor definitions for PETSc." )
        
        mark_as_advanced(MORIS_PETSC_INCLUDE_DIRS
        	MORIS_PETSC_LIBRARY_DIR
        	MORIS_PETSC_LIBRARIES
        	MORIS_PETSC_DEFINITIONS )
    endif()
    message(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
	message(STATUS "PETSC_LIBRARY_DIR: ${PETSC_LIBRARY_DIR}")
	message(STATUS "PETSC_LIBRARY_RELEASE: ${PETSC_LIBRARY_RELEASE}")
endif()

add_definitions("-DMORIS_HAVE_PETSC")
#include_directories(${PETSC_INCLUDE_DIR} ${PETSC_ARCH_INCLUDE_DIR})
#set(MORIS_PETSC_LIBS ${PETSC_LIBRARY_RELEASE})
include_directories(${MORIS_PETSC_INCLUDE_DIRS})
set(MORIS_PETSC_LIBS ${MORIS_PETSC_LIBRARIES})
