#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------

# -------------------------------------------------------------------------
# Kokkos libraries and includes -----------------------------------------
# -------------------------------------------------------------------------

if(NOT KOKKOS_FOUND_ONCE)

     find_package(Kokkos REQUIRED HINTS $ENV{Kokkos_DIR})
 
     if(Kokkos_FOUND)
        set(KOKKOS_FOUND_ONCE TRUE)
        
        set(MORIS_KOKKOS_INCLUDE_DIRS $ENV{Kokkos_DIR}/include CACHE INTERNAL "Kokkos include dir" )
	         
        mark_as_advanced(
            MORIS_KOKKOS_INCLUDE_DIRS
            Kokkos_DIR )
    endif()

endif()

if(NOT TARGET ${MORIS}::kokkos)
	add_library(${MORIS}::kokkos INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::kokkos INTERFACE Kokkos::kokkos)
endif()

# -------------------------------------------------------------------------
# Trilinos libraries and includes -----------------------------------------
# -------------------------------------------------------------------------

if(NOT TRILINOS_FOUND_ONCE) 
    set(TRILINOS_FILE "TrilinosConfig.cmake")

    set(TRILINOS_ENV_VARS
        $ENV{TRILINOSDIR}
        $ENV{TRILINOS_DIR}
        $ENV{Trilinos_DIR}
        $ENV{TRILINOS_ROOT}
        $ENV{Trilinos_ROOT}
        $ENV{TRILINOS_PATH}
        $ENV{Trilinos_PATH} )

    find_path(TRILINOS_DIR 
        NAMES include/${TRILINOS_FILE}
        HINTS
        ${TRILINOS_ENV_VARS}
        PATHS
        /usr/lib/trilinos/gcc-openmpi
        PATH_SUFFIXES
        gcc-openmpi )

    if(NOT TRILINOS_DIR)
        message(FATAL_ERROR 
            "\nPlease set the Trilinos_DIR environment variable. It should be the absolute path to the Trilinos library (ex: /lib/trilinos/gcc-openmpi).\n" )
    endif()

    set(TRILINOS_DEBUG_ENV_VARS
        $ENV{TRILINOSDEBUGDIR}
        $ENV{TRILINOS_DEBUG_DIR}
        $ENV{Trilinos_DEBUG_DIR}
        $ENV{TRILINOS_DEBUG_ROOT}
        $ENV{Trilinos_DEBUG_ROOT}
        $ENV{TRILINOS_DEBUG_PATH}
        $ENV{Trilinos_DEBUG_PATH} )

    find_path(TRILINOS_DEBUG_DIR 
        NAMES include/${TRILINOS_FILE}
        HINTS
        ${TRILINOS_DEBUG_ENV_VARS}
        PATHS
        /usr/lib/trilinos-dbg/gcc-openmpi
        PATH_SUFFIXES
        gcc-openmpi )

    if ( NOT MORIS_HAVE_DEBUG )
        set(TRILINOS_PATH ${TRILINOS_DIR})
    else()
        if(TRILINOS_DEBUG_DIR)
            set(TRILINOS_PATH "${TRILINOS_DEBUG_DIR}")
        else()
            set(TRILINOS_PATH "${TRILINOS_DIR}")
        endif()
    endif()

    MESSAGE(STATUS "\nLooking for ${TRILINOS_PATH}\n\n")

    if(NOT "${TRILINOS_PATH}" MATCHES "${TRILINOS_DIR}")
        set(Trilinos_DIR Trilinos_DIR-NOTFOUND CACHE PATH 
            "Force CMake to find_package(Trilinos ...)"
            FORCE )
    endif()

    FIND_PACKAGE(Trilinos HINTS ${TRILINOS_PATH}/lib/cmake/Trilinos ${TRILINOS_PATH} NO_CMAKE_ENVIRONMENT_PATH NO_SYSTEM_ENVIRONMENT_PATH)
    
    # If FIND_PACKAGE successfully found your Trilinos install, it will
    # set the Boolean flag Trilinos_FOUND.  The following IF statement
    # fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
    # found, it prints out the values of some Trilinos configuration
    # details.  You may find them useful for building your application
    # that uses Trilinos.

    IF(Trilinos_FOUND)

        #  # If Trilinos is included, it may use Pardiso.
        #  if (${MORIS_USE_PARDISO})
        #      if(NOT ${MORIS_USE_MKL})
        # 	       include(${MORIS_CMAKE_DIR}/find_modules/FindMKL.cmake)
        # 	       set(PARDISO_LIBS ${MKL_LIBRARIES})
        # 	       set(PARDISO_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
        #      endif()
        #  endif()
        # 
        # MESSAGE("\nFound Trilinos!  Here are the details: ")
        # MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
        # MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
        # MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
        # MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
        # MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
        # MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
        # MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
        # MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
        # MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
        # MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
        # MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
        # MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
        # MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
        # MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
        # MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
        # MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
        # MESSAGE("   Trilinos_ALL_PACKAGES_TARGETS          = ${Trilinos_ALL_PACKAGES_TARGETS}")
        # MESSAGE("   Trilinos_ALL_SELECTED_PACKAGES_TARGETS = ${Trilinos_ALL_SELECTED_PACKAGES_TARGETS}")
        # MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
        # MESSAGE("End of Trilinos details\n")
        
    ELSE()
        MESSAGE(FATAL_ERROR "Could not find Trilinos!")
    ENDIF()

    MESSAGE("Trilinos was found.")
    set(TRILINOS_FOUND_ONCE TRUE)

    # -------------------------------------------------------------------------

    mark_as_advanced(TRILINOS_DIR Trilinos_DIR TRILINOS_DEBUG_DIR)
    
    # -------------------------------------------------------------------------

    set(MORIS_TRILINOS_LIBS ${MORIS_T_LIBS})

    set(MORIS_TRILINOS_INCLUDE_DIRS
        ${Trilinos_INCLUDE_DIRS}
        ${Trilinos_TPL_INCLUDE_DIRS}
	    ${MORIS_KOKKOS_INCLUDE_DIRS}
        CACHE INTERNAL "Directories included by Trilinos. Very long." )

    mark_as_advanced(
    	MORIS_TRILINOS_LIBS
        MORIS_TRILINOS_INCLUDE_DIRS )
    
    list(REVERSE MORIS_TRILINOS_LIBS)
    list(REMOVE_DUPLICATES MORIS_TRILINOS_LIBS)
    list(REVERSE MORIS_TRILINOS_LIBS)

    set(MORIS_TRILINOS_LIBRARIES ${MORIS_TRILINOS_LIBS}
        CACHE INTERNAL "Trilinos Libraries.")
    
    message(STATUS "TRILINOS_PATH: ${TRILINOS_PATH}")
    
    mark_as_advanced(MORIS_TRILINOS_LIBRARIES)
    
endif()

# -------------------------------------------------------------------------

if(NOT TARGET ${MORIS}::trilinos)
    set(MORIS_TRILINOS_TPLS
        "arpack"
        )
    
    foreach(TPL ${MORIS_TRILINOS_TPLS})
        include(${TPL}_new)
        list(APPEND MORIS_TRILINOS_TPLS_TARGETS ${MORIS}::${TPL})
    endforeach()

    _import_libraries(TRILINOS_LIBRARY_TARGETS "${MORIS_TRILINOS_LIBRARIES}")
    
    _link_each_target("${TRILINOS_LIBRARY_TARGETS}" "${MORIS_TRILINOS_TPLS_TARGETS}")

    add_library(${MORIS}::trilinos INTERFACE IMPORTED GLOBAL)
    target_link_libraries(${MORIS}::trilinos INTERFACE ${TRILINOS_LIBRARY_TARGETS})

    MESSAGE("   Trilinos_ALL_SELECTED_PACKAGES_TARGETS = ${Trilinos_ALL_SELECTED_PACKAGES_TARGETS}")

    target_link_libraries(${MORIS}::trilinos INTERFACE  ${MORIS_TRILINOS_TPLS_TARGETS})
    target_link_libraries(${MORIS}::trilinos INTERFACE  Trilinos::all_libs)
endif()
