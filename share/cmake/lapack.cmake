# -------------------------------------------------------------------------
# LAPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LAPACK_FOUND_ONCE)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        set(LAPACK_FOUND_ONCE TRUE CACHE INTERNAL "LAPACK was found.")
    endif()
    message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
endif()

add_definitions("-DMORIS_HAVE_LAPACK")
set(MORIS_LAPACK_LIBS ${LAPACK_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_LIBS ${LAPACK_LIBRARIES})
