# -------------------------------------------------------------------------
# LAPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LAPACK_FOUND_ONCE)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        set(LAPACK_FOUND_ONCE TRUE CACHE INTERNAL "LAPACK was found.")
    endif()
endif()

message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")

# list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_LAPACK")
# list(APPEND MORIS_LDLIBS ${LAPACK_LIBRARIES})

add_definitions("-DMORIS_HAVE_LAPACK")
set(MORIS_LAPACK_LIBS ${LAPACK_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_LIBS ${LAPACK_LIBRARIES})
