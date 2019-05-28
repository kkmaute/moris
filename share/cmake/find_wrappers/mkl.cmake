# -------------------------------------------------------------------------
# MKL libraries -----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT MKL_FOUND_ONCE)
    find_package(MKL)
    if(MKL_FOUND)
        set(MKL_FOUND_ONCE TRUE CACHE INTERNAL "MKL was found.")
    endif()
    message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
endif()

if(MORIS_USE_MKL)
    add_definitions("-DMORIS_HAVE_MKL")
endif()
include_directories(${MKL_INCLUDE_DIRS})
set(MORIS_MKL_LIBS ${MKL_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_LIBS ${MKL_LIBRARIES})
