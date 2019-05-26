# -------------------------------------------------------------------------
# ACML libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ACML_FOUND_ONCE)
    find_package(ACML)
    
    if(ACML_FOUND)
        set(ACML_FOUND_ONCE TRUE CACHE INTERNAL "ACML was found." FORCE)
    endif()
    
    message(STATUS "ACML_LIBRARIES: ${ACML_LIBRARIES}")
endif()

add_definitions("-DMORIS_HAVE_ACML")
set(MORIS_ACML_LIBS ${LAPACK_LIBRARIES})
set(MORIS_ACML_LAPACK_MKL_LIBS ${ACML_LIBRARIES})
