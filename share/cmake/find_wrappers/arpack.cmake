# -------------------------------------------------------------------------
# ARPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ARPACK_FOUND_ONCE)
    find_package(ARPACK)
    
    if(ARPACK_FOUND)
        set(ARPACK_FOUND_ONCE TRUE CACHE INTERNAL "ARPACK was found.")
    endif()
    
    message(STATUS "ARPACK_LIBRARIES: ${ARPACK_LIBRARIES}")
endif()

set(MORIS_ARPACK_LIBS ${ARPACK_LIBRARIES})
