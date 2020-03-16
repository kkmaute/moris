# -------------------------------------------------------------------------
# LBFGSB libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LBFGSB_FOUND_ONCE)
    find_package(LBFGSB)
    
    if(LBFGSB_FOUND)
        set(LBFGSB_FOUND_ONCE TRUE CACHE INTERNAL "LBFGSB was found.")
    endif()
    message(STATUS "LBFGSB_LIBRARIES: ${LBFGSB_LIBRARIES}")
endif()

link_directories(${MORIS_LBFGSB_LIBRARY_DIRS})
set(MORIS_LBFGSB_LIBS ${MORIS_LBFGSB_LIBRARIES})
