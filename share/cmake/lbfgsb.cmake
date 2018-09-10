# -------------------------------------------------------------------------
# LBFGSB libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LBFGSB_FOUND_ONCE)
    find_package(LBFGSB)
    
    if(LBFGSB_FOUND)
        set(LBFGSB_FOUND_ONCE TRUE CACHE INTERNAL "LBFGSB was found.")
    endif()
endif()

message(STATUS "LBFGSB_LIBRARIES: ${LBFGSB_LIBRARIES}")

# list(APPEND MORIS_LDFLAGS "${LBFGSB_LIBRARY_DIRS}")
# list(APPEND MORIS_LDLIBS "${LBFGSB_LIBRARIES}")

link_directories(${LBFGSB_LIBRARY_DIRS})
set(MORIS_LBFGSB_LIBS ${LBFGSB_LIBRARIES})
