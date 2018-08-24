# -------------------------------------------------------------------------
# GCMMA libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GCMMA_FOUND_ONCE)
    find_package(GCMMA)

    message(STATUS "GCMMA_LIBRARIES: ${GCMMA_LIBRARIES}")
    
    if(GCMMA_FOUND)
        set(GCMMA_FOUND_ONCE TRUE CACHE INTERNAL "GCMMA was found.")
    endif()
endif()

# list(APPEND MORIS_INCDIRS "${GCMMA_INCLUDE_DIRS}")
# list(APPEND MORIS_LDFLAGS "${GCMMA_LIBRARY_DIRS}")
# list(APPEND MORIS_LDLIBS "${GCMMA_LIBRARIES}")

include_directories(${GCMMA_INCLUDE_DIRS})
link_directories(${GCMMA_LIBRARY_DIRS})
set(MORIS_GCMMA_LIBS ${GCMMA_LIBRARIES})
