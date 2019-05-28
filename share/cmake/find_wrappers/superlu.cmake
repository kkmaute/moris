# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

if(NOT SUPERLU_FOUND_ONCE)
    find_package(SuperLU)
    if(SUPERLU_LIBRARIES) # SuperLU_FOUND should be used here*
        set(SUPERLU_FOUND_ONCE TRUE CACHE INTERNAL "SuperLU was found.")
    endif()
    message(STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
	message(STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")
endif()

include_directories(${SUPERLU_INCLUDES})
set(MORIS_SUPERLU_LIBS ${SUPERLU_LIBRARIES})

# -------------------------------------------------------------------------
# SuperLU_DIST

if(NOT SUPERLU_DIST_FOUND_ONCE)
    find_package(SuperLU_DIST)
    if(SUPERLU_DIST_LIBRARIES) # SuperLU_DIST_FOUND should be used here*
        set(SUPERLU_DIST_FOUND_ONCE TRUE 
            CACHE INTERNAL "SuperLU_DIST was found." )
    endif()
    message(STATUS "SUPERLU_DIST_LIBRARIES: ${SUPERLU_DIST_LIBRARIES}")
endif()

list(APPEND MORIS_SUPERLU_LIBS ${SUPERLU_DIST_LIBRARIES})

# *variable only gets set for one or the other for some reason... 
# Using the librairies variable as a workaround
