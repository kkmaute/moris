# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

if(NOT SUPERLU_FOUND_ONCE)
    find_package(SuperLU)
    if(SUPERLU_FOUND)
        set(SUPERLU_FOUND_ONCE TRUE CACHE INTERNAL "SuperLU was found.")
    endif()
endif()

message(STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
message(STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")

# # list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_SUPERLU")
# list(APPEND MORIS_INCDIRS ${SUPERLU_INCLUDES})
# list(APPEND MORIS_LDLIBS ${SUPERLU_LIBRARIES})

include_directories(${SUPERLU_INCLUDES})
set(MORIS_SUPERLU_LIBS ${SUPERLU_LIBRARIES})

# -------------------------------------------------------------------------
# SuperLU_DIST

if(NOT SUPERLU_DIST_FOUND_ONCE)
    find_package(SuperLU_DIST)
    if(SUPERLU_DIST_FOUND)
        set(SUPERLU_DIST_FOUND_ONCE TRUE 
            CACHE INTERNAL "SuperLU_DIST was found." )
    endif()
endif()

message(STATUS "SUPERLU_DIST_LIBRARIES: ${SUPERLU_DIST_LIBRARIES}")

# list(APPEND MORIS_LDLIBS ${SUPERLU_DIST_LIBRARIES})

list(APPEND MORIS_SUPERLU_LIBS ${SUPERLU_DIST_LIBRARIES})
