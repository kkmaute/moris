# -----------------------------------------------------------------------------
# SuperLU libraries and includes ----------------------------------------------
# -----------------------------------------------------------------------------

# if(${MORIS_HAVE_PARALLEL})
#     if(${MORIS_USE_MPI} STREQUAL "OPENMPI")
        find_package(SuperLU)

        message(STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
        message(STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")

        list(APPEND MORIS_DEFINITIONS "-DMORIS_HAVE_SUPERLU")
        list(APPEND MORIS_INCDIRS ${SUPERLU_INCLUDES})
        list(APPEND MORIS_LDLIBS ${SUPERLU_LIBRARIES})
#     else()
#         message(FATAL_ERROR "MORIS_USE_SUPERLU supported packages: ${MORIS_SUPERLU_LIBS}")
#     endif()
# endif()