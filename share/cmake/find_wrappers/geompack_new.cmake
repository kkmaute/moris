#
# Copyright ...
#

# -------------------------------------------------------------------------
# GeomPack libraries ------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GEOPK_FOUND_ONCE)

    find_package(GEOPK)

    if(GEOPK_FOUND)

        set(GEOPK_FOUND_ONCE TRUE)

        set(MORIS_GEOMPACK_INCLUDE_DIRS ${GEOPK_LIBRARY_DIRS})
        set(MORIS_GEOMPACK_LIBRARIES ${GEOPK_LIBRARIES})

        mark_as_advanced(
            MORIS_GEOMPACK_INCLUDE_DIRS
            MORIS_GEOMPACK_LIBRARIES
        )

    endif()

    message(STATUS "GEOPK_LIBRARIES: ${GEOPK_LIBRARIES}")

endif()

if(GEOPK_FOUND AND NOT TARGET ${MORIS}::geompack)

    _import_libraries(
        GEOPK_LIBRARY_TARGETS
        ${GEOPK_LIBRARIES}
    )

    add_library(
        ${MORIS}::geompack
        INTERFACE
        IMPORTED
        GLOBAL
    )

    target_link_libraries(
        ${MORIS}::geompack
        INTERFACE
        ${GEOPK_LIBRARY_TARGETS}
    )

endif()