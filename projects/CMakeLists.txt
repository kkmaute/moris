# projects ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directories
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    set(${MORIS_SOURCE_DIR}_CONFIG_BUILD_DIR
        "${CONFIG_BUILD_DIR}/${MORIS_SOURCE_DIR}" )
    add_subdirectory(${MORIS_SOURCE_DIR})
endforeach()
