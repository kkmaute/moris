# projects ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directories
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    string(REGEX REPLACE "([^\/]*)\/" "" DIR_NAME ${MORIS_SOURCE_DIR})
    set(${DIR_NAME}_CONFIG_BUILD_DIR
        "${CONFIG_BUILD_DIR}/${DIR_NAME}" )
    add_subdirectory(${MORIS_SOURCE_DIR})
endforeach()
