# projects ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directories
foreach(MORIS_SRC_DIR ${MORIS_SRC_DIRS})
    set(${MORIS_SRC_DIR}_CONFIG_BUILD_DIR
        "${CONFIG_BUILD_DIR}/${MORIS_SRC_DIR}" )
    add_subdirectory(${MORIS_SRC_DIR})
endforeach()

# Add header directories
foreach(MORIS_HEADER_DIR ${MORIS_HEADER_DIRS})
    add_subdirectory(${MORIS_HEADER_DIR})
endforeach()

# Add mains
# add_subdirectory(mains)
