# projects ----------------------------------------------------------------
# -------------------------------------------------------------------------

# Add source directories
foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
    set(${MORIS_SOURCE_DIR}_CONFIG_BUILD_DIR
        "${CONFIG_BUILD_DIR}/${MORIS_SOURCE_DIR}" )
    add_subdirectory(${MORIS_SOURCE_DIR})
endforeach()
<<<<<<< HEAD

# Add header directories
foreach(MORIS_HEADER_DIR ${MORIS_HEADER_DIRS})
    add_subdirectory(${MORIS_HEADER_DIR})
endforeach()

# Add mains
#add_subdirectory(mains)
=======
>>>>>>> c9d0c4bac78ed9c4039e5ee70a9397e4f568c955
