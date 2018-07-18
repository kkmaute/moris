# -------------------------------------------------------------------------
# ViennaCL libraries ------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT VIENNACL_FOUND_ONCE)
    find_package(ViennaCL)
    mark_as_advanced(ViennaCL_DIR)
    if(ViennaCL_FOUND)
        set(VIENNACL_FOUND_ONCE TRUE CACHE INTERNAL 
            "ViennaCL was found." FORCE)
        set(VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIRS} CACHE INTERNAL "ViennaCL include directories" FORCE)
    endif()
endif()

message(STATUS "VIENNACL_LIBRARIES: ${VIENNACL_LIBRARIES}")

# list(APPEND MORIS_VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIRS})

include_directories(${VIENNACL_INCLUDE_DIRS})

