# -------------------------------------------------------------------------
# LAPACK libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT LAPACK_FOUND_ONCE)
    find_package(LAPACK)
    if(LAPACK_FOUND)
        set(LAPACK_FOUND_ONCE TRUE CACHE INTERNAL "LAPACK was found.")
        
        set(MORIS_LAPACK_LIBRARIES ${LAPACK_LIBRARIES}
        	CACHE INTERNAL "LAPACK libraries.")
        set(MORIS_LAPACK_DEFINTIONS "-DMORIS_HAVE_LAPACK"
        	CACHE INTERNAL "Moris preprocessor definitions for LAPACK.")
        
        mark_as_advanced(MORIS_LAPACK_LIBRARIES
        	MORIS_LAPACK_DEFINTIONS )
    endif()
    message(STATUS "LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
endif()

if(NOT TARGET lapack)
	_import_libraries(LAPACK_LIBRARY_TARGETS "${MORIS_LAPACK_LIBRARIES}")

	add_library(lapack INTERFACE IMPORTED GLOBAL)
	target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARY_TARGETS})
endif()

#add_definitions("-DMORIS_HAVE_LAPACK")
#set(MORIS_LAPACK_LIBS ${LAPACK_LIBRARIES})
#set(MORIS_ACML_LAPACK_MKL_LIBS ${LAPACK_LIBRARIES})
