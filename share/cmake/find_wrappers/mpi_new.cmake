# -----------------------------------------------------------------------------
# MPI libraries and includes --------------------------------------------------
# -----------------------------------------------------------------------------

if(${MORIS_HAVE_PARALLEL})
    if(${MORIS_USE_MPI} STREQUAL "OPENMPI")
        if(NOT MPI_FOUND_ONCE)
            find_package(MPI)
            if(MPI_FOUND)
                set(MPI_FOUND_ONCE TRUE CACHE INTERNAL "MPI was found.")
                
                set(MORIS_MPI_INCLUDE_DIRS ${MPI_CXX_INCLUDE_PATH}
                	CACHE PATH "OPENMPI include directories." )
                set(MORIS_MPI_LIBRARIES ${MPI_CXX_LIBRARIES}
                	CACHE INTERNAL "OPENMPI libraries." )
                set(MORIS_MPI_DEFINITIONS "-DMORIS_HAVE_PARALLEL"
                	CACHE INTERNAL "Moris preprocessor definitions for OPENMPI." )
                
                mark_as_advanced(MPI_LIBRARY
                	MPI_EXTRA_LIBRARY
                	MORIS_MPI_INCLUDE_DIRS
                	MORIS_MPI_LIBRARIES
                	MORIS_MPI_DEFINITIONS )
                
                set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
        		set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
            endif()
            
            message(STATUS "MPI_C_FOUND: ${MPI_C_FOUND}")
	        message(STATUS "MPI_C_COMPILER: ${MPI_C_COMPILER}")
	        message(STATUS "MPI_C_COMPILE_FLAGS: ${MPI_C_COMPILE_FLAGS}")
	        message(STATUS "MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
	        message(STATUS "MPI_C_LINK_FLAGS: ${MPI_C_LINK_FLAGS}")
	        message(STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
	        
	        message(STATUS "MPI_CXX_FOUND: ${MPI_CXX_FOUND}")
	        message(STATUS "MPI_CXX_COMPILER: ${MPI_CXX_COMPILER}")
	        message(STATUS "MPI_CXX_COMPILE_FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
	        message(STATUS "MPI_CXX_INCLUDE_PATH: ${MPI_CXX_INCLUDE_PATH}")
	        message(STATUS "MPI_CXX_LINK_FLAGS: ${MPI_CXX_LINK_FLAGS}")
	        message(STATUS "MPI_CXX_LIBRARIES: ${MPI_CXX_LIBRARIES}")
        endif()

        if((NOT ${MPI_C_FOUND}) AND (NOT ${MPI_CXX_FOUND}))
            message(FATAL_ERROR "A library with MPI API not found.")
        endif()

        #add_definitions("-DMORIS_HAVE_PARALLEL")
        #include_directories(${MPI_CXX_INCLUDE_PATH})
        #set(MORIS_MPI_LIBS ${MPI_CXX_LIBRARIES})
    else()
        message(FATAL_ERROR "MORIS_USE_MPI supported packages: ${MORIS_MPI_LIBS}")
    endif()
endif()

#mark_as_advanced(MPI_LIBRARY MPI_EXTRA_LIBRARY)
