#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# HDF5 libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT HDF5_FOUND_ONCE)
	set(HDF5_ENV_VARS
		$ENV{HDF5DIR}
		$ENV{HDF5_DIR}
		$ENV{hdf5_DIR}
		$ENV{HDF5_ROOT}
		$ENV{hdf5_ROOT}
		$ENV{HDF5_PATH}
		$ENV{hdf5_PATH} )
	
	find_package(HDF5)
	
	if(HDF5_FOUND)
		#set(MORIS_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS}
		#    CACHE PATH "HDF5 include directories." )
		#set(MORIS_HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS}
		#    CACHE PATH "HDF5 library directories." )
		#set(MORIS_HDF5_LIBRARIES ${HDF5_LIBRARIES}
		#    CACHE FILEPATH "HDF5 libraries." )
		set(MORIS_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
		set(MORIS_HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS})
		set(MORIS_HDF5_LIBRARIES ${HDF5_LIBRARIES})
		
	    #set(HDF5_FOUND_ONCE TRUE CACHE INTERNAL "HDF5 was found.")
	    set(HDF5_FOUND_ONCE TRUE)
	endif()
	
	mark_as_advanced(#HDF5_DIR
		MORIS_HDF5_INCLUDE_DIRS
		MORIS_HDF5_LIBRARY_DIRS
		MORIS_HDF5_LIBRARIES )
	
	message(STATUS "HDF5_INCLUDE_DIRS: ${HDF5_INCLUDE_DIRS}")
	message(STATUS "HDF5_LIBRARY_DIRS: ${HDF5_LIBRARY_DIRS}")
	message(STATUS "HDF5_LIBRARIES: ${HDF5_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::hdf5)
	add_library(${MORIS}::hdf5 INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::hdf5 INTERFACE ${MORIS_HDF5_LIBRARIES})
endif()

