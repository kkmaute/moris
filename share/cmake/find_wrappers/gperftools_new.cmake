#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# Gperftools find wrappper ------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GPERFTOOLS_FOUND_ONCE)
	find_package(Gperftools)
	
	if(GPERFTOOLS_FOUND)
		set(MORIS_GPERFTOOLS_INCLUDE_DIRS ${GPERFTOOLS_INCLUDE_DIRS})
		set(MORIS_GPERFTOOLS_LIBRARY_DIRS ${GPERFTOOLS_LIBRARY_DIRS})
		set(MORIS_GPERFTOOLS_LIBRARIES ${GPERFTOOLS_LIBRARIES})
		
		set(GPERFTOOLS_FOUND TRUE)
	endif()
	
	message(STATUS "GPERFTOOLS_FOUND: ${GPERFTOOLS_FOUND}")
	message(STATUS "GPERFTOOLS_LIBRARIES: ${MORIS_GPERFTOOLS_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::gperftools)
	_import_libraries(GPERFTOOLS_LIBRARY_TARGETS "${MORIS_GPERFTOOLS_LIBRARIES}")
	
	add_library(${MORIS}::gperftools INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::gperftools INTERFACE ${GPERFTOOLS_LIBRARY_TARGETS})
endif()
