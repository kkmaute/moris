#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# Macros for the moris build system

# creates a list of imported targets for each library in a list
macro(_import_libraries OUTPUT_VAR LIBRARY_LIST)
	set(${OUTPUT_VAR})
	
	foreach(LIB ${LIBRARY_LIST})
        
		string(REGEX REPLACE "^.*/([^/]+).(a|so)" "subtarget_\\1.\\2" LIB_TARGET ${LIB})
		string(REGEX REPLACE "^.*/([^/]+).(a|so)" "\\2" LIB_TYPE ${LIB})

		if(NOT TARGET ${LIB_TARGET})
			if(${LIB_TYPE} STREQUAL "a")
				add_library(${LIB_TARGET} STATIC IMPORTED GLOBAL)
			else()
				add_library(${LIB_TARGET} SHARED IMPORTED GLOBAL)
			endif()
			
			set_target_properties(${LIB_TARGET} PROPERTIES IMPORTED_LOCATION ${LIB})
		endif()
		
		list(APPEND ${OUTPUT_VAR} ${LIB_TARGET})
	endforeach()
endmacro()

# links all targets in A_TARGETS to all targets in B_TARGETS
macro(_link_each_target A_TARGETS B_TARGETS)
	foreach(A_TARGET ${A_TARGETS})
		target_link_libraries(${A_TARGET} INTERFACE ${B_TARGETS})
	endforeach()
endmacro()
