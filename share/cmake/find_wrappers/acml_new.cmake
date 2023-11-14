#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# ACML libraries ----------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT ACML_FOUND_ONCE)
    find_package(ACML)
	
    if(ACML_FOUND)
        set(ACML_FOUND_ONCE TRUE)
        
        set(MORIS_ACML_DEFINITIONS "-DMORIS_HAVE_ACML")
        set(MORIS_ACML_LIBRARIES ACML::acml)
        
        mark_as_advanced(MORIS_ACML_LIBRARIES
        	MORIS_ACML_DEFINITIONS )
    endif()
    
    message(STATUS "ACML_LIBRARIES: ${ACML_LIBRARIES}")
endif()

if(NOT TARGET ${MORIS}::acml)
	add_library(${MORIS}::acml INTERFACE IMPORTED GLOBAL)
	target_link_libraries(${MORIS}::acml INTERFACE ${MORIS_ACML_LIBRARIES})
endif()
