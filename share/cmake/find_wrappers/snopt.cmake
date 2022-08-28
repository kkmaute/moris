#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# SNOPT libraries ---------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT SNOPT_FOUND_ONCE)
    find_package(SNOPT)
    
    if(SNOPT_FOUND)
        set(SNOPT_FOUND_ONCE TRUE CACHE INTERNAL "SNOPT was found.")
    endif()
    message(STATUS "SNOPT_LIBRARIES: ${SNOPT_LIBRARIES}")
endif()

link_directories(${MORIS_SNOPT_LIBRARY_DIRS})
set(MORIS_SNOPT_LIBS ${MORIS_SNOPT_LIBRARIES})

