#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# -------------------------------------------------------------------------
# GCMMA libraries --------------------------------------------------------
# -------------------------------------------------------------------------

if(NOT GCMMA_FOUND_ONCE)
    find_package(GCMMA)

	if(GCMMA_FOUND)
        set(GCMMA_FOUND_ONCE TRUE CACHE INTERNAL "GCMMA was found.")
    endif()

    message(STATUS "GCMMA_LIBRARIES: ${GCMMA_LIBRARIES}")
endif()

include_directories(${GCMMA_INCLUDE_DIRS})
link_directories(${GCMMA_LIBRARY_DIRS})
set(MORIS_GCMMA_LIBS ${GCMMA_LIBRARIES})

