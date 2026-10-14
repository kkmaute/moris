# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# UMFPACK may be distributed as a static or shared library; do not assume either.
if(NOT TARGET ${MORIS}::suitesparse)
    find_package(SuiteSparse REQUIRED)
    set(MORIS_SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_INCLUDE_DIRS}")
    set(MORIS_SUITESPARSE_LIBRARIES "${SUITESPARSE_LIBRARIES}")
    add_library(${MORIS}::suitesparse INTERFACE IMPORTED GLOBAL)
    target_include_directories(${MORIS}::suitesparse INTERFACE ${MORIS_SUITESPARSE_INCLUDE_DIRS})
    target_link_libraries(${MORIS}::suitesparse INTERFACE ${MORIS_SUITESPARSE_LIBRARIES})
endif()
