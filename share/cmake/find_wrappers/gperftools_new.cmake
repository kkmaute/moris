# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# A single imported target carries profiling requirements into installed consumers.
if(NOT TARGET ${MORIS}::gperftools)
    find_package(Gperftools REQUIRED)
    set(MORIS_GPERFTOOLS_INCLUDE_DIRS "${GPERFTOOLS_INCLUDE_DIRS}")
    set(MORIS_GPERFTOOLS_LIBRARY_DIRS "${GPERFTOOLS_LIBRARY_DIRS}")
    set(MORIS_GPERFTOOLS_LIBRARIES "${GPERFTOOLS_LIBRARIES}")
    add_library(${MORIS}::gperftools INTERFACE IMPORTED GLOBAL)
    target_include_directories(${MORIS}::gperftools INTERFACE ${GPERFTOOLS_INCLUDE_DIRS})
    target_compile_options(
        ${MORIS}::gperftools INTERFACE -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc
                                       -fno-builtin-free
    )
    # Keep the profiler loaded even when only its environment-driven hooks are used.
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_libraries(
            ${MORIS}::gperftools INTERFACE "-Wl,--no-as-needed" ${GPERFTOOLS_LIBRARIES}
                                           "-Wl,--as-needed"
        )
    else()
        target_link_libraries(${MORIS}::gperftools INTERFACE ${GPERFTOOLS_LIBRARIES})
    endif()
endif()
