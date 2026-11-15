# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Keep one-process coverage when multi-process testing is disabled.
function(moris_register_unit_tests target)
    set(process_counts 1)
    if(MORIS_HAVE_PARALLEL_TESTS)
        set(process_counts ${ARGN})
    endif()
    foreach(procs IN LISTS process_counts)
        add_test(NAME ${target}-${procs}-procs
                 COMMAND ${MORIS_EXECUTE_COMMAND} -n ${procs} ${VALGRIND} ${VALGRIND_OPTIONS}
                         $<TARGET_FILE:${target}>
        )
        set_tests_properties(${target}-${procs}-procs PROPERTIES PROCESSORS ${procs})
    endforeach()
endfunction()
