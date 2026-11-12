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

# Build a sequence of isolated launches without registering extra CTest entries.
# Bracket quoting preserves spaces, backslashes, and literal CMake variables in
# arguments; generator expressions are resolved by file(GENERATE) below.
function(moris_append_test_group_case script_variable name directory timeout)
    set(invocation "moris_run_test_case(")
    foreach(argument IN ITEMS "${name}" "${directory}" "${timeout}" ${ARGN})
        set(equals "=")
        string(FIND "${argument}" "]${equals}]" closing_bracket)
        while(NOT closing_bracket EQUAL -1)
            string(APPEND equals "=")
            string(FIND "${argument}" "]${equals}]" closing_bracket)
        endwhile()
        string(APPEND invocation " [${equals}[${argument}]${equals}]")
    endforeach()
    string(APPEND invocation ")\n")
    set(${script_variable} "${${script_variable}}${invocation}" PARENT_SCOPE)
endfunction()

function(moris_register_test_group name case_script)
    cmake_parse_arguments(GROUP "" "PROCESSORS;TIMEOUT" "LABELS" ${ARGN})
    set(script "${CMAKE_CURRENT_BINARY_DIR}/${name}-$<CONFIG>.cmake")
    file(
        GENERATE
        OUTPUT "${script}"
        CONTENT
            "include([==[${CMAKE_CURRENT_FUNCTION_LIST_DIR}/moris_run_test_group.cmake]==])\n${case_script}moris_finish_test_group()\n"
    )
    add_test(NAME ${name} COMMAND "${CMAKE_COMMAND}" -P "${script}")
    set_tests_properties(
        ${name}
        PROPERTIES PROCESSORS
                   "${GROUP_PROCESSORS}"
                   TIMEOUT
                   "${GROUP_TIMEOUT}"
                   LABELS
                   "${GROUP_LABELS}"
                   WORKING_DIRECTORY
                   "${CMAKE_CURRENT_BINARY_DIR}"
    )
endfunction()
