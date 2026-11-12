# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.

set_property(GLOBAL PROPERTY MORIS_FAILED_TEST_CASES "")

function(moris_run_test_case name directory timeout)
    message(STATUS "Running ${name}")
    execute_process(
        COMMAND ${ARGN}
        WORKING_DIRECTORY "${directory}"
        TIMEOUT "${timeout}"
        RESULT_VARIABLE result
    )
    if(NOT "${result}" STREQUAL "0")
        message(STATUS "FAILED ${name}: ${result}")
        set_property(GLOBAL APPEND PROPERTY MORIS_FAILED_TEST_CASES "${name}")
    else()
        message(STATUS "Passed ${name}")
    endif()
endfunction()

function(moris_finish_test_group)
    get_property(failed_cases GLOBAL PROPERTY MORIS_FAILED_TEST_CASES)
    if(failed_cases)
        list(JOIN failed_cases ", " failures)
        message(FATAL_ERROR "The following cases failed: ${failures}")
    endif()
endfunction()
