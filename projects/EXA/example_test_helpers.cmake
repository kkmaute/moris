# Copyright (c) 2026 University of Colorado
# Licensed under the MIT license. See LICENSE.txt in the MORIS root.

# Call after registering this directory's tests and input-library targets.
function(moris_isolate_example_tests driver)
    target_include_directories(${driver} PRIVATE "${MORIS_PACKAGE_DIR}/EXA/src")

    foreach(input IN LISTS ARGN)
        # Drivers and optimization callbacks must find inputs independently of
        # the run directory, including configuration subdirectories.
        target_compile_definitions(
            ${driver} PRIVATE MORIS_EXAMPLE_INPUT_${input}="$<TARGET_FILE:${input}>"
        )
        target_compile_definitions(
            ${input} PRIVATE MORIS_EXAMPLE_INPUT_LIBRARY="$<TARGET_FILE:${input}>"
        )
        add_dependencies(${driver} ${input})
    endforeach()

    get_property(example_tests DIRECTORY PROPERTY TESTS)
    foreach(example_test IN LISTS example_tests)
        set(run_directory "${CMAKE_CURRENT_BINARY_DIR}/runs/$<CONFIG>/${example_test}")
        file(GENERATE OUTPUT "${run_directory}/.moris-test-directory" CONTENT "")
        set_tests_properties(${example_test} PROPERTIES WORKING_DIRECTORY "${run_directory}")

        if(example_test MATCHES "-([0-9]+)-procs$")
            set_tests_properties(${example_test} PROPERTIES PROCESSORS "${CMAKE_MATCH_1}")
        endif()
    endforeach()
endfunction()
