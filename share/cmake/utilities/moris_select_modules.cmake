# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

include_guard(GLOBAL)
include("${CMAKE_CURRENT_LIST_DIR}/../moris_modules.cmake")

# Preserve cache descriptions while applying the one-shot ALL/NONE controls.
function(moris_set_build_option name value)
    get_property(description CACHE "${name}" PROPERTY HELPSTRING)
    set(${name} ${value} CACHE BOOL "${description}" FORCE)
endfunction()

macro(moris_select_modules)
    if(BUILD_ALL AND BUILD_NONE)
        message(FATAL_ERROR "BUILD_ALL and BUILD_NONE cannot both be ON.")
    endif()

    # BUILD_EXA is a legacy alias. A changed canonical option takes precedence
    # if both are edited; otherwise allow an alias edit on an existing cache.
    if(DEFINED MORIS_LAST_EXAMPLES_VALUE AND DEFINED BUILD_EXA)
        if("${MORIS_USE_EXAMPLES}" STREQUAL "${MORIS_LAST_EXAMPLES_VALUE}"
           AND NOT "${BUILD_EXA}" STREQUAL "${MORIS_LAST_EXAMPLES_VALUE}"
        )
            moris_set_build_option(MORIS_USE_EXAMPLES "${BUILD_EXA}")
        endif()
    endif()

    if(BUILD_ALL)
        foreach(component IN LISTS MORIS_BUILD_MODULES)
            if(NOT component IN_LIST MORIS_ALL_EXCLUSIONS)
                moris_set_build_option(BUILD_${component} ON)
            elseif(NOT component STREQUAL "GUI")
                moris_set_build_option(BUILD_${component} OFF)
            endif()
        endforeach()
        moris_set_build_option(BUILD_MAIN ON)
        moris_set_build_option(MORIS_USE_EXAMPLES ON)
        message(
            STATUS
                "BUILD_ALL: maintained modules and examples enabled; GUI and tutorials remain opt-in; MOD/STK unavailable."
        )
    elseif(BUILD_NONE)
        foreach(component IN LISTS MORIS_BUILD_MODULES)
            moris_set_build_option(BUILD_${component} OFF)
        endforeach()
        moris_set_build_option(BUILD_MAIN OFF)
        moris_set_build_option(BUILD_TUTORIALS OFF)
        moris_set_build_option(MORIS_USE_EXAMPLES OFF)
    endif()

    # Do not silently accept retired selectors as working build targets.
    foreach(component IN ITEMS MOD STK)
        if(BUILD_${component})
            message(
                WARNING
                    "BUILD_${component} has no active standalone module in this checkout; no ${component} executable will be built."
            )
        endif()
    endforeach()
    moris_set_build_option(BUILD_ALL OFF)
    moris_set_build_option(BUILD_NONE OFF)
    set(BUILD_EXA "${MORIS_USE_EXAMPLES}" CACHE BOOL "Legacy alias for MORIS_USE_EXAMPLES." FORCE)
    mark_as_advanced(BUILD_EXA)
    set(MORIS_LAST_EXAMPLES_VALUE "${MORIS_USE_EXAMPLES}"
        CACHE INTERNAL "Last synchronized examples selection" FORCE
    )
endmacro()
