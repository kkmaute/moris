# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

include_guard(GLOBAL)
include(CMakePackageConfigHelpers)

# All modules use the same installed-package contract; target/header installation
# stays beside each module's target definition. Dependencies may contain cycles.
function(moris_install_component_package component)
    cmake_parse_arguments(PACKAGE "" "TARGET" "DEPENDENCIES;TPLS" ${ARGN})
    if(PACKAGE_UNPARSED_ARGUMENTS OR NOT TARGET "${PACKAGE_TARGET}")
        message(FATAL_ERROR "Invalid package arguments for ${component}: ${ARGN}")
    endif()
    set(MORIS_COMPONENT "${component}")
    set(MORIS_HAS_RUNTIME_TARGET FALSE)
    if(TARGET moris_runtime)
        set(MORIS_HAS_RUNTIME_TARGET TRUE)
    endif()
    set(MORIS_COMPONENT_CONFIG_DIR "${${component}_CONFIG_INSTALL_DIR}")
    set(MORIS_COMPONENT_DEPENDENCIES "")
    set(MORIS_COMPONENT_TPLS "")
    foreach(tpl IN LISTS PACKAGE_TPLS)
        string(TOLOWER "${tpl}" tpl)
        string(
            APPEND
            MORIS_COMPONENT_TPLS
            "    include(\"\${_${component}_PREFIX}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n"
        )
    endforeach()

    # A quiet pass creates targets around static dependency cycles before the
    # required pass validates them. Do not interpret raw linker flags as packages.
    foreach(mode IN ITEMS QUIET REQUIRED)
        foreach(dependency IN LISTS PACKAGE_DEPENDENCIES)
            if(dependency MATCHES "^([A-Za-z0-9_]+)-lib$")
                set(dependency "${CMAKE_MATCH_1}")
                if(NOT DEFINED ${dependency}_CONFIG_INSTALL_DIR)
                    message(
                        FATAL_ERROR
                            "Missing package directory for ${dependency} (required by ${component})"
                    )
                endif()
                string(
                    APPEND
                    MORIS_COMPONENT_DEPENDENCIES
                    "    find_package(${dependency} ${mode} CONFIG PATHS \"\${_${component}_PREFIX}/${${dependency}_CONFIG_INSTALL_DIR}\" NO_DEFAULT_PATH)\n"
                )
            endif()
        endforeach()
    endforeach()

    get_target_property(MORIS_COMPONENT_DEFINITIONS "${PACKAGE_TARGET}" COMPILE_DEFINITIONS)
    if(NOT MORIS_COMPONENT_DEFINITIONS)
        set(MORIS_COMPONENT_DEFINITIONS "")
    endif()
    configure_package_config_file(
        "${MORIS_CMAKE_DIR}/morisComponentConfig.cmake.in"
        "${${component}_CONFIG_BUILD_DIR}/${component}Config.cmake"
        INSTALL_DESTINATION "${MORIS_COMPONENT_CONFIG_DIR}"
    )
    write_basic_package_version_file(
        "${${component}_CONFIG_BUILD_DIR}/${component}ConfigVersion.cmake"
        VERSION "${PROJECT_VERSION}" COMPATIBILITY AnyNewerVersion
    )
    install(FILES "${${component}_CONFIG_BUILD_DIR}/${component}Config.cmake"
                  "${${component}_CONFIG_BUILD_DIR}/${component}ConfigVersion.cmake"
            DESTINATION "${MORIS_COMPONENT_CONFIG_DIR}"
    )
    set_property(GLOBAL APPEND PROPERTY MORIS_INSTALLED_COMPONENTS "${component}")
    set_property(
        GLOBAL APPEND PROPERTY MORIS_INSTALLED_COMPONENT_DIRS "${MORIS_COMPONENT_CONFIG_DIR}"
    )
endfunction()
