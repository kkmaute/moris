# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Assert Source

# Public headers.
set(HEADERS assert.hpp fn_assert.hpp fn_log_assert.hpp)

# Assert Library
set(LIB_DEPENDENCIES ${IOS}-lib ${EXC}-lib ${COR}-lib)

add_library(${ASR}-lib INTERFACE)
target_include_directories(
    ${ASR}-lib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                         $<INSTALL_INTERFACE:${${ASR}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${ASR}-lib INTERFACE ${LIB_DEPENDENCIES})

# Assert Executable

# Generate and install the common component package configuration.
moris_install_component_package(
    ASR
    TARGET
    ${ASR}-lib
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${ASR_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${ASR}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${ASR}-lib DESTINATION ${${ASR}_LIB_INSTALL_DIR} EXPORT ${ASR}Targets)

# Install exported target definitions.
install(EXPORT ${ASR}Targets DESTINATION ${${ASR}_CONFIG_INSTALL_DIR})
