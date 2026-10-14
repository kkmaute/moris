# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

# Core Source

# Public headers.
set(HEADERS core.hpp banner.hpp common.hpp moris_typedefs.hpp)

configure_file(paths.hpp.in ${CMAKE_BINARY_DIR}/generated/paths.hpp)
include_directories(${CMAKE_BINARY_DIR}/generated/) # Make sure it can be included...

# Core Library
add_library(${COR}-lib INTERFACE)
target_link_libraries(${COR}-lib INTERFACE moris_runtime)
if(MORIS_USE_GPERFTOOLS)
    target_link_libraries(${COR}-lib INTERFACE ${MORIS}::gperftools)
endif()
target_include_directories(
    ${COR}-lib
    INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
              $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/>
              $<INSTALL_INTERFACE:${${COR}_HEADER_INSTALL_DIR}>
)

# Add memeory profile as a static library to the interface to be linked against
if(MORIS_USE_CHECK_MEMORY)
    add_library(MEM-lib STATIC memoryprofiler.cpp)
    target_link_libraries(${COR}-lib INTERFACE MEM-lib)
endif()
# Core Executable

# Generate and install the common component package configuration.
moris_install_component_package(COR TARGET ${COR}-lib DEPENDENCIES TPLS ${COR_TPL_DEPENDENCIES})

# Install public headers.
install(FILES ${HEADERS} "${CMAKE_BINARY_DIR}/generated/paths.hpp"
        DESTINATION ${${COR}_HEADER_INSTALL_DIR}
)

# Install library
if(MORIS_USE_CHECK_MEMORY)
    install(TARGETS ${COR}-lib MEM-lib DESTINATION ${${COR}_LIB_INSTALL_DIR} EXPORT ${COR}Targets)
else()
    install(TARGETS ${COR}-lib DESTINATION ${${COR}_LIB_INSTALL_DIR} EXPORT ${COR}Targets)
endif()

# Install exported target definitions.
install(EXPORT ${COR}Targets DESTINATION ${${COR}_CONFIG_INSTALL_DIR})
