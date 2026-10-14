# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.

#  SDF Source

# Package version.
set(${SDF}_VERSION ${MORIS_VERSION})

# Public headers.
set(HEADERS
    cl_SDF_Arguments.hpp
    cl_SDF_Cell.hpp
    cl_SDF_Core.hpp
    cl_SDF_Facet_Vertex.hpp
    cl_SDF_Facet.hpp
    cl_SDF_Field.hpp
    cl_SDF_Generator.hpp
    cl_SDF_Line.hpp
    cl_SDF_Mesh.hpp
    cl_SDF_Object.hpp
    cl_SDF_Parameters.hpp
    cl_SDF_State.hpp
    cl_SDF_STK.hpp
    cl_SDF_Triangle.hpp
    cl_SDF_Vertex.hpp
    fn_SDF_Raycast.hpp
    SDF_Tools.hpp
)

# SDF Library

# Library sources.
set(LIB_SOURCES
    cl_SDF_Arguments.cpp
    cl_SDF_Cell.cpp
    cl_SDF_Core.cpp
    cl_SDF_Facet_Vertex.cpp
    cl_SDF_Facet.cpp
    cl_SDF_Field.cpp
    cl_SDF_Generator.cpp
    cl_SDF_Line.cpp
    cl_SDF_Mesh.cpp
    cl_SDF_Object.cpp
    cl_SDF_Parameters.cpp
    cl_SDF_STK.cpp
    cl_SDF_Triangle.cpp
    cl_SDF_Vertex.cpp
    fn_SDF_Raycast.cpp
)

# MORIS library dependencies.
set(LIB_DEPENDENCIES ${LINALG}-lib ${COM}-lib ${CNT}-lib ${MTK}-lib)

# Library output directories.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Define the library and its public build requirements.
moris_add_module(${SDF}-lib SDF_BUILD_TARGET ${LIB_SOURCES})
target_include_directories(
    ${SDF_BUILD_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                               $<INSTALL_INTERFACE:${${SDF}_HEADER_INSTALL_DIR}>
)
target_link_libraries(${SDF_BUILD_TARGET} PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${SDF_BUILD_TARGET} PROPERTIES OUTPUT_NAME ${SDF})

foreach(TPL IN LISTS SDF_TPL_DEPENDENCIES)
    string(TOLOWER ${TPL} tpl)
    include("${MORIS_TPL_DIR}/${tpl}_new.cmake")
    target_link_libraries(${SDF_BUILD_TARGET} PUBLIC ${MORIS}::${tpl})

    string(TOUPPER ${TPL} TPL)
    target_include_directories(
        ${SDF_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS}
    )
    target_compile_definitions(${SDF_BUILD_TARGET} PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# SDF Executable

# List executable source files
set(EXE_SOURCES main.cpp)

# Executable dependencies.
set(EXE_DEPENDENCIES ${CHR}-lib ${SDF}-lib ${MTK}-lib ${COM}-lib ${MORIS_BASE_LIBS})

# Set the output path for executable
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Define the executable and its dependencies.
add_executable(${SDF}-exe ${EXE_SOURCES})
target_link_libraries(${SDF}-exe PRIVATE ${EXE_DEPENDENCIES})
set_target_properties(${SDF}-exe PROPERTIES OUTPUT_NAME sdf)

# Generate and install the common component package configuration.
moris_install_component_package(
    SDF
    TARGET
    ${SDF_BUILD_TARGET}
    DEPENDENCIES
    ${LIB_DEPENDENCIES}
    TPLS
    ${SDF_TPL_DEPENDENCIES}
)

# Install public headers.
install(FILES ${HEADERS} DESTINATION ${${SDF}_HEADER_INSTALL_DIR})

# Install library
install(TARGETS ${SDF}-lib DESTINATION ${${SDF}_LIB_INSTALL_DIR} EXPORT ${SDF}Targets)

# Install executable
install(TARGETS ${SDF}-exe DESTINATION ${${SDF}_EXE_INSTALL_DIR} EXPORT ${SDF}Targets)

# Export library into *Targets file
install(EXPORT ${SDF}Targets DESTINATION ${${SDF}_CONFIG_INSTALL_DIR})
