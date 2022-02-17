# Parameters Source -------------------------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${WRK}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_WRK_Performer_Manager.hpp
    cl_WRK_GEN_Performer.hpp
    cl_WRK_Performer.hpp
    cl_WRK_Workflow.hpp
    cl_WRK_Workflow_Factory.hpp
    cl_WRK_Workflow_HMR_XTK.hpp
    cl_WRK_Workflow_STK_XTK.hpp
    cl_WRK_Workflow_STK_FEM.hpp
    cl_WRK_perform_refinement.hpp
    cl_WRK_perform_remeshing.hpp
    cl_WRK_DataBase_Performer.hpp
    )
 
set(GEN_INCLUDES
    "field"
    "field/geometry"
    "field/property"
    "pdv"
    )

# Include depedency directories
foreach(GEN_INCLUDES ${GEN_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/GEN/GEN_MAIN/src/${GEN_INCLUDES})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Library

# List library source files
set(LIB_SOURCES
    cl_WRK_Performer_Manager.cpp
    cl_WRK_Performer.cpp
    cl_WRK_GEN_Performer.cpp
    cl_WRK_Workflow_Factory.cpp
    cl_WRK_Workflow_HMR_XTK.cpp
    cl_WRK_Workflow_STK_XTK.cpp
    cl_WRK_Workflow_STK_FEM.cpp
    fn_WRK_Workflow_Main_Interface.cpp
    cl_WRK_perform_refinement.cpp
    cl_WRK_perform_remeshing.cpp
    cl_WRK_DataBase_Performer.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${INT}-lib 
    ${XTK}-lib
    ${HMR}-lib
    ${GEN_MAIN}-lib
    ${MDL}-lib
    ${OPT}-lib
    ${MIG}-lib
     )

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${WRK}-lib STATIC ${LIB_SOURCES})
target_include_directories(${WRK}-lib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${${WRK}_HEADER_INSTALL_DIR}> )
target_link_libraries(${WRK}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${WRK}-lib PROPERTIES OUTPUT_NAME ${WRK})

# Link third party libraries
foreach(TPL ${WRK_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${WRK}-lib
        PUBLIC ${MORIS}::${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${WRK}-lib 
        PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${WRK}-lib
        PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Config Files

# Makes sure all tpl dependencies are searched for in the config file
foreach(TPL ${WRK_TPL_DEPENDENCIES})
    string(TOLOWER ${TPL} tpl)
    string(APPEND WRK_CONFIG_TPL_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${CONFIG_INSTALL_DIR}/find_wrappers/${tpl}_new.cmake\")\n" )
endforeach()

# Makes sure all moris library dependencies are included in the config file
foreach(LIB_TARGET ${LIB_DEPENDENCIES})
    string(REGEX REPLACE "([^-]+).*" "\\1" LIB ${LIB_TARGET})
    string(APPEND WRK_CONFIG_MORIS_DEPENDENCIES
        "include(\"\${PACKAGE_PREFIX_DIR}/${${LIB}_CONFIG_INSTALL_DIR}/${LIB}Config.cmake\")\n" )
endforeach()

# Properties to be copied into config file
get_target_property(WRK_COMPILE_DEFINITIONS ${WRK}-lib COMPILE_DEFINITIONS)

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${WRK}Config.cmake.in
    ${${WRK}_CONFIG_BUILD_DIR}/${WRK}Config.cmake
    INSTALL_DESTINATION
    ${CONFIG_INSTALL_DIR}/${WRK}
    PATH_VARS
    ${WRK}_CONFIG_INSTALL_DIR
    )

write_basic_package_version_file(
    ${${WRK}_CONFIG_BUILD_DIR}/${WRK}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Parameters Installation

# Install *Config and *ConfigVersion files
install(FILES ${${WRK}_CONFIG_BUILD_DIR}/${WRK}Config.cmake 
    ${${WRK}_CONFIG_BUILD_DIR}/${WRK}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${WRK} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${WRK} )

# Install library and executable
install(TARGETS ${WRK}-lib
    EXPORT ${WRK}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${WRK}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${WRK}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${WRK} )

# Export library in to *Targets file
install(EXPORT ${WRK}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${WRK} )
