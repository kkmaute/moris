#  Model Source ------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${MDL}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
    cl_MDL_Model.hpp
    cl_MDL_Mesh_Model_Helper.hpp
    )

# List header dependencies
set(MDL_INCLUDES
    #${SOL}/${DLA}
    #${SOL}/${NLA}
    #${SOL}/${TSA}
    #${SOL}/${SOL_CORE}
    #${LINALG}
    #${IOS}
    #${MTK}
    #${COM}
    #${HMR}  # for executable
    #${ALG}
    #${TOL}
    #${FEM}/${INT}
    #${FEM}/${MSI} 
    )

# Include depedency directories
foreach(MDL_INCLUDE ${MDL_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MDL_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Library

# List library source files
set(LIB_SOURCES
    cl_MDL_Model.cpp
    cl_MDL_Mesh_Model_Helper.cpp
    )

# List library dependencies
set(LIB_DEPENDENCIES
	${INT}-lib
	${MTK}-lib
    ${NLA}-lib
    #${SOL_CORE}-lib
    ${TSA}-lib
    #${MSI}-lib
    ${SDF}-lib
    
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_BASE_LIBS})

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${MDL}-lib STATIC ${LIB_SOURCES})
target_include_directories(${MDL}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries(${MDL}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${MDL}-lib PROPERTIES OUTPUT_NAME ${MDL})

foreach(TPL ${MDL_LIB_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${MDL}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${MDL}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${MDL}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Config Files

# # Build *Config and *ConfigVersion files
# include(CMakePackageConfigHelpers)
# 
# configure_package_config_file(
#     ${MDL}Config.cmake.in
#     ${${MDL}_CONFIG_BUILD_DIR}/${MDL}Config.cmake
#     INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${MDL}
#     PATH_VARS LIB_INSTALL_DIR )
# 
# write_basic_package_version_file(
#     ${${MDL}_CONFIG_BUILD_DIR}/${MDL}ConfigVersion.cmake
#     COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# fem/mdl Installation

# # Install *Config and *ConfigVersion files
# install(FILES ${${MDL}_CONFIG_BUILD_DIR}/${MDL}Config.cmake 
#     ${${MDL}_CONFIG_BUILD_DIR}/${MDL}ConfigVersion.cmake
#     DESTINATION ${CONFIG_INSTALL_DIR}/${MDL} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MDL} )

# Install library and executable
install(TARGETS ${MDL}-lib
    EXPORT ${MDL}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MDL}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${MDL}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${MDL} )

# Export library in to *Targets file
install(EXPORT ${MDL}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${MDL} )
