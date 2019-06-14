# NonLinear Algebra Source ---------------------------------------
# -------------------------------------------------------------------------

# Set version
set(${NLA}_VERSION ${MORIS_VERSION})

# List header files
set(HEADERS
cl_NLA_Arc_Length.hpp
cl_NLA_Convergence.hpp
cl_NLA_Newton_Solver.hpp
cl_NLA_NLBGS.hpp
cl_NLA_Nonlinear_Algorithm.hpp
cl_NLA_Nonlinear_Problem.hpp
cl_NLA_Nonlinear_Solver_Enums.hpp
cl_NLA_Nonlinear_Solver_Factory.hpp
cl_NLA_Nonlinear_Solver.hpp
)
    
# List header dependencies
set(NLA_INCLUDES
    #${LINALG}
    #${SOL}/${DLA}
    #${SOL}/${SOL_CORE}
    #${NLA}
    #${IOS}
    #${COM}
    #${FEM}/${MSI}
    )

# Include depedency directories
foreach(NLA_INCLUDE ${NLA_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${NLA_INCLUDE}/src)
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Library

# List library source files
set(LIB_SOURCES
    cl_NLA_Arc_Length.cpp   
    cl_NLA_Newton_Solver.cpp
    cl_NLA_Nonlinear_Solver_Factory.cpp
    cl_NLA_Nonlinear_Algorithm.cpp
    cl_NLA_Nonlinear_Solver.cpp
    cl_NLA_Nonlinear_Problem.cpp
    cl_NLA_Convergence.cpp
    cl_NLA_NLBGS.cpp )

# List library dependencies
set(LIB_DEPENDENCIES
    ${LINALG}-lib
    ${DLA}-lib
    ${COM}-lib
    ${MSI}-lib
    ${SOL_CORE}-lib
    #${MORIS_PETSC_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS}
    #${MORIS_ARPACK_LIBS}
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE LIB_DEPENDENCIES)
list(REMOVE_DUPLICATES LIB_DEPENDENCIES)
list(REVERSE LIB_DEPENDENCIES)

# Set the output path for library
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${LIB})

# Create library
add_library(${NLA}-lib STATIC ${LIB_SOURCES})
target_include_directories(${NLA}-lib PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
target_link_libraries(${NLA}-lib PUBLIC ${LIB_DEPENDENCIES})
set_target_properties(${NLA}-lib PROPERTIES OUTPUT_NAME ${NLA})

foreach(TPL ${NLA_LIB_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${NLA}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${NLA}-lib
    	PUBLIC ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${NLA}-lib
    	PUBLIC ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Config Files

# Build *Config and *ConfigVersion files
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${NLA}Config.cmake.in
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake
    INSTALL_DESTINATION ${CONFIG_INSTALL_DIR}/${NLA}
    PATH_VARS LIB_INSTALL_DIR )

write_basic_package_version_file(
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    COMPATIBILITY AnyNewerVersion )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NonLinear Installation

# Install *Config and *ConfigVersion files
install(FILES ${${NLA}_CONFIG_BUILD_DIR}/${NLA}Config.cmake 
    ${${NLA}_CONFIG_BUILD_DIR}/${NLA}ConfigVersion.cmake
    DESTINATION ${CONFIG_INSTALL_DIR}/${NLA} )

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${NLA} )

# Install library and executable
install(TARGETS ${NLA}-lib
    EXPORT ${NLA}Targets
    ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${NLA}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${NLA}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${NLA} )

# Export library in to *Targets file
install(EXPORT ${NLA}Targets
    DESTINATION ${CONFIG_INSTALL_DIR}/${NLA} )
