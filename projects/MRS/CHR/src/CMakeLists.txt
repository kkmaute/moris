# Chronos Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    chronos.hpp
    cl_Stopwatch.hpp
    cl_Waitbar.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Library

add_library(${CHR}-lib INTERFACE)
target_include_directories(${CHR}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

foreach(TPL ${CHR_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${CHR}-lib
    	INTERFACE ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${CHR}-lib 
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_compile_definitions(${CHR}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Chronos Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${CHR} )

install(TARGETS ${CHR}-lib
	EXPORT ${CHR}Targets
	ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MRS}/${CHR}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${MRS}/${CHR}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${MRS}/${CHR} )

install(EXPORT ${CHR}Targets 
	DESTINATION ${CONFIG_INSTALL_DIR}/${MRS}/${CHR} )
	