# IOS --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    ios.hpp
    cl_Logger.hpp
    fn_Logger.hpp
    fn_to_stdio.hpp
    IO_Tools.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Target
add_library(${IOS}-lib INTERFACE)
target_include_directories(${IOS}-lib INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

foreach(TPL ${IOS_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    target_link_libraries(${IOS}-lib
    	INTERFACE ${tpl})
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${IOS}-lib
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    #target_link_libraries(${IOS}-lib
    #	INTERFACE ${MORIS_${TPL}_LIBRARIES})
    target_compile_definitions(${IOS}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Library

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IOS Installation
# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${IOS} )

