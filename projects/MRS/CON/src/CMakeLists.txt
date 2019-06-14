# Containers Source --------------------------------------------------------------
# -------------------------------------------------------------------------

# List header files
set(HEADERS
    cl_Array.hpp
    cl_Bi_Map.hpp
    cl_Bitbool.hpp
    cl_Bitset.hpp
    cl_BoostBitset.hpp
    cl_Cell.hpp
    cl_Dist_Map.hpp
    cl_Map.hpp
    cl_Param_List.hpp
    cl_Tuple.hpp
    containers.hpp
    fn_zip.hpp )

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Library

add_library(${CON}-lib INTERFACE)
target_include_directories(${CON}-lib INTERFACE
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

foreach(TPL ${CON_TPL_DEPENDENCIES})
	string(TOLOWER ${TPL} tpl)
    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
    
    string(TOUPPER ${TPL} TPL)
    target_include_directories(${CON}-lib 
    	INTERFACE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
    target_link_libraries(${CON}-lib
    	INTERFACE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
    target_compile_definitions(${CON}-lib
    	INTERFACE ${MORIS_${TPL}_DEFINITIONS})
endforeach()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Executable

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Config Files

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Containers Installation

# Install header files
install(FILES ${HEADERS}
    DESTINATION ${HEADER_INSTALL_DIR}/${MRS}/${CON} )

install(TARGETS ${CON}-lib
	EXPORT ${CON}Targets
	ARCHIVE DESTINATION ${LIB_INSTALL_DIR}/${MRS}/${CON}
    LIBRARY DESTINATION ${LIB_INSTALL_DIR}/${MRS}/${CON}
    RUNTIME DESTINATION ${EXE_INSTALL_DIR}/${MRS}/${CON} )

install(EXPORT ${CON}Targets 
	DESTINATION ${CONFIG_INSTALL_DIR}/${MRS}/${CON} )
