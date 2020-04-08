## This CMake function is used to dynamically link the input file
## Name : dynamic_link_input
## Params: path to input file
## Output: .so

function(dynamic_link_input
		 base_name 
		 cpp_name 
		 so_name
		 so_includes)

add_library(${so_name} SHARED ${cpp_name})
target_include_directories(${so_name} PRIVATE ${SO_INCLUDES})    
target_link_libraries(${so_name} ${SO_LIB_REQS})
target_compile_definitions(${so_name} INTERFACE ${MORIS_DEFINITIONS})                                                                                                                                                                                                       
set_target_properties(${so_name} PROPERTIES OUTPUT_NAME ${base_name})

endfunction()




