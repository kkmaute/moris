#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(EXE_INCLUDES
    )

foreach(EXE_INCLUDE ${EXE_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
endforeach()# Include depedency directories   

set(XTK_INCLUDES                                                                                                                                                                                                                                                               
    "assert"                                                                                                                                                                                                                                                                   
    "containers"                                                                                                                                                                                                                                                               
    "core"                                                                                                                                                                                                                                                                     
    "geomeng"                                                                                                                                                                                                                                                                  
    "geometry"                                                                                                                                                                                                                                                                 
    "ios"                                                                                                                                                                                                                                                                      
    "linalg"                                                                                                                                                                                                                                                                   
    "mesh"                                                                                                                                                                                                                                                                     
    "tools"                                                                                                                                                                                                                                                                    
    "topology"                                                                                                                                                                                                                                                                 
    "xtk"                                                                                                                                                                                                                                                                      
    )  
                                                                                                                                                                                                                                             
foreach(XTK_INCLUDE ${XTK_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})                                                                                                                                                                                                           
endforeach()   

set(EXE_DEPENDENCIES
    ${COM}-lib
    ${IOS}-lib
    ${INT}-lib
    ${MDL}-lib
    ${HMR}-lib
    ${MTK}-lib
    ${MAP}-lib
    ${SDF}-lib
    ${XTK}-lib
    ${MORIS_BASE_LIBS}
    ${DLA}-lib #> reeeeeee
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE EXE_DEPENDENCIES)
list(REMOVE_DUPLICATES EXE_DEPENDENCIES)
list(REVERSE EXE_DEPENDENCIES)
    
# Create Tutorial 1
add_executable(tutorial_NASA-exe tutorial_NASA.cpp )
target_link_libraries(tutorial_NASA-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(tutorial_NASA-exe PROPERTIES OUTPUT_NAME tutorial_NASA.exe )

set(MDL_TUTORIAL_LIST
    tutorial_NASA-exe
	)

foreach(MDL_EXE ${MDL_TUTORIAL_LIST})
	foreach(TPL ${MDL_EXE_TPL_DEPENDENCIES})
		string(TOLOWER ${TPL} tpl)
	    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
	    
	    string(TOUPPER ${TPL} TPL)
	    target_include_directories(${MDL_EXE}
	    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
	    target_link_libraries(${MDL_EXE}
	    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
	    target_compile_definitions(${MDL_EXE}
	    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
	endforeach()
endforeach()

