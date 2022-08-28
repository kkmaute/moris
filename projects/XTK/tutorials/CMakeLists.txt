#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(XTK_TUTORIAL_INCLUDES
    #${LINALG}
    #${ALG}
    #${MTK}
    #${MTK}/stk_impl
    #${COM}
    )
    
    
foreach(XTK_TUTORIAL_INCLUDES ${XTK_TUTORIAL_INCLUDES})
#	MESSAGE(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
    include_directories(${MORIS_PACKAGE_DIR}/${XTK_TUTORIAL_INCLUDES}/src)
endforeach()

# List additional includes
include_directories(../src)
include_directories(../include)

# List header dependencies                                                                                                                                                                                                                                                     
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

# Include depedency directories                                                                                                                                                                                                                                                
foreach(XTK_INCLUDE ${XTK_INCLUDES})  
    include_directories(${MORIS_PACKAGE_DIR}/XTK/src/${XTK_INCLUDE})                                                                                                                                                                                                           
endforeach()                                                                                                                                                                                                                                                                   

include_directories(${MORIS_PACKAGE_DIR}/XTK/include)    

    
set(TUTORIAL_DEPENDENCIES
    ${COM}-lib
    ${MTK}-lib
    ${ALG}-lib
    ${XTK}-lib
	${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_ARPACK_LIBS}
    )    
    
# Create Tutorial 1
add_executable(XTK_Tutorial_Multi-exe XTK_Tutorial_Multi.cpp )
target_link_libraries(XTK_Tutorial_Multi-exe ${TUTORIAL_DEPENDENCIES} )
set_target_properties(XTK_Tutorial_Multi-exe PROPERTIES OUTPUT_NAME XTK_Tutorial_Multi.exe )


