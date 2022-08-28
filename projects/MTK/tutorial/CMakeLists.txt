#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# List test includes
set(MTK_TUTORIAL_INCLUDES
    #${LINALG}
    #${COM}
    #${ALG}
    #${SOL}/${DLA}
    #${MTK}
    #${TOL}
    )
    
    
foreach(MTK_TUTORIAL_INCLUDES ${MTK_TUTORIAL_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${MTK_TUTORIAL_INCLUDES}/src)
endforeach()

    
set(TUTORIAL_DEPENDENCIES
    ${COM}-lib
    ${MTK}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    )    
    
# Create Tutorial 1
add_executable(MTK_Tutorial_Data-exe MTK_Tutorial_Data_Mesh.cpp )
target_link_libraries(MTK_Tutorial_Data-exe ${TUTORIAL_DEPENDENCIES} )
set_target_properties(MTK_Tutorial_Data-exe PROPERTIES OUTPUT_NAME MTK_Tutorial_Data.exe )



