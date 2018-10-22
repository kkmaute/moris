# List test includes
set(MTK_TUTORIAL_INCLUDES
    ${LINALG}
    ${COM}
    ${ALG}
    ${DLA}
    ${MTK}
    ${TOL}
    )
    
    
foreach(MTK_TUTORIAL_INCLUDES ${MTK_TUTORIAL_INCLUDES})
	MESSAGE(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
    include_directories(${MORIS_PACKAGE_DIR}/${MTK_TUTORIAL_INCLUDES}/src)
endforeach()

    
set(TUTORIAL_DEPENDENCIES
    ${COM}-lib
    ${MTK}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    )    
    
# Create Tutorial 1
add_executable(MTK_Tutorial_Multi-exe MTK_Tutorial_Data_Mesh.cpp )
target_link_libraries(MTK_Tutorial_Multi-exe ${TUTORIAL_DEPENDENCIES} )
set_target_properties(MTK_Tutorial_Multi-exe PROPERTIES OUTPUT_NAME MTK_Tutorial_Multi.exe )


