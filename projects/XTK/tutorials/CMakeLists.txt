# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(XTK_TUTORIAL_INCLUDES
    ${LINALG}
    ${MTK}
    ${MTK}/stk_impl
    ${COM}
    )
    
    
foreach(XTK_TUTORIAL_INCLUDES ${XTK_TUTORIAL_INCLUDES})
	MESSAGE(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
    include_directories(${MORIS_PACKAGE_DIR}/${XTK_TUTORIAL_INCLUDES}/src)
endforeach()

# List additional includes
include_directories(../src)
include_directories(../include)

    
set(TUTORIAL_DEPENDENCIES
    ${MTK}-lib
    ${COM}-lib
	${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    )    
    
# Create Tutorial 1
add_executable(XTK_Tutorial_Multi-exe XTK_Tutorial_Multi.cpp )
target_link_libraries(XTK_Tutorial_Multi-exe ${TUTORIAL_DEPENDENCIES} )
set_target_properties(XTK_Tutorial_Multi-exe PROPERTIES OUTPUT_NAME XTK_Tutorial_Multi.exe )


