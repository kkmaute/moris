# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(EXE_INCLUDES
    ${LINALG}
    ${COM}
    ${HMR}
    ${MTK}
    ${ALG} #MTK
    ${TOL} #MTK
    ${XTK} #MTK
    )

foreach(EXE_INCLUDE ${EXE_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
endforeach()

set(EXE_DEPENDENCIES
    ${HMR}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
    ${ALG}-lib #MTK
    ${TOL}-lib #MTK
#    ${XTK}-lib #MTK
#    ${STK}-lib #MTK
    ${COM}-lib
    ${MORIS_BASE_LIBS}
    ${MORIS_ACML_LAPACK_MKL_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_ARPACK_LIBS}
    ${MORIS_BOOST_LIBS}
    )
    
# Create Tutorial 1
add_executable(HMR_Tutorial_1-exe HMR_Tutorial_1.cpp )
target_link_libraries(HMR_Tutorial_1-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_1-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_1.exe )

# Create Tutorial 2
add_executable(HMR_Tutorial_2-exe HMR_Tutorial_2.cpp )
target_link_libraries(HMR_Tutorial_2-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_2-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_2.exe )

# Create Tutorial 3
add_executable(HMR_Tutorial_3-exe HMR_Tutorial_3.cpp )
target_link_libraries(HMR_Tutorial_3-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_3-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_3.exe )

# Create Tutorial 4
add_executable(HMR_Tutorial_4-exe HMR_Tutorial_4.cpp )
target_link_libraries(HMR_Tutorial_4-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_4-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_4.exe )

# Create Facets 1
add_executable(HMR_Facets_1-exe HMR_Facets_1.cpp )
target_link_libraries(HMR_Facets_1-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Facets_1-exe PROPERTIES OUTPUT_NAME HMR_Facets_1.exe )
