# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(EXE_INCLUDES
    ${LINALG}
    ${COM}
    ${HMR}
    ${MTK}
    ${MTK}/${MAP}
#    ${ALG} #MTK
#    ${TOL} #MTK
    ${XTK} #MTK
    ${GEN}/${SDF}
    )

foreach(EXE_INCLUDE ${EXE_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
endforeach()

set(EXE_DEPENDENCIES
    ${HMR}-lib
    ${MAP}-lib
    ${MTK}-lib
    ${MDL}-lib # for L2
    ${INT}-lib # for L2
    ${MSI}-lib # for L2
    ${NLA}-lib # for L2
    ${SDF}-lib # for SDF
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

# Create Tutorial 5
add_executable(HMR_Tutorial_5-exe HMR_Tutorial_5.cpp )
target_link_libraries(HMR_Tutorial_5-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_5-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_5.exe )

# For Christian
add_executable(HMR_Christian-exe HMR_Christian.cpp )
target_link_libraries(HMR_Christian-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Christian-exe PROPERTIES OUTPUT_NAME HMR_Christian.exe )

# For Christian
add_executable(HMR_Test-exe HMR_Test.cpp )
target_link_libraries(HMR_Test-exe ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Test-exe PROPERTIES OUTPUT_NAME HMR_Test.exe )

# Create Facets 1 ( uncommented. This needs to go into a test)
#add_executable(HMR_Facets_1-exe HMR_Facets_1.cpp )
#target_link_libraries(HMR_Facets_1-exe ${EXE_DEPENDENCIES} )
#set_target_properties(HMR_Facets_1-exe PROPERTIES OUTPUT_NAME HMR_Facets_1.exe )
