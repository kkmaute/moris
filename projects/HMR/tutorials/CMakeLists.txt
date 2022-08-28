#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# List test includes
set(EXE_INCLUDES
    #${LINALG}
    #${COM}
    #${HMR}
    #${MTK}
    #${MTK}/${MAP}
#    ${ALG} #MTK
#    ${TOL} #MTK
    #${XTK} #MTK
    #${GEN}/${SDF}
    #${FEM}/${INT}
    #${FEM}/${MSI}
    #${FEM}/${MDL}
    )

foreach(EXE_INCLUDE ${EXE_INCLUDES})
    include_directories(${MORIS_PACKAGE_DIR}/${EXE_INCLUDE}/src)
endforeach()

set(EXE_DEPENDENCIES
	${COM}-lib
	${IOS}-lib
    ${HMR}-lib
    ${MTK}-lib
    ${MAP}-lib
    ${SDF}-lib
    #${MDL}-lib # for L2
    #${INT}-lib # for L2
    #${MSI}-lib # for L2
    #${NLA}-lib # for L2
    #${SDF}-lib # for SDF
    #${COM}-lib
    ${MORIS_BASE_LIBS}
    #${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    #${MORIS_ARMADILLO_EIGEN_LIBS}
    #${MORIS_SUPERLU_LIBS}
    #${MORIS_TRILINOS_LIBS}
    #${MORIS_ARPACK_LIBS}
    #${MORIS_BOOST_LIBS}
    ${DLA}-lib #> reeeeeee
    )

# Remove duplicate dependencies
# REVERSE is required to maintain proper linking order
list(REVERSE EXE_DEPENDENCIES)
list(REMOVE_DUPLICATES EXE_DEPENDENCIES)
list(REVERSE EXE_DEPENDENCIES)
    
# Create Tutorial 1
add_executable(HMR_Tutorial_1-exe HMR_Tutorial_1.cpp )
target_link_libraries(HMR_Tutorial_1-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_1-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_1.exe )

# Create Tutorial 2
add_executable(HMR_Tutorial_2-exe HMR_Tutorial_2.cpp )
target_link_libraries(HMR_Tutorial_2-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_2-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_2.exe )

# Create Tutorial 3
add_executable(HMR_Tutorial_3-exe HMR_Tutorial_3.cpp )
target_link_libraries(HMR_Tutorial_3-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_3-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_3.exe )

# Create Tutorial 4
add_executable(HMR_Tutorial_4-exe HMR_Tutorial_4.cpp )
target_link_libraries(HMR_Tutorial_4-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_4-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_4.exe )

# Create Tutorial 5
add_executable(HMR_Tutorial_5-exe HMR_Tutorial_5.cpp )
target_link_libraries(HMR_Tutorial_5-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Tutorial_5-exe PROPERTIES OUTPUT_NAME HMR_Tutorial_5.exe )

# For Christian
add_executable(HMR_Christian-exe HMR_Christian.cpp )
target_link_libraries(HMR_Christian-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Christian-exe PROPERTIES OUTPUT_NAME HMR_Christian.exe )

# For Mathias
add_executable(HMR_Mathias-exe HMR_Mathias.cpp )
target_link_libraries(HMR_Mathias-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Mathias-exe PROPERTIES OUTPUT_NAME HMR_Mathias.exe )

# For Christian
add_executable(HMR_Test-exe HMR_Test.cpp )
target_link_libraries(HMR_Test-exe PRIVATE ${EXE_DEPENDENCIES} )
set_target_properties(HMR_Test-exe PROPERTIES OUTPUT_NAME HMR_Test.exe )

# Create Facets 1 ( uncommented. This needs to go into a test)
#add_executable(HMR_Facets_1-exe HMR_Facets_1.cpp )
#target_link_libraries(HMR_Facets_1-exe PRIVATE ${EXE_DEPENDENCIES} )
#set_target_properties(HMR_Facets_1-exe PROPERTIES OUTPUT_NAME HMR_Facets_1.exe )

set(HMR_TUTORIAL_LIST
	HMR_Tutorial_1-exe
	HMR_Tutorial_2-exe
	HMR_Tutorial_3-exe
	HMR_Tutorial_4-exe
	HMR_Tutorial_5-exe
	HMR_Christian-exe
	HMR_Mathias-exe
	HMR_Test-exe
	)

foreach(HMR_EXE ${HMR_TUTORIAL_LIST})
	foreach(TPL ${HMR_EXE_TPL_DEPENDENCIES})
		string(TOLOWER ${TPL} tpl)
	    include(${MORIS_TPL_DIR}/${tpl}_new.cmake)
	    
	    string(TOUPPER ${TPL} TPL)
	    target_include_directories(${HMR_EXE}
	    	PRIVATE ${MORIS_${TPL}_LIBRARY_DIRS} ${MORIS_${TPL}_INCLUDE_DIRS})
	    target_link_libraries(${HMR_EXE}
	    	PRIVATE ${MORIS_${TPL}_LIBRARIES} ${MORIS_${TPL}_TARGETS})
	    target_compile_definitions(${HMR_EXE}
	    	PRIVATE ${MORIS_${TPL}_DEFINITIONS})
	endforeach()
endforeach()

