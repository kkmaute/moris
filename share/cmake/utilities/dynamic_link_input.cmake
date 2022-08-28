#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

## This CMake function is used to dynamically link the input file
## Name : dynamic_link_input
## Params: path to input file
## Output: .so

function(dynamic_link_input
    	 target_name
		 base_name 
		 cpp_name 
		 so_includes)

set(SO_LIB_REQS
#    ${HMR}-lib
#    ${COM}-lib
#    ${WRK}-lib
#    ${MORIS_PETSC_LIBS}
#    ${MORIS_BOOST_LIBS}
#    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
#    ${MORIS_MPI_LIBS}
#    ${MORIS_ARMADILLO_EIGEN_LIBS}
#    ${MORIS_SUPERLU_LIBS}
#    ${MORIS_LDLIBS}
#    ${MORIS_TRILINOS_LIBS}
#    ${MORIS_BASE_LIBS}
    )


add_library(${target_name} SHARED ${cpp_name})

target_include_directories(${target_name} PRIVATE ${SO_INCLUDES})    

target_link_libraries(${target_name} ${SO_LIB_REQS})

target_compile_definitions(${target_name} INTERFACE ${MORIS_DEFINITIONS})                          

set_target_properties(${target_name} PROPERTIES OUTPUT_NAME ${base_name}
                                                PREFIX      ""          )

endfunction()





