#
# Copyright (c) 2022 University of Colorado Licensed under the MIT license. See
# LICENSE.txt file in the MORIS root for details.
#
# ------------------------------------------------------------------------------------
#

set(EXAMPLE_MAIN Main_Toss_Ruler)
set(EXAMPLE_FILES Toss_Ruler)

# List source files
set(EXAMPLE_SOURCES ../../../src/example_main.cpp example_test_case.cpp)

# List include directories
set(EXAMPLE_INCLUDES ${COM} ${IOS} ${WRK} ${LINALG})

foreach(EXAMPLE_INCLUDE ${EXAMPLE_INCLUDES})
  include_directories(../${EXAMPLE_INCLUDE}/src)
endforeach()

# List executable dependencies
set(EXAMPLE_DEPENDENCIES
    ${HMR}-lib
    ${COM}-lib
    ${WRK}-lib
    ${MORIS_PETSC_LIBS}
    ${MORIS_BOOST_LIBS}
    ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS}
    ${MORIS_SUPERLU_LIBS}
    ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS}
    ${MORIS_BASE_LIBS})

# Set the output path for example
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

# Create executable
add_executable(${EXAMPLE_MAIN} ${EXAMPLE_SOURCES})
target_link_libraries(${EXAMPLE_MAIN} PRIVATE ${EXAMPLE_DEPENDENCIES})
set_target_properties(${EXAMPLE_MAIN} PROPERTIES OUTPUT_NAME
                                                 ${EXAMPLE_MAIN}${EXE_EXT})

# create dynamically linked .so files
set(SO_TPLS "trilinos" ${ARMADILLO_EIGEN})

# get list of INT subfolders
get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)

set(SO_INCLUDES
    ${MORIS_BOOST_INCLUDE_DIRS}
    ${MORIS_TRILINOS_INCLUDE_DIRS}
    ${MORIS_PACKAGE_DIR}/ALG/src
    ${INT_SRC_LIST}
    ${MTK_SRC_LIST}
    ${MORIS_PACKAGE_DIR}/FEM/MSI/src
    ${MORIS_PACKAGE_DIR}/FEM/VIS/src
    ${MORIS_PACKAGE_DIR}/SOL/DLA/src
    ${MORIS_PACKAGE_DIR}/SOL/TSA/src
    ${MORIS_PACKAGE_DIR}/SOL/NLA/src
    ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src
    ${MORIS_PACKAGE_DIR}/LINALG/src
    ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
    ${MORIS_PACKAGE_DIR}/COM/src
    ${MORIS_PACKAGE_DIR}/PRM/ENM/src
    ${MORIS_PACKAGE_DIR}/MTK/src
    ${MORIS_PACKAGE_DIR}/HMR/src
    ${MORIS_PACKAGE_DIR}/XTK/src
    ${MORIS_PACKAGE_DIR}/PRM/src
    ${MORIS_PACKAGE_DIR}/MRS/COR/src)

foreach(TPL ${SO_TPLS})
  string(TOUPPER ${TPL} TPL)
  list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
endforeach()

# create dynamically linked .so files
foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
  dynamic_link_input(${EXAMPLE_FILE} ${EXAMPLE_FILE} ${EXAMPLE_FILE}.cpp
                     ${SO_INCLUDES})
endforeach()

target_compile_definitions(
  ${EXAMPLE_MAIN}
  PRIVATE MORIS_TOSS_LIBRARY="${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Toss_Ruler.so")
target_compile_definitions(Toss_Ruler PRIVATE MORIS_TOSS_DRIVER)
add_dependencies(${EXAMPLE_MAIN} Toss_Ruler)

set(TOSS_PROCESS_COUNTS 1)
if(MORIS_HAVE_PARALLEL_TESTS)
  list(APPEND TOSS_PROCESS_COUNTS 2)
endif()
foreach(PROCS ${TOSS_PROCESS_COUNTS})
  set(TOSS_RUN_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/runs/regression-${PROCS}")
  file(MAKE_DIRECTORY "${TOSS_RUN_DIRECTORY}")
  add_test(
    NAME Toss_Ruler-PlaneStress-${PROCS}-procs
    WORKING_DIRECTORY "${TOSS_RUN_DIRECTORY}"
    COMMAND
      ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS_EXA}
      $<TARGET_FILE:${EXAMPLE_MAIN}> "[TossRuler]")
endforeach()
# The longer, finer 40 ms flight is explicitly selected through Catch, keeping
# normal CTest regression cost separate from the demonstration configuration.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/runs/full")
