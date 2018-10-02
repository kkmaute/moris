/*
 * cl_Model_Solver_Interface.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt */

#include "cl_Model_Solver_Interface_Solver.hpp"
#include "cl_Solver_Input.hpp"

using namespace moris;

Model_Solver_Interface::Model_Solver_Interface()
{
}

//---------------------------------------------------------------------------------------------------------
Model_Solver_Interface::Model_Solver_Interface( moris::Linear_Solver * aLin,
                                                moris::Solver_Input  * aInput,
                                                moris::Sparse_Matrix * aMat,
                                                moris::Dist_Vector   * aVectorRHS )
{
    // Get local number of elements
    moris::uint numLocElements = aInput->get_num_my_elements();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology(Ii, tElementTopology );

        aMat->build_graph( tElementTopology.length(), tElementTopology );
    }
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_asembly();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        aInput->get_element_matrix( Ii, tElementMatrix );

        Matrix< DDRMat > tElementRHS;
        aInput->get_element_rhs( Ii, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix( tElementTopology.length(),
                                 tElementMatrix,
                                 tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values( tElementTopology.length(),
                                            tElementTopology,
                                            tElementRHS );
    }
    // global assembly to switch entries to the right proceccor
    aVectorRHS->vector_global_asembly();
    aMat->matrix_global_asembly();

    // build linear system on solver class
    //aLin->build_linear_system();
}

//---------------------------------------------------------------------------------------------------------
void Model_Solver_Interface::build_graph( moris::Solver_Input  * aInput,
                                          moris::Sparse_Matrix * aMat )
{
    // Get local number of elements
    moris::uint numLocElements = aInput->get_num_my_elements();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology(Ii, tElementTopology );

        aMat->build_graph( tElementTopology.length(), tElementTopology );
    }
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_asembly();
}

//---------------------------------------------------------------------------------------------------------
void Model_Solver_Interface::fill_matrix_and_RHS( moris::Linear_Solver * aLin,
                                                  moris::Solver_Input  * aInput,
                                                  moris::Sparse_Matrix * aMat,
                                                  moris::Dist_Vector   * aVectorRHS,
                                                  moris::Dist_Vector   * aFullSolutionVector )
{
    aInput->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numLocElements = aInput->get_num_my_elements();

//    // Loop over all local elements to build matrix graph
//    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
//    {
//        Matrix< DDSMat > tElementTopology;
//        aInput->get_element_topology(Ii, tElementTopology );
//
//        aMat->build_graph( tElementTopology.length(), tElementTopology );
//    }
//    // global assembly to switch entries to the right proceccor
//    aMat->matrix_global_asembly();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        aInput->get_element_matrix( Ii, tElementMatrix );

        Matrix< DDRMat > tElementRHS;
        aInput->get_element_rhs( Ii, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix( tElementTopology.length(),
                                 tElementMatrix,
                                 tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values( tElementTopology.length(),
                                            tElementTopology,
                                            tElementRHS );
    }
    // global assembly to switch entries to the right proceccor
    aVectorRHS->vector_global_asembly();
    aMat->matrix_global_asembly();

    // build linear system on solver class
    //aLin->build_linear_system();
}

//---------------------------------------------------------------------------------------------------------
void Model_Solver_Interface::fill_matrix_and_RHS( moris::Linear_Solver * aLin,
                                                  moris::Solver_Input  * aInput,
                                                  moris::Sparse_Matrix * aMat,
                                                  moris::Dist_Vector   * aVectorRHS )
{
    // Get local number of elements
    moris::uint numLocElements = aInput->get_num_my_elements();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
    {
        Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology(Ii, tElementTopology );

        aMat->build_graph( tElementTopology.length(), tElementTopology );
    }
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_asembly();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        aInput->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        aInput->get_element_matrix( Ii, tElementMatrix );

        Matrix< DDRMat > tElementRHS;
        aInput->get_element_rhs( Ii, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix( tElementTopology.length(),
                                 tElementMatrix,
                                 tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values( tElementTopology.length(),
                                            tElementTopology,
                                            tElementRHS );
    }
    // global assembly to switch entries to the right proceccor
    aVectorRHS->vector_global_asembly();
    aMat->matrix_global_asembly();

    // build linear system on solver class
    //aLin->build_linear_system();
}
//---------------------------------------------------------------------------------------------------------
Model_Solver_Interface::~Model_Solver_Interface()
{
}

