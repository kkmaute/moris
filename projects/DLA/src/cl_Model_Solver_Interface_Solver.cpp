/*
 * cl_Model_Solver_Interface.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt */

#include "cl_Model_Solver_Interface_Solver.hpp"
#include "cl_Solver_Input.hpp"

using namespace moris;

//Model_Solver_Interface::Model_Solver_Interface( Linear_Solver*     aLin,
//                                                Solver_Input*      aInput )
//{
//    // Get local number of elements
//    moris::uint numLocElements               =    aInput->get_num_my_elements();
//    // Get number of local dofs
//    moris::uint aNumMyDofs                   =    aInput->get_num_my_dofs();
//
//    // Create matrix and RHS temporary in Solver-Model-Interface
//    Sparse_Matrix_Factory<Epetra_FECrsMatrix>      tMatFactory;
//    std::shared_ptr<Sparse_Matrix<Epetra_FECrsMatrix>> tEpetraMat = tMatFactory.create( aNumMyDofs,
//                                                                                        aInput->get_my_local_global_map(),
//                                                                                        aInput->get_constr_dof());
//
//    Dist_Vector_Factory<Epetra_FEVector>      tVecFactory;
//    std::shared_ptr<Dist_Vector<Epetra_FEVector>> tEpetraVectorRHS = tVecFactory.create( aNumMyDofs,
//                                                                                      aInput->get_my_local_global_map(),
//                                                                                      aInput->get_constr_dof());
//
//    std::shared_ptr<Dist_Vector<Epetra_FEVector>> tEpetraVectorLHS = tVecFactory.create( aNumMyDofs,
//                                                                                      aInput->get_my_local_global_map(),
//                                                                                      aInput->get_constr_dof());
//
//    // Loop over all local elements to build matrix graph
//    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
//    {
//        Mat< int > tElementTopology;
//        aInput->get_element_topology(Ii, tElementTopology );
//
//        tEpetraMat->build_graph( tElementTopology.length(), tElementTopology );
//    }
//    // global assembly to switch entries to the right process
//    tEpetraMat->matrix_global_asembly();
//
//    // Loop over all local elements to fill matrix and RHS
//    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
//    {
//        Mat< int > tElementTopology;
//        aInput->get_element_topology(Ii, tElementTopology );
//
//        Mat< real > tElementMatrix;
//        aInput->get_element_matrix(Ii, tElementMatrix );
//
//        Mat< real > tElementRHS;
//        aInput->get_element_rhs(Ii, tElementRHS );
//
//        // Fill element in distributed matrix
//        tEpetraMat->fill_matrix( tElementTopology.length(),
//                                 tElementMatrix,
//                                 tElementTopology);
//
//        // Fill elementRHS in distributed RHS
//        tEpetraVectorRHS->sum_into_global_values( tElementTopology.length(),
//                                               tElementTopology,
//                                               tElementRHS);
//    }
//    // global assembly to switch entries to the right proceccor
//    tEpetraMat->matrix_global_asembly();
//
//    // build linear system on solver class
//    aLin->build_linear_system(tEpetraMat->get_matrix(), tEpetraVectorRHS->get_vector(), tEpetraVectorLHS->get_vector());
//}

//---------------------------------------------------------------------------------------------------------
Model_Solver_Interface::Model_Solver_Interface( moris::Linear_Solver * aLin,
                                                moris::Solver_Input  * aInput,
                                                       Sparse_Matrix * aMat,
                                                moris::Dist_Vector   * aVectorRHS )
{
    // Get local number of elements
    moris::uint numLocElements = aInput->get_num_my_elements();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
    {
        Mat< int > tElementTopology;
        aInput->get_element_topology(Ii, tElementTopology );

        aMat->build_graph( tElementTopology.length(), tElementTopology );
    }
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_asembly();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Mat< int > tElementTopology;
        aInput->get_element_topology( Ii, tElementTopology );

        moris::Mat< moris::real > tElementMatrix;
        aInput->get_element_matrix( Ii, tElementMatrix );

        moris::Mat< moris::real > tElementRHS;
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
    aLin->build_linear_system();
}

//---------------------------------------------------------------------------------------------------------

//Model_Solver_Interface::Model_Solver_Interface( moris::Linear_Solver* aLin,
//                                                Solver_Input*         aInput,
//                                                Matrix_PETSc*         aMarix,
//                                                Dist_Vector *         aVectorRHS)
//{
//    // Get local number of elements
//    moris::uint numLocElements = aInput->get_num_my_elements();
//
//    // Loop over all local elements to fill matrix and RHS
//    for ( moris::uint Ii=0; Ii< numLocElements; Ii++ )
//    {
//        moris::Mat< int > tElementTopology;
//        aInput->get_element_topology( Ii, tElementTopology );
//
//        moris::Mat< moris::real > tElementMatrix;
//        aInput->get_element_matrix( Ii, tElementMatrix );
//
//        moris::Mat< moris::real > tElementRHS;
//        aInput->get_element_rhs( Ii, tElementRHS );
//
//        // Fill element in distributed matrix
//        aMarix->fill_matrix( tElementTopology.length(),
//                             tElementMatrix,
//                             tElementTopology );
//
//        // Fill elementRHS in distributed RHS
//        aVectorRHS->sum_into_global_values( tElementTopology.length(),
//                                            tElementTopology,
//                                            tElementRHS );
//    }
//
//    // global assembly to switch entries to the right proceccor
//    aMarix->matrix_global_asembly();
//    aVectorRHS->vector_global_asembly();
//
//    // build linear system on solver class
//    aLin->build_linear_system();
//}

//---------------------------------------------------------------------------------------------------------

Model_Solver_Interface::~Model_Solver_Interface()
{
}

