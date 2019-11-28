/*
 * cl_DLA_Solver_Interface.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt */

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_Sparse_Matrix.hpp"
#include "cl_Vector.hpp"

using namespace moris;

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::build_graph( moris::Sparse_Matrix * aMat )
{
    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationOnjOnBlock = this->get_num_my_elements_on_block( Ii );

        for ( moris::uint Ik=0; Ik < tNumEquationOnjOnBlock; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

            aMat->build_graph( tElementTopology.length(), tElementTopology );
        }
    }
//    aMat->print();
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::fill_matrix_and_RHS( moris::Sparse_Matrix * aMat,
                                            moris::Dist_Vector   * aVectorRHS,
                                            moris::Dist_Vector   * aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numLocElements = this->get_num_my_elements();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        this->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        this->get_element_matrix( Ii, tElementMatrix );

        Matrix< DDRMat > tElementRHS;
        this->get_element_rhs( Ii, tElementRHS );

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
    aMat->matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::assemble_RHS( moris::Dist_Vector * aVectorRHS,
                                     moris::Dist_Vector * aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

//    aFullSolutionVector->print();

    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationOnjOnBlock = this->get_num_my_elements_on_block( Ii );

        this->initialize_block( Ii );

        for ( moris::uint Ik=0; Ik < tNumEquationOnjOnBlock; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

//            print(tElementTopology,"tElementTopology");

            Matrix< DDRMat > tElementRHS;
            this->get_element_rhs( Ii, Ik, tElementRHS );

//            print(tElementRHS,"tElementRHS");

            // Fill elementRHS in distributed RHS
            aVectorRHS->sum_into_global_values( tElementTopology.length(),
                                                tElementTopology,
                                                tElementRHS );
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right proceccor
    aVectorRHS->vector_global_asembly();

//    std::cout<<"Assembled Residual Vector"<<std::endl;
//    aVectorRHS->print();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::assemble_jacobian( moris::Sparse_Matrix * aMat,
                                          moris::Dist_Vector   * aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

//#ifdef WITHGPERFTOOLS
//     ProfilerStart("./main.prof");
//#endif

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        std::cout<<"Block "<<Ii<<std::endl;
        moris::uint tNumEquationOnjOnBlock = this->get_num_my_elements_on_block( Ii );

        this->initialize_block( Ii );

        for ( moris::uint Ik=0; Ik < tNumEquationOnjOnBlock; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

//            print(tElementTopology,"tElementTopology");

            Matrix< DDRMat > tElementMatrix;
            this->get_element_matrix( Ii, Ik, tElementMatrix );

//            print(tElementMatrix,"tElementMatrix");

            // Fill element in distributed matrix
            aMat->fill_matrix( tElementTopology.length(),
                               tElementMatrix,
                               tElementTopology );
        }
        aMat->matrix_global_assembly();
        this->free_block_memory( Ii );
    }
    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_assembly();

    aMat->save_matrix_to_matlab_file( "Matrix.dat");

//#ifdef WITHGPERFTOOLS
//    ProfilerStop();
//#endif

//    aMat->print();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::fill_matrix_and_RHS( moris::Sparse_Matrix * aMat,
                                            moris::Dist_Vector   * aVectorRHS )
{
    // Get local number of elements
    moris::uint numLocElements = this->get_num_my_elements();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        this->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        this->get_element_matrix( Ii, tElementMatrix );

        Matrix< DDRMat > tElementRHS;
        this->get_element_rhs( Ii, tElementRHS );

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
    aMat->matrix_global_assembly();
    aVectorRHS->vector_global_asembly();
}


