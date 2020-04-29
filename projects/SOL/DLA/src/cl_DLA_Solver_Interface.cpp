/*
 * cl_DLA_Solver_Interface.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: schmidt */

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

using namespace moris;

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::build_graph( moris::Dist_Matrix * aMat )
{
    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
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
void Solver_Interface::fill_matrix_and_RHS( moris::Dist_Matrix * aMat,
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
        this->get_equation_object_operator( Ii, tElementMatrix );

        Cell< Matrix< DDRMat > > tElementRHS;
        this->get_equation_object_rhs( Ii, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix( tElementTopology.length(),
                           tElementMatrix,
                           tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values( tElementTopology,
                                            tElementRHS(0) );
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
    moris::uint tNumBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < tNumBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, true );

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

//            print(tElementTopology,"tElementTopology");

            Cell< Matrix< DDRMat > > tElementRHS;
            this->get_equation_object_rhs( Ii, Ik, tElementRHS );

//            print(tElementRHS,"tElementRHS");

            for ( moris::uint Ia=0; Ia < tNumRHS; Ia++ )
            {
                // Fill elementRHS in distributed RHS
                aVectorRHS->sum_into_global_values( tElementTopology,
                                                    tElementRHS( Ia ),
                                                    Ia );
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right proceccor
    aVectorRHS->vector_global_asembly();

//    std::cout<<"Assembled Residual Vector"<<std::endl;
//    aVectorRHS->print();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::assemble_jacobian( moris::Dist_Matrix * aMat,
                                          moris::Dist_Vector   * aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
//        std::cout<<"Block "<<Ii<<std::endl;
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false );

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

//            print(tElementTopology,"tElementTopology");

            Matrix< DDRMat > tElementMatrix;
            this->get_equation_object_operator( Ii, Ik, tElementMatrix );

//            print(tElementMatrix,"tElementMatrix");

            // Fill element in distributed matrix
            aMat->fill_matrix( tElementTopology.length(),
                               tElementMatrix,
                               tElementTopology );
        }
        aMat->matrix_global_assembly();
        this->free_block_memory( Ii );
    }
    // global assembly to switch entries to the right processor
    aMat->matrix_global_assembly();

//    aMat->save_matrix_to_matlab_file( "Matrix.dat");
//    std::cout<<"-----------------Matrix output finished---------------------"<<std::endl;

//    aMat->print();

}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::fill_matrix_and_RHS( moris::Dist_Matrix * aMat,
                                            moris::Dist_Vector * aVectorRHS )
{
    // Get local number of elements
    moris::uint numLocElements = this->get_num_my_elements();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        this->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat > tElementMatrix;
        this->get_equation_object_operator( Ii, tElementMatrix );

        Cell< Matrix< DDRMat > >tElementRHS;
        this->get_equation_object_rhs( Ii, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix( tElementTopology.length(),
                           tElementMatrix,
                           tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values( tElementTopology,
                                            tElementRHS(0) );
    }

    // global assembly to switch entries to the right proceccor
    aMat->matrix_global_assembly();
    aVectorRHS->vector_global_asembly();
}

void Solver_Interface::get_adof_ids_based_on_criteria( moris::Cell< moris::Matrix< IdMat > > & aCriteriaIds,
                                                       const moris::real                       aThreshold  )       // FIXME find better name
{
    // Get number of Sets
    moris::uint tNumSets = this->get_num_my_blocks();

    uint tCounter = 0;
    moris::real tMinVolVraction = 1.0;

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < tNumSets; Ii++ )
    {
    	// only check bulk sets
        if( this->get_set_type( Ii ) == fem::Element_Type::BULK )
        {
            // get number of equations on set
            moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

            // resize adof id vector
            aCriteriaIds.resize( aCriteriaIds.size() + tNumEquationObjectOnSet );

            moris::Cell< moris::Cell< enum fem::IQI_Type > > tRequestedIQITypes( 1 );
            tRequestedIQITypes( 0 ).resize( 1, fem::IQI_Type::VOLUME_FRACTION );

            this->set_requested_IQI_type( Ii, tRequestedIQITypes );

            // initialize set
            this->initialize_set( Ii, false );

            // loop over equation objects on set
            for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
            {
                // calcuate criteria
                this->calculate_criteria( Ii, Ik );

                // get criteria
                const moris::Cell< moris::Matrix< DDRMat> > & tCriteria = this->get_criteria( Ii );

                // if criteria mets requirement
                if( tCriteria( 0 )( 0 ) < aThreshold )
                {
                    // get adof ids of this equation object
                    moris::Matrix< DDSMat > tMat;
                    this->get_element_topology( Ii, Ik, tMat );

                    // store ids in cell
                    aCriteriaIds( tCounter++ ) = tMat;
                }

                tMinVolVraction = std::min( tMinVolVraction, tCriteria( 0 )( 0 ) );
            }
        }

        // resize cell to number of triggered equation objects
        aCriteriaIds.resize( tCounter );
    }
}


