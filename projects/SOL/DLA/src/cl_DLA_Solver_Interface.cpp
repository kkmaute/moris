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
void Solver_Interface::build_graph( moris::sol::Dist_Matrix * aMat )
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

    // global assembly to communicate entries
    aMat->matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------
void Solver_Interface::fill_matrix_and_RHS(
        moris::sol::Dist_Matrix * aMat,
        moris::sol::Dist_Vector * aVectorRHS,
        moris::sol::Dist_Vector * aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numLocElements = this->get_num_my_elements();

    // Loop over all local elements to fill matrix and RHS
    for (moris::uint Ii=0; Ii< numLocElements; Ii++)
    {
        moris::Matrix< DDSMat > tElementTopology;
        this->get_element_topology( Ii, tElementTopology );

        Matrix< DDRMat >         tElementMatrix;
        Cell< Matrix< DDRMat > > tElementRHS;
        this->get_equation_object_operator_and_rhs( Ii, tElementMatrix, tElementRHS );

        // Fill element in distributed matrix
        aMat->fill_matrix(
                tElementTopology.length(),
                tElementMatrix,
                tElementTopology );

        // Fill elementRHS in distributed RHS
        aVectorRHS->sum_into_global_values(
                tElementTopology,
                tElementRHS(0) );
    }

    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_asembly();
    aMat->matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_RHS( moris::sol::Dist_Vector * aVectorRHS )
{
    // Get local number of elements
    moris::uint tNumBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < tNumBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false );                     // FIXME FIXME shoudl be true. this is a brutal hack and will be changed in a few days

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

            Cell< Matrix< DDRMat > > tElementRHS;
            this->get_equation_object_rhs( Ii, Ik, tElementRHS );

            for ( moris::uint Ia=0; Ia < tNumRHS; Ia++ )
            {
                // Fill elementRHS in distributed RHS
                aVectorRHS->sum_into_global_values(
                        tElementTopology,
                        tElementRHS( Ia ),
                        Ia );
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_asembly();

//    aVectorRHS->print();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_additional_DqDs_RHS_contribution( moris::sol::Dist_Vector * aVectorRHS )
{
    if( !this->get_is_forward_analysis() )
    {
        // Get local number of elements
        moris::uint tNumBlocks = this->get_num_my_blocks();

        moris::uint tNumRHS = this->get_num_rhs();

        // Loop over all local elements to build matrix graph
        for ( moris::uint Ii=0; Ii < tNumBlocks; Ii++ )
        {
            moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

            this->initialize_set( Ii, false, true);                     // FIXME FIXME shoudl be true. this is a brutal hack and will be changed in a few days

            for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
            {
                Matrix< DDSMat > tElementTopology;
                this->get_element_topology(Ii, Ik, tElementTopology );

                Cell< Matrix< DDRMat > > tElementRHS;
                this->get_equation_object_off_diag_rhs( Ii, Ik, tElementRHS );

                for ( moris::uint Ia=0; Ia < tNumRHS; Ia++ )
                {
                    // Fill elementRHS in distributed RHS
                    aVectorRHS->sum_into_global_values(
                            tElementTopology,
                            tElementRHS( Ia ),
                            Ia );
                }
            }

            this->free_block_memory( Ii );
        }

        // global assembly to switch entries to the right processor
        aVectorRHS->vector_global_asembly();

//        aVectorRHS->print();
    }
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_jacobian( moris::sol::Dist_Matrix * aMat )
{
    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false );

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

            Matrix< DDRMat > tElementMatrix;
            this->get_equation_object_operator( Ii, Ik, tElementMatrix );

            // Fill element in distributed matrix
            aMat->fill_matrix( tElementTopology.length(),
                    tElementMatrix,
                    tElementTopology );
        }
        aMat->matrix_global_assembly();   // Mathias: check whether this needs to be done for each block
        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aMat->matrix_global_assembly();

    //aMat->print();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::fill_matrix_and_RHS( moris::sol::Dist_Matrix * aMat,
        moris::sol::Dist_Vector * aVectorRHS )
{
    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii=0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false );

        for ( moris::uint Ik=0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology(Ii, Ik, tElementTopology );

            Matrix< DDRMat >         tElementMatrix;
            Cell< Matrix< DDRMat > > tElementRHS;
            this->get_equation_object_operator_and_rhs( Ii, Ik, tElementMatrix, tElementRHS );

            // Fill element in distributed matrix
            aMat->fill_matrix(
                    tElementTopology.length(),
                    tElementMatrix,
                    tElementTopology );

            // Loop over all RHS vectors
            for ( moris::uint Ia=0; Ia < tNumRHS; Ia++ )
            {
                // Fill elementRHS in distributed RHS
                aVectorRHS->sum_into_global_values(
                        tElementTopology,
                        tElementRHS( Ia ),
                        Ia );
            }
        }
        aMat->matrix_global_assembly();
        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aMat->matrix_global_assembly();
    //    aMat->print();
	//aVectorRHS->print();

    aVectorRHS->vector_global_asembly();
}

void Solver_Interface::get_adof_ids_based_on_criteria( 
        moris::Cell< moris::Matrix< IdMat > > & aCriteriaIds,
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

            // hardcoded for right now. FIXME the moment I have time
            this->set_requested_IQI_names( {"IQIBulkVolumeFraction"} );

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
