/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Solver_Interface.cpp

*/

#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Warehouse.hpp"

using namespace moris;

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::build_graph( moris::sol::Dist_Matrix* aMat, bool aUseSparsityPattern )
{
    if ( aUseSparsityPattern )
    {
        // compute the sparsity pattern
        this->compute_sparsity_pattern();

        // build the sparsity pattern based on the
        aMat->build_graph( mNonZeroDigonal, mNonZeroOffDigonal );

        return;
    }

    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology( Ii, Ik, tElementTopology );

            aMat->build_graph( tElementTopology.length(), tElementTopology );
        }
    }

    // global assembly to communicate entries
    aMat->initial_matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::fill_matrix_and_RHS(
        moris::sol::Dist_Matrix* aMat,
        moris::sol::Dist_Vector* aVectorRHS,
        moris::sol::Dist_Vector* aFullSolutionVector )
{
    this->set_solution_vector( aFullSolutionVector );

    // Get local number of elements
    moris::uint numLocElements = this->get_num_my_elements();

    moris::Matrix< DDSMat >    tElementTopology;
    Matrix< DDRMat >           tElementMatrix;
    Vector< Matrix< DDRMat > > tElementRHS;

    // Loop over all local elements to fill matrix and RHS
    for ( moris::uint Ii = 0; Ii < numLocElements; Ii++ )
    {
        this->get_element_topology( Ii, tElementTopology );

        // compute linear operator and RHS
        this->get_equation_object_operator_and_rhs( Ii, tElementMatrix, tElementRHS );

        // Fill element in distributed matrix
        if ( tElementMatrix.numel() > 0 )
        {
            aMat->fill_matrix(
                    tElementTopology.length(),
                    tElementMatrix,
                    tElementTopology );
        }

        // Fill elementRHS in distributed RHS
        if ( tElementRHS.size() > 0 )
        {
            if ( tElementRHS( 0 ).numel() > 0 )
            {
                aVectorRHS->sum_into_global_values(
                        tElementTopology,
                        tElementRHS( 0 ) );
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_assembly();
    aMat->matrix_global_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_RHS(
        moris::sol::Dist_Vector*        aVectorRHS,
        const fem::Time_Continuity_Flag aTimeContinuityOnlyFlag )
{
    // Get local number of elements
    moris::uint tNumBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    Vector< Matrix< DDRMat > > tElementRHS;
    Matrix< DDSMat >           tElementTopology;

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < tNumBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false, aTimeContinuityOnlyFlag );

        for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            this->get_element_topology( Ii, Ik, tElementTopology );

            // compute RHS
            this->get_equation_object_rhs( Ii, Ik, tElementRHS );

            // Fill elementRHS in distributed RHS
            if ( tElementRHS.size() > 0 )
            {
                MORIS_ASSERT( tElementRHS.size() == tNumRHS,
                        "Number of RHS does not match cell with RHS vectors.\n" );

                for ( moris::uint Ia = 0; Ia < tNumRHS; Ia++ )
                {
                    if ( tElementRHS( Ia ).numel() > 0 )
                    {
                        aVectorRHS->sum_into_global_values(
                                tElementTopology,
                                tElementRHS( Ia ),
                                Ia );
                    }
                }
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_staggered_RHS_contribution( moris::sol::Dist_Vector* aVectorRHS )
{
    // Get local number of elements
    moris::uint tNumBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    Matrix< DDSMat >           tElementTopology;
    Vector< Matrix< DDRMat > > tElementRHS;

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < tNumBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, true );    // FIXME FIXME should be true. this is a brutal hack and will be changed in a few days

        for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            this->get_element_topology( Ii, Ik, tElementTopology );

            this->get_equation_object_staggered_rhs( Ii, Ik, tElementRHS );

            // Fill elementRHS in distributed RHS
            if ( tElementRHS.size() > 0 )
            {
                MORIS_ASSERT( tElementRHS.size() == tNumRHS,
                        "Number of RHS does not match cell with RHS vectors.\n" );

                for ( moris::uint Ia = 0; Ia < tNumRHS; Ia++ )
                {
                    if ( tElementRHS( Ia ).numel() > 0 )
                    {
                        aVectorRHS->sum_into_global_values(
                                tElementTopology,
                                tElementRHS( Ia ),
                                Ia );
                    }
                }
            }
        }

        this->free_block_memory( Ii );
    }
    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_additional_DqDs_RHS_contribution( moris::sol::Dist_Vector* aVectorRHS )
{
    // Get local number of elements
    moris::uint tNumBlocks = this->get_num_my_blocks();

    moris::uint tNumRHS = this->get_num_rhs();

    Matrix< DDSMat >           tElementTopology;
    Vector< Matrix< DDRMat > > tElementRHS;

    this->report_beginning_of_assembly();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < tNumBlocks; Ii++ )
    {
        // only check bulk sets
        if ( this->get_set_type( Ii ) == fem::Element_Type::TIME_SIDESET )
        {
            moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

            this->initialize_set( Ii, false, fem::Time_Continuity_Flag::DEFAULT, true );

            for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
            {
                this->get_element_topology( Ii, Ik, tElementTopology );

                this->get_equation_object_off_diag_rhs( Ii, Ik, tElementRHS );

                // Fill elementRHS in distributed RHS
                if ( tElementRHS.size() > 0 )
                {
                    MORIS_ASSERT( tElementRHS.size() == tNumRHS,
                            "Number of RHS does not match cell with RHS vectors.\n" );

                    for ( moris::uint Ia = 0; Ia < tNumRHS; Ia++ )
                    {
                        if ( tElementRHS( Ia ).numel() > 0 )
                            aVectorRHS->sum_into_global_values(
                                    tElementTopology,
                                    tElementRHS( Ia ),
                                    Ia );
                    }
                }
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aVectorRHS->vector_global_assembly();

    this->report_end_of_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::assemble_jacobian(
        moris::sol::Dist_Matrix*        aMat,
        const fem::Time_Continuity_Flag aTimeContinuityOnlyFlag )
{
    // Get local number of elements
    moris::uint numBlocks = this->get_num_my_blocks();

    Matrix< DDSMat > tElementTopology;
    Matrix< DDRMat > tElementMatrix;

    this->report_beginning_of_assembly();

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < numBlocks; Ii++ )
    {
        moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

        this->initialize_set( Ii, false, aTimeContinuityOnlyFlag );

        for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
        {
            this->get_element_topology( Ii, Ik, tElementTopology );

            this->get_equation_object_operator( Ii, Ik, tElementMatrix );

            // Fill element in distributed matrix
            if ( tElementMatrix.numel() > 0 )
            {
                aMat->fill_matrix(
                        tElementTopology.length(),
                        tElementMatrix,
                        tElementTopology );
            }
        }

        this->free_block_memory( Ii );
    }

    // global assembly to switch entries to the right processor
    aMat->matrix_global_assembly();

    this->report_end_of_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::fill_matrix_and_RHS(
        moris::sol::Dist_Matrix* aMat,
        moris::sol::Dist_Vector* aVectorRHS )
{
    this->report_beginning_of_assembly();

    uint const tNumBlocks = this->get_num_my_blocks();
    uint const tNumRHS    = this->get_num_rhs();
    // Loop over all local elements to build matrix graph
    for ( uint iBlock = 0; iBlock < tNumBlocks; iBlock++ )
    {
        this->initialize_set( iBlock );

        uint const tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( iBlock );

        for ( moris::uint iEquationObject = 0; iEquationObject < tNumEquationObjectOnSet; iEquationObject++ )
        {
            Matrix< DDSMat > tElementTopology;
            this->get_element_topology( iBlock, iEquationObject, tElementTopology );

            Matrix< DDRMat >           tElementMatrix;
            Vector< Matrix< DDRMat > > tElementRHS;
            this->get_equation_object_operator_and_rhs( iBlock, iEquationObject, tElementMatrix, tElementRHS );

            // Fill element in distributed matrix
            if ( tElementMatrix.numel() > 0 )
            {
                aMat->fill_matrix(
                        tElementTopology.length(),
                        tElementMatrix,
                        tElementTopology );
            }

            // Loop over all RHS vectors
            if ( tElementRHS.size() > 0 )
            {
                MORIS_ASSERT( tElementRHS.size() == tNumRHS,
                        "Number of RHS does not match cell with RHS vectors.\n" );

                for ( moris::uint Ia = 0; Ia < tNumRHS; Ia++ )
                {
                    if ( tElementRHS( Ia ).numel() > 0 )
                    {
                        // Fill elementRHS in distributed RHS
                        aVectorRHS->sum_into_global_values(
                                tElementTopology,
                                tElementRHS( Ia ),
                                Ia );
                    }
                }
            }
        }

        this->free_block_memory( iBlock );
    }

    // global assembly to switch entries to the right processor
    aMat->matrix_global_assembly();
    aVectorRHS->vector_global_assembly();

    this->report_end_of_assembly();
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::get_adof_ids_based_on_criteria(
        Vector< moris::Matrix< IdMat > >& aCriteriaIds,
        const moris::real                        aThreshold )    // FIXME find better name
{
    // Get number of Sets
    moris::uint tNumSets = this->get_num_my_blocks();

    uint        tCounter        = 0;
    moris::real tMinVolFraction = 1.0;

    // Loop over all local elements to build matrix graph
    for ( moris::uint Ii = 0; Ii < tNumSets; Ii++ )
    {
        // only check bulk sets
        if ( this->get_set_type( Ii ) == fem::Element_Type::BULK )
        {
            // get number of equations on set
            moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( Ii );

            // resize adof id vector
            aCriteriaIds.resize( aCriteriaIds.size() + tNumEquationObjectOnSet );

            // hard-coded for right now. FIXME the moment I have time
            this->set_requested_IQI_names( { "IQIBulkVolumeFraction" } );

            // initialize set
            this->initialize_set( Ii );

            // loop over equation objects on set
            for ( moris::uint Ik = 0; Ik < tNumEquationObjectOnSet; Ik++ )
            {
                // compute criteria
                this->calculate_criteria( Ii, Ik );

                // get criteria
                const moris::Vector< moris::Matrix< DDRMat > >& tCriteria = this->get_criteria( Ii );

                // if criteria meets requirement
                if ( tCriteria( 0 )( 0 ) < aThreshold )
                {
                    // get adof ids of this equation object
                    moris::Matrix< DDSMat > tMat;
                    this->get_element_topology( Ii, Ik, tMat );

                    // store ids in cell
                    aCriteriaIds( tCounter++ ) = tMat;
                }

                tMinVolFraction = std::min( tMinVolFraction, tCriteria( 0 )( 0 ) );
            }
        }

        // resize cell to number of triggered equation objects
        aCriteriaIds.resize( tCounter );
    }
}

//---------------------------------------------------------------------------------------------------------

void Solver_Interface::set_solver_warehouse( std::shared_ptr< sol::SOL_Warehouse > aSolverWarehouse )
{
    MORIS_ERROR( false, "Solver_Interface::set_solver_warehouse(), not implemented for base class" );
}
