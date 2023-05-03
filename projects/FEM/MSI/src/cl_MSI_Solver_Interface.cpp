/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Solver_Interface.cpp

 */
#include "cl_MSI_Solver_Interface.hpp"
#include "cl_MDL_Model.hpp"
#include "cl_MSI_Equation_Model.hpp"

#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Map.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Warehouse.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace MSI
    {

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::get_exact_solution_from_hdf5_and_calculate_error( const char* aFilename )
        {
            mPrevSolutionVector->read_vector_from_HDF5( aFilename );

            mSolutionVector->vec_plus_vec( 1.0, *mPrevSolutionVector, -1.0 );
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::get_residual_vector_for_output( const char* aFilename )
        {
            //        mPrevSolutionVector->read_vector_from_HDF5( aFilename );
            mPrevSolutionVector->vec_put_scalar( 5.0 );

            mSolutionVector->vec_plus_vec( 1.0, *mPrevSolutionVector, 0.0 );
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::write_solution_to_hdf5_file( const char* aFilename )
        {
            mSolutionVector->save_vector_to_HDF5( aFilename );
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::initiate_output(
                const uint aOutputIndex,
                const real aTime,
                const bool aEndOfTimeIteration )
        {
            // end of time iteration that the exodus file should be closed
            mModel->output_solution( aOutputIndex, aTime, aEndOfTimeIteration );
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            mSolutionVector = aSolutionVector;
            mMSI->mEquationModel->set_solution_vector( mSolutionVector );
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_eigen_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            mEigenSolutionVector = aSolutionVector;
            mMSI->mEquationModel->set_eigen_solution_vector( mEigenSolutionVector );
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::postmultiply_implicit_dQds()
        {
            mMSI->mEquationModel->compute_explicit_and_implicit_dQIdp();
            // mMSI->mEquationModel->compute_explicit_dQIdp();
            // mMSI->mEquationModel->compute_implicit_dQIdp();
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::compute_IQI()
        {
            mMSI->mEquationModel->compute_IQIs();
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_solution_vector_prev_time_step( sol::Dist_Vector* aSolutionVector )
        {
            mPrevSolutionVector = aSolutionVector;
            mMSI->mEquationModel->set_previous_solution_vector( mPrevSolutionVector );
        }

        //------------------------------------------------------------------------------

        sol::Dist_Vector*
        MSI_Solver_Interface::get_solution_vector_prev_time_step()
        {
            return mPrevSolutionVector;
        }

        //------------------------------------------------------------------------------

        sol::Dist_Vector*
        MSI_Solver_Interface::get_solution_vector(
                const moris::Cell< enum MSI::Dof_Type >& aListOfDofTypes,
                moris::Cell< moris_index > const &       aLocalCoefficientsIndices )
        {
            // create the factory based on the tpl
            sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

            // get the map that goes from local adof mesh index to adof id
            moris::Matrix< DDSMat > tLocalAdofIds = this->get_my_local_global_overlapping_map( aListOfDofTypes );

            // pick the indices that are provided in the input list
            moris::Matrix< DDSMat > tRequiredAdofIds( aLocalCoefficientsIndices.size(), 1 );

            // fill out the ids based on the provided local indices
            std::transform( aLocalCoefficientsIndices.begin(), aLocalCoefficientsIndices.end(), tRequiredAdofIds.begin(),    //
                    [ &tLocalAdofIds ]( moris_index const & aCoeffIndex ) { return tLocalAdofIds( aCoeffIndex ); } );

            // create a map based on the dof type that is specified
            sol::Dist_Map* tDofDMap = tMatFactory.create_map( tRequiredAdofIds );

            // create a dist. vector
            sol::Dist_Vector* tDofDVec = tMatFactory.create_vector( this, tDofDMap, 1 );

            // populate the vector based on the map
            tDofDVec->import_local_to_global( *mSolutionVector );

            // return the dist vector
            return tDofDVec;
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            mAdjointSolutionVector = aSolutionVector;
            mMSI->mEquationModel->set_adjoint_solution_vector( mAdjointSolutionVector );
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_previous_adjoint_solution_vector( sol::Dist_Vector* aSolutionVector )
        {
            mPreviousAdjointSolutionVector = aSolutionVector;
            mMSI->mEquationModel->set_previous_adjoint_solution_vector( mPreviousAdjointSolutionVector );
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_time( const Matrix< DDRMat >& aTime )
        {
            mTime = aTime;
            mMSI->mEquationModel->set_time( mTime );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        MSI_Solver_Interface::get_time()
        {
            return mTime;
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_previous_time( const Matrix< DDRMat >& aTime )
        {
            mPrevTime = aTime;
            mMSI->mEquationModel->set_previous_time( mPrevTime );
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat >
        MSI_Solver_Interface::get_previous_time()
        {
            return mPrevTime;
        }

        //------------------------------------------------------------------------------

        const moris::Cell< moris::Matrix< DDRMat > >&
        MSI_Solver_Interface::get_criteria( const moris::uint& aMySetInd )
        {
            return mMSI->get_equation_set( aMySetInd )->get_QI();
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_requested_IQI_names( const moris::Cell< std::string >& aIQINames )
        {
            mMSI->get_equation_model()->set_requested_IQI_names( aIQINames );
        }

        //-------------------------------------------------------------------------------------------------------

        moris::uint
        MSI_Solver_Interface::get_num_rhs()
        {
            return mMSI->get_equation_model()->get_num_rhs();
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::initialize_set(
                const uint                      aMyEquSetInd,
                const bool                      aIsStaggered,
                const fem::Time_Continuity_Flag aTimeContinuityOnlyFlag,
                const bool                      aIsAdjointOffDiagonalTimeContribution )
        {
            mMSI->get_equation_model()->set_is_adjoint_off_diagonal_time_contribution( aIsAdjointOffDiagonalTimeContribution );

            mMSI->get_equation_set( aMyEquSetInd )->initialize_set( aIsStaggered, aTimeContinuityOnlyFlag );
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::report_beginning_of_assembly()
        {
            mMSI->get_equation_model()->reset();
        }

        //------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::report_end_of_assembly()
        {
            mMSI->get_equation_model()->report_on_assembly();
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::free_block_memory( const uint aMyEquSetInd )
        {
            mMSI->get_equation_model()->set_is_adjoint_off_diagonal_time_contribution( false );

            mMSI->get_equation_set( aMyEquSetInd )->free_matrix_memory();
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::set_solver_warehouse( std::shared_ptr< sol::SOL_Warehouse > aSolverWarehouse )
        {
            mSolverWarehouse = aSolverWarehouse;
        }

        //-------------------------------------------------------------------------------------------------------

        void
        MSI_Solver_Interface::compute_sparsity_pattern()
        {
            // This is meant for the mapper solver since it is working with Trillions
            if ( !mSolverWarehouse ) return;

            /* ---------------------------------------------------------------------------------------- */
            /* Step 0: build the comm map  */

            // get the communication table and build a communication cell
            Matrix< IdMat > tCommTable = mDofMgn->get_comm_table();

            // convert it to a cell
            moris::Cell< moris_index > tCommCell;
            tCommCell.insert( tCommCell.size(), tCommTable.begin(), tCommTable.end() );

            // build a map between neighboring processors and index in the comm table
            std::map< moris_id, moris_index > tProcIdToCommTableIndex;
            for ( uint iLocalProcIndex = 0; iLocalProcIndex < tCommCell.size(); iLocalProcIndex++ )
            {
                tProcIdToCommTableIndex.insert( std::make_pair( tCommCell( iLocalProcIndex ), iLocalProcIndex ) );
            }

            /* ---------------------------------------------------------------------------------------- */
            /* Step 1: get the necessary information and build initial data to store the desired data structures  */

            // get the owned and owned/shared adof ids on this processor
            Matrix< DDSMat >      tOwnedAdofIndexToIdMap = this->get_my_local_global_map();
            Matrix< DDSMat >      tAdofIndexToIdMap      = this->get_my_local_global_overlapping_map();
            moris::Cell< Adof* >& tAdofs                 = mDofMgn->get_adofs();

            // Cache for reassembled entries on owned and shared rows
            // Outer cell is the the row number
            // Inner cell is the list of columns
            // This data will process in process connectivity
            moris::Cell< moris::Cell< uint > > tOwnedRows( tOwnedAdofIndexToIdMap.numel() );
            moris::Cell< moris::Cell< uint > > tSharedRows( tOwnedAdofIndexToIdMap.numel() );

            // this is for cross-processor coupling
            moris::size_t tNumSharedRows = tAdofs.size() - tOwnedAdofIndexToIdMap.numel();


            // list of all shared adof ids in this processor
            moris::Cell< moris::Cell< uint > > tSharedRowsIds( tCommCell.size(), moris::Cell< uint >( tNumSharedRows ) );

            // initialize the shared rows off processors
            // the outer cell is the processor rank
            // the inner cells are lits of adofs that are created
            moris::Cell< moris::Cell< moris::Cell< uint > > > tSharedRowsOffProc( tCommCell.size(), moris::Cell< moris::Cell< uint > >( tNumSharedRows ) );


            /* ---------------------------------------------------------------------------------------- */
            /* Step 2: build the transpose of the maps and fill out the data strcuture needed to process  */

            // adof id to adof local index map
            std::unordered_map< moris_id, moris_index > tOwnedAofIdToIndexMap;

            // shared adof id to {local index(0-based), neighboring processor index}
            std::unordered_map< moris_id, std::pair< moris_index, moris_id > > tSharedAofIdToIndexMap;

            // reserve space in the map
            tOwnedAofIdToIndexMap.reserve( tOwnedAdofIndexToIdMap.numel() );
            tSharedAofIdToIndexMap.reserve( tNumSharedRows );

            // transpose map for the owned adofs , from adof id to local adof index
            for ( uint iAdofIndex = 0; iAdofIndex < tOwnedAdofIndexToIdMap.numel(); iAdofIndex++ )
            {
                tOwnedAofIdToIndexMap[ tOwnedAdofIndexToIdMap( iAdofIndex ) ] = iAdofIndex;
            }

            // transpose map for all of shared adofs on this processor
            // counter for local adof shared index
            uint iCounter = 0;

            for ( uint iAdofIndex = 0; iAdofIndex < tAdofIndexToIdMap.numel(); iAdofIndex++ )
            {
                // if it does not belong to owned then it must be shared
                if ( tOwnedAofIdToIndexMap.find( tAdofIndexToIdMap( iAdofIndex ) ) == tOwnedAofIdToIndexMap.end() )
                {
                    // get the owning processor
                    moris_id tOwner = tAdofs( iAdofIndex )->get_adof_owning_processor();

                    // check if it belongs to a neighboring processor
                    auto tIter = tProcIdToCommTableIndex.find( tOwner );
                    MORIS_ASSERT(
                            tIter != tProcIdToCommTableIndex.end(),
                            "MSI_Solver_Interface::compute_sparsity_pattern() - "
                            "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                            tOwner,
                            par_rank(),
                            ios::stringify_log( tCommTable ).c_str() );
                    moris_index tProcDataIndex = tIter->second;

                    // add it to the list of shared adof ids in the correct neighboring processor
                    tSharedRowsIds( tProcDataIndex )( iCounter ) = tAdofIndexToIdMap( iAdofIndex );

                    // add to shared adof map
                    tSharedAofIdToIndexMap[ tAdofIndexToIdMap( iAdofIndex ) ] = { iCounter, tProcDataIndex };

                    // increase the local index of the shared adof
                    iCounter++;
                }
            }


            // get initial estimate
            uint tNumNonzero = this->estimate_number_of_nonzero_columns();

            //  reserve out the data related to on-processor based data
            for ( uint iRowNum = 0; iRowNum < tOwnedRows.size(); iRowNum++ )
            {
                tOwnedRows( iRowNum ).reserve( tNumNonzero * 2 );
                tSharedRows( iRowNum ).reserve( tNumNonzero );
            }

            //  reserve out the data related to off-processor based data
            for ( auto& iAofListsOnProc : tSharedRowsOffProc )
            {
                for ( auto& iAodList : iAofListsOnProc )
                {
                    iAodList.reserve( tNumNonzero * 2 );
                }
            }

            /* ---------------------------------------------------------------------------------------- */
            /* Step 3: fill the on-processor part of the sparsity pattern and prepare the data for the off-processor  */

            // get number of blocks
            moris::uint tNumEqBlocks = this->get_num_my_blocks();

            // loop over number of blocks to the get the element topology (adof connectivity) and generate sparsity pattern
            for ( size_t iEqBlock = 0; iEqBlock < tNumEqBlocks; iEqBlock++ )
            {
                // get number of equation objects(lagrange IP elements)
                moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( iEqBlock );

                // loop over the lagrange IP elements to find the basis that are interpolating into them
                for ( moris::uint iEqObject = 0; iEqObject < tNumEquationObjectOnSet; iEqObject++ )
                {
                    // get the element topology
                    Matrix< DDSMat > tElementTopology;
                    this->get_element_topology( iEqBlock, iEqObject, tElementTopology );

                    // loop over the global adof ids and cross-couple them toghter
                    for ( uint i = 0; i < tElementTopology.numel(); i++ )
                    {
                        // get ownership information of adof
                        auto tIterator = tOwnedAofIdToIndexMap.find( tElementTopology( i ) );

                        // decide based on the ownership of adofs where to put it
                        //  the if blocks pertain to the on-processor part of the data
                        if ( tIterator != tOwnedAofIdToIndexMap.end() )
                        {
                            // sort through element topology and decide what adofs are in the diagonal block and what are in the off-diagonal block
                            std::partition_copy( tElementTopology.begin(), tElementTopology.end(),    //
                                    std::back_inserter( tOwnedRows( tIterator->second ) ),
                                    std::back_inserter( tSharedRows( tIterator->second ) ),           //
                                    [ &tOwnedAofIdToIndexMap ]( int aLocalAodfId ) {                  //
                                        // lambda function to check the the membership of the adof
                                        auto tIterator = tOwnedAofIdToIndexMap.find( aLocalAodfId );
                                        return tIterator != tOwnedAofIdToIndexMap.end();
                                    } );
                        }
                        // Belong to another processor, fill out the data to be communicated later
                        else
                        {
                            // get index, ownership pair tp decide what processor to put the data and what index
                            auto tIndexOwnershipPair = tSharedAofIdToIndexMap.find( tElementTopology( i ) );

                            // get the relevant part of the data that will be communictaed based on the ownership and index
                            moris::Cell< uint >& tConnectedAdofs = tSharedRowsOffProc( tIndexOwnershipPair->second.second )( tIndexOwnershipPair->second.first );

                            // put all the connected adofs without deciding on-off diagonal
                            // this data will be processed by the owning processor of this adof
                            tConnectedAdofs.insert( tConnectedAdofs.size(), tElementTopology.begin(), tElementTopology.end() );
                        }
                    }
                }
            }

            /* ---------------------------------------------------------------------------------------- */
            /* Step 4: communicate the shared adof data to the owning processor  */

            // outer cell: received from neighbour processor index, inner cells: list of all adofs that are connoted to shared adofs consecutively
            moris::Cell< moris::Cell< uint > > tAdofConnectivityReceive;

            // outer cell: received from neighbour processor index, inner cells: offset indicating position of the connected adofs
            moris::Cell< moris::Cell< uint > > tAdofConnectivityOffsetReceive;

            this->communicate_shared_adof_connectivity( tSharedRowsOffProc, tAdofConnectivityReceive, tAdofConnectivityOffsetReceive, tCommCell );

            // communicate shared row ids to all neighboring processors
            moris::Cell< moris::Cell< uint > > tSharedRowsIdsReceive;
            communicate_cells( tCommCell, tSharedRowsIds, tSharedRowsIdsReceive );

            /* ---------------------------------------------------------------------------------------- */
            /* Step 5: analyze the communicated data based and put in the on and off diagonal parts  */

            // local neighboring processor index
            uint iProc = 0;

            // loop over the adof connectivity based on the processor index
            for ( auto const & iAdofConnOffset : tAdofConnectivityOffsetReceive )
            {

                // get the list of all connected dofs got from this processor
                moris::Cell< uint >& tAdofConn = tAdofConnectivityReceive( iProc );

                // check if there are any connected adofs to this adof
                if ( iAdofConnOffset.size() < 2 )
                {
                    iProc++;
                    continue;
                }

                // if there is connectivity data loop over the offset data
                for ( uint i = 0; i < iAdofConnOffset.size() - 1; i++ )
                {

                    // if there are connected adofs, if there are adofs connected to this ad
                    if ( ( iAdofConnOffset( i + 1 ) - iAdofConnOffset( i ) ) != 0 )
                    {
                        // the global adof id
                        uint tAdofGlobalId = tSharedRowsIdsReceive( iProc )( i );

                        // find the corresponding index of the adof
                        auto tIterator = tOwnedAofIdToIndexMap.find( tAdofGlobalId );

                        // reserve enough space on the connected adofs
                        tOwnedRows( tIterator->second ).reserve( tOwnedRows.size() + iAdofConnOffset( i + 1 ) - iAdofConnOffset( i ) );
                        tSharedRows( tIterator->second ).reserve( tSharedRows.size() + iAdofConnOffset( i + 1 ) - iAdofConnOffset( i ) );

                        // grab the part in the connectivity adof list and parrtion it based on the diagonal and off-digonal data
                        std::partition_copy( tAdofConn.begin() + iAdofConnOffset( i ), tAdofConn.begin() + iAdofConnOffset( i + 1 ),    //
                                std::back_inserter( tOwnedRows( tIterator->second ) ),
                                std::back_inserter( tSharedRows( tIterator->second ) ),                                                 //
                                [ &tOwnedAofIdToIndexMap ]( int aLocalAodfId ) {                                                        //
                                    auto tIterator = tOwnedAofIdToIndexMap.find( aLocalAodfId );
                                    return tIterator != tOwnedAofIdToIndexMap.end();
                                } );
                    }
                }

                // increment the processor count
                iProc++;
            }

            /* ---------------------------------------------------------------------------------------- */
            /* Step 6: determine the sparsity pattern with the on-processor + analyzed off-processor data  */

            // set the size of data that belong to DLA solver interface
            mNonZeroDigonal.resize( tOwnedRows.size() );
            mNonZeroOffDigonal.resize( tOwnedRows.size() );

            // loop over the owned rows
            for ( uint iRowNum = 0; iRowNum < tOwnedRows.size(); iRowNum++ )
            {
                // sort in place
                std::sort( tOwnedRows( iRowNum ).begin(), tOwnedRows( iRowNum ).end() );
                std::sort( tSharedRows( iRowNum ).begin(), tSharedRows( iRowNum ).end() );

                // find the first unique
                auto tOwnedUnique  = std::unique( tOwnedRows( iRowNum ).begin(), tOwnedRows( iRowNum ).end() );
                auto tSharedUnique = std::unique( tSharedRows( iRowNum ).begin(), tSharedRows( iRowNum ).end() );

                // find the number of non-zeros
                mNonZeroDigonal( iRowNum )    = std::distance( tOwnedRows( iRowNum ).begin(), tOwnedUnique );
                mNonZeroOffDigonal( iRowNum ) = std::distance( tSharedRows( iRowNum ).begin(), tSharedUnique );
            }
        }

        //-------------------------------------------------------------------------------------------------------

        uint
        MSI_Solver_Interface::estimate_number_of_nonzero_columns()
        {
            // get spatial dimension from the model
            uint tSpatialDim = mModel->get_spatial_dim();

            uint tNumNonZero = 0;
            // get number of blocks
            moris::uint tNumEqBlocks = this->get_num_my_elements();

            // loop over number of blocks
            for ( size_t iEqBlock = 0; iEqBlock < tNumEqBlocks; iEqBlock++ )
            {
                // get number of equation objects
                moris::uint tNumEquationObjectOnSet = this->get_num_equation_objects_on_set( iEqBlock );

                // get the first equation object
                for ( moris::uint iEqObject = 0; iEqObject < tNumEquationObjectOnSet; iEqObject++ )
                {
                    // obtain element topology
                    Matrix< DDSMat > tElementTopology;
                    this->get_element_topology( iEqBlock, iEqObject, tElementTopology );

                    // get connectivity
                    tNumNonZero = tElementTopology.numel();

                    // quit as soon as a non-empty element is found
                    goto tReturnStatement;
                }
            }

        // go to tge return statement
        tReturnStatement:
            return tNumNonZero * tSpatialDim;
        }


        void
        MSI_Solver_Interface::communicate_shared_adof_connectivity(
                moris::Cell< moris::Cell< moris::Cell< uint > > > const & aSharedAdofConn,
                moris::Cell< moris::Cell< uint > >&                       aAdofConnectivityReceive,
                moris::Cell< moris::Cell< uint > >&                       aAdofConnectivityOffsetReceive,
                moris::Cell< moris_index > const &                        aCommCell )
        {
            // convert 3 times nested cell to 2 cells (data cell + offset cell)
            moris::Cell< moris::Cell< uint > > tAdofConnectivityOffsetSend( aSharedAdofConn.size() );

            // convert each cell of cells to one cell with connectivity and count the size of the data cell
            for ( uint iProc = 0; iProc < tAdofConnectivityOffsetSend.size(); iProc++ )
            {
                // resize with the proper number
                tAdofConnectivityOffsetSend( iProc ).resize( aSharedAdofConn( iProc ).size() + 1, 0 );

                for ( size_t idAofId = 0; idAofId < aSharedAdofConn( iProc ).size(); idAofId++ )
                {
                    tAdofConnectivityOffsetSend( iProc )( idAofId + 1 ) = tAdofConnectivityOffsetSend( iProc )( idAofId ) + aSharedAdofConn( iProc )( idAofId ).size();
                }
            }

            // prepare the data cell
            moris::Cell< moris::Cell< uint > > tAdofConnectivitySend( aSharedAdofConn.size() );

            // loop over the data cell sizes
            for ( uint iProc = 0; iProc < tAdofConnectivitySend.size(); iProc++ )
            {
                // get the size of data cell from the the offset , the last element
                tAdofConnectivitySend( iProc ).reserve( *( tAdofConnectivityOffsetSend( iProc ).end() - 1 ) );

                // accumulate the data cell by inserting all the values in it for the relevant processor
                for ( size_t idAofId = 0; idAofId < aSharedAdofConn( iProc ).size(); idAofId++ )
                {
                    tAdofConnectivitySend( iProc ).insert( tAdofConnectivitySend( iProc ).size(), aSharedAdofConn( iProc )( idAofId ).begin(), aSharedAdofConn( iProc )( idAofId ).end() );
                }
            }

            // let all the processors reach this point
            barrier();

            // communicate cells
            communicate_cells( aCommCell, tAdofConnectivitySend, aAdofConnectivityReceive );
            communicate_cells( aCommCell, tAdofConnectivityOffsetSend, aAdofConnectivityOffsetReceive );
        }

    }    // namespace MSI
}    // namespace moris
