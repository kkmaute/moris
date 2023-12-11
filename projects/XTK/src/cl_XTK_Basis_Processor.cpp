/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Basis_Processor.cpp
 *
 */

#include "cl_XTK_Basis_Processor.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Subphase_Group.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "fn_XTK_convert_cell_to_map.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "fn_stringify_matrix.hpp"


#include "cl_XTK_HMR_Helper.hpp"    // hmr helper

using namespace moris;

namespace xtk
{
    Basis_Processor::Basis_Processor( xtk::Model* aXTKModelPtr )
            : mXTKModelPtr( aXTKModelPtr )
            , mMeshIndices( aXTKModelPtr->mEnrichment->mMeshIndices )
            , mParameterList( &aXTKModelPtr->mParameterList )
            , mSpatialDim( mXTKModelPtr->get_spatial_dim() )

    {
        // resize the basis data with the correct size
        mBasisData.resize( mMeshIndices.numel() );
        mHMRHelper.reserve( mMeshIndices.numel() );

        // create hmr helper objects with std::back_inserter and the correct mesh index
        std::transform( mMeshIndices.begin(), mMeshIndices.end(), std::back_inserter( mHMRHelper ), [ &aXTKModelPtr ]( moris_index iMeshIndex )    //
                { return new HMR_Helper( aXTKModelPtr, iMeshIndex ); } );
    }

    // ----------------------------------------------------------------------------------

    Basis_Processor::~Basis_Processor()
    {
        for ( auto& iHMRHelper : mHMRHelper )
        {
            delete iHMRHelper;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_basis_in_subspace( enum Dependency_Criteria aCriteria )
    {
        // choose the dependency criteria for blocking
        switch ( aCriteria )
        {
            case Dependency_Criteria::VOLUME:
            {
                this->construct_volumetric_clustering_of_basis();
                break;
            }


            default:
            {
                MORIS_ERROR( 0, "Basis_Processor::construct_basis_in_subspace, the given criteria does not exist" );
                break;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_basis_list( enum Dependency_Criteria aCriteria, moris_index aMeshIndex )
    {
        switch ( aCriteria )
        {
            case Dependency_Criteria::VOLUME:
            {
                this->construct_follower_basis_using_volume( aMeshIndex );
                break;
            }

            case Dependency_Criteria::BASIS_SUPPORT:
            {
                this->construct_follower_basis_using_basis_support( aMeshIndex );
                break;
            }

            default:
            {
                MORIS_ERROR( 0, " Basis_Processor::construct_follower_basis_list, the given criteria does not exist" );
                break;
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_volumetric_clustering_of_basis()
    {
        // generate the subphase to enriched basis map
        mXTKModelPtr->mEnrichment->construct_enriched_basis_in_subphase_group_map();

        // loop over the mesh indices
        for ( const auto iMeshIndex : mMeshIndices )
        {
            // creat a cell to keep track of the basis that has been grouped
            moris::Cell< uint > tBasisHasBeenUsed( mXTKModelPtr->mEnrichment->mEnrichmentData( iMeshIndex ).mEnrLvlOfEnrBf.size(), 0 );

            moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( iMeshIndex ).mEnrichedBasisInSubphaseGroup;

            Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( iMeshIndex );

            // Get the number of subphase groups (on the current proc) and resize
            moris_id tNumSubphaseGroups = tBsplineMeshInfo->get_num_SPGs();

            // loop over the SPGs and integrate over the IG cells to compute the volume of the SPG
            for ( size_t iSPGIndex = 0; iSPGIndex < (size_t)tNumSubphaseGroups; iSPGIndex++ )
            {
                // get the subphase group based on the index
                Subphase_Group* tSubphaseGroup = tBsplineMeshInfo->mSubphaseGroups( iSPGIndex );

                // check if the B-spline element has no cuts in it. get the subphase group based on the index
                moris_index tBsplineCellIndex = tSubphaseGroup->get_bspline_cell_index();

                // check the size of the SPGs that are defined on the B-spline Cell
                uint tNumSPGsInCell = tBsplineMeshInfo->mSpgIndicesInBsplineCells( tBsplineCellIndex ).size();

                // B-spline element is not cut thus does not need to be grouped
                if ( 1 == tNumSPGsInCell )
                {
                    continue;
                }

                // otherwise get the IG cells in that SPG and integrate over the volume of those
                const moris::Cell< moris_index >& tIGCellsInGroup = tSubphaseGroup->get_ig_cell_indices_in_group();

                // compute the volume of the cell
                real tVolume = std::accumulate( tIGCellsInGroup.begin(), tIGCellsInGroup.end(), 0.0, [ & ]( real aVol, int index ) { 
                                    const mtk::Cell & tCell = mXTKModelPtr->mCutIntegrationMesh->get_mtk_cell(index) ;
                                    return aVol + tCell.compute_cell_measure(); } );

                // If volume is less than a specific threshold then add it to the list of
                if ( tVolume < 0.1 )
                {
                    moris::Cell< moris_index > const & tEnrichedBFIndices = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                    mBasisData( iMeshIndex ).mBasisInSubspace.push_back( tEnrichedBFIndices );

                    std::for_each( tEnrichedBFIndices.cbegin(), tEnrichedBFIndices.cend(), [ &tBasisHasBeenUsed ]( int aIndex )    //
                            {
                                tBasisHasBeenUsed( aIndex ) += 1;
                            } );
                }
            }

            // Now loop over the dofs that have not been used and group them as points
            for ( moris_index iBasisIndex = 0; iBasisIndex < (moris_index)tBasisHasBeenUsed.size(); iBasisIndex++ )
            {
                // when this b-spline node has not been attached to it.
                if ( 0 == tBasisHasBeenUsed( iBasisIndex ) )
                {
                    mBasisData( iMeshIndex ).mBasisInSubspace.push_back( { iBasisIndex } );
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_basis_using_basis_support( moris_index aMeshIndex )
    {
        // get the enriched to spg map and vice versa
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup     = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        Cell< moris::Matrix< IndexMat > > const &         tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;

        // get the b-spline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the number of active B-spline elements
        uint tNumBspElems = tBsplineMeshInfo->get_num_Bspline_cells();

        // get the basis data and resize the follower basis as the number of enriched basis
        mBasisData( aMeshIndex ).mFollowerBasis.resize( tSubphaseGroupIndsInEnrichedBasis.size(), 1 );

        // allocate root spgs index for the spgs, initialize with -1 to catch the errors later
        mRootSPGIds.resize( tBsplineMeshInfo->get_num_SPGs(), gNoID );
        mRootSPGOwner.resize( tBsplineMeshInfo->get_num_SPGs(), gNoID );

        // loop over the B-spline elements to determine the cut elements
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            moris::Cell< moris_index > const & tSPGIndicesInBsplineCell = tBsplineMeshInfo->get_SPG_indices_in_bspline_cell( iBspElem );

            // B-spline element is not cut thus does not need to be grouped
            if ( 1 == tSPGIndicesInBsplineCell.size() )
            {
                // assign as the root as itself because it is not cut
                mRootSPGIds( tSPGIndicesInBsplineCell( 0 ) ) = tSPGIndicesInBsplineCell( 0 );

                // find the basis that are on this non-cut cell
                moris::Cell< moris_index > const & tEnrichedBasis = tEnrichedBasisInSubphaseGroup( tSPGIndicesInBsplineCell( 0 ) );

                // loop over the basis and integrate the support of each basis function
                for ( const auto& iEnrBasisIndex : tEnrichedBasis )
                {
                    mBasisData( aMeshIndex ).mFollowerBasis( iEnrBasisIndex ) = 0;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_dependent_basis_to_root_map( moris_index aMeshIndex )
    {
        // get number of SPGs
        uint tNumSubphaseGroups = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex )->get_num_SPGs();

        // get the enriched basis to subphase group map
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;

        // resize the dependent basis to root map
        mBasisData( aMeshIndex ).mFollowerBasisOwningCell.resize( mBasisData( aMeshIndex ).mFollowerBasis.size(), gNoIndex );

        // loop over the SPGs
        for ( moris_index iSPGIndex = 0; iSPGIndex < (int)tNumSubphaseGroups; iSPGIndex++ )
        {
            // if it is agglomerated cell
            if ( iSPGIndex != mRootSPGIds( iSPGIndex ) )
            {
                // get the enriched basis indices present in the SPGs
                moris::Cell< moris_index > const & tEnrichedBasis = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                // loop over the enriched basis and if there a is a follower then assign the root as the owner, might be overwritten later
                for ( const auto& iEnrichedBase : tEnrichedBasis )
                {
                    // if it a bad basis then this is true
                    if ( mBasisData( aMeshIndex ).mFollowerBasis( iEnrichedBase ) == 1 )
                    {
                        mBasisData( aMeshIndex ).mFollowerBasisOwningCell( iEnrichedBase ) = mRootSPGIds( iSPGIndex );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::visualize_cell_aggregates( moris_index aMeshIndex )
    {
        // Get the Bspline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // Get the number of subphase groups
        uint tNumSubphaseGroups = tBsplineMeshInfo->get_num_SPGs();

        uint tNumIGCells = mXTKModelPtr->mEnrichedIntegMesh( 0 )->get_num_elems();
        // get the cells that exist in the

        // Fields constructed here
        Cell< std::string > tCellFields = { "Root_Cell_Index_B" + std::to_string( aMeshIndex ) };

        // set field index
        moris_index tFieldIndex = mXTKModelPtr->mEnrichedIntegMesh( 0 )->create_field( tCellFields( 0 ), mtk::EntityRank::ELEMENT, MORIS_INDEX_MAX );

        // create the field data that will store SPG ids
        moris::Matrix< moris::DDRMat > tCellIdField( 1, tNumIGCells, -1.0 );

        // loop over the SPGs and assign dpg id to each cell index
        for ( uint iSPG = 0; iSPG < tNumSubphaseGroups; iSPG++ )
        {
            const moris::Cell< moris_index >& tIGCellIndices = tBsplineMeshInfo->mSubphaseGroups( iSPG )->get_ig_cell_indices_in_group();

            std::for_each( tIGCellIndices.begin(), tIGCellIndices.end(), [ &tCellIdField, this, iSPG ]( moris_index aIGCellIndex )    //
                    { tCellIdField( aIGCellIndex ) = mRootSPGIds( iSPG ); } );
        }

        // add the field data to the mesh
        mXTKModelPtr->mEnrichedIntegMesh( 0 )->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, tCellIdField );

        //---------------------------------------------------------------------------------- */

        // write on IP mesh
        uint tNumIPCells = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_num_elems();

        // set field index
        tFieldIndex                = mXTKModelPtr->mEnrichedInterpMesh( 0 )->create_field( tCellFields( 0 ), mtk::EntityRank::ELEMENT, 0 );
        moris_index tFieldIndexSPG = mXTKModelPtr->mEnrichedInterpMesh( 0 )->create_field( "SPG_index", mtk::EntityRank::ELEMENT, 0 );

        moris::Matrix< moris::DDRMat > tCellIPField( 1, tNumIPCells, -1.0 );
        moris::Matrix< moris::DDRMat > tCellSPGField( 1, tNumIPCells, -1.0 );

        // loop over the SPGs and assign dpg id to each cell index
        for ( uint iSPG = 0; iSPG < tNumSubphaseGroups; iSPG++ )
        {
            // get the all the lagrange cells of the that SPG index
            moris::Cell< moris_index > const & tUIPCIndices = mXTKModelPtr->mEnrichment->get_UIPC_indices_on_SPG( aMeshIndex, iSPG );

            // add the root SPG and SPG index to the field
            std::for_each( tUIPCIndices.begin(), tUIPCIndices.end(), [ &tCellIPField, this, iSPG, &tCellSPGField ]( moris_index aIPCellIndex )    //
                    { tCellIPField( aIPCellIndex ) = mRootSPGIds( iSPG );
                     tCellSPGField( aIPCellIndex ) =  iSPG ; } );
        }

        // add the field data to the mesh
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->add_field_data( tFieldIndex, mtk::EntityRank::ELEMENT, tCellIPField );
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->add_field_data( tFieldIndexSPG, mtk::EntityRank::ELEMENT, tCellSPGField );

        // write this field to the exodus file if there is a cell that is not grouped
        auto it = std::find( mRootSPGIds.begin(), mRootSPGIds.end(), -1 );
        if ( it != mRootSPGIds.end() )
        {
            // output the xtk mesh for debugging
            mXTKModelPtr->mEnrichedIntegMesh( 0 )->write_mesh( mParameterList );

            mXTKModelPtr->mEnrichedInterpMesh( 0 )->write_mesh( mParameterList );

            MORIS_ASSERT( false, "The SPG with the index %d is not grouped, refer to XTK IG mesh for more info", *it );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::perform_basis_extention()
    {
        Tracer tTracer( "XTK", "Basis Processor", "Basis Extension" );

        // generate the sub-phase to enriched basis map
        mXTKModelPtr->mEnrichment->construct_enriched_basis_in_subphase_group_map();

        for ( const auto& iMeshIndex : mMeshIndices )
        {
            this->construct_follower_cells_using_volume( iMeshIndex );

            this->construct_cell_aggregates( iMeshIndex );

            this->construct_follower_basis_using_volume( iMeshIndex );

            this->compute_averaging_weights( iMeshIndex );

            this->visualize_cell_aggregates( iMeshIndex );

            if ( par_size() == 1 )
            {
                this->construct_follower_to_leader_basis_weights_indices( iMeshIndex );
            }
            else
            {
                this->construct_comm_table();
                this->construct_follower_to_leader_basis_weights_indices_mine( iMeshIndex );
                this->construct_follower_to_leader_basis_weights_indices_neighbors( iMeshIndex );
                this->communicate_shared_basis( iMeshIndex );
            }

            this->replace_t_matrices( iMeshIndex );
        }

        this->update_comm_table();
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::perform_cell_agglomeration()
    {
        Tracer tTracer( "XTK", "Basis Processor", "Cell Agglomeration" );

        // generate the sub-phase to enriched basis map
        mXTKModelPtr->mEnrichment->construct_enriched_basis_in_subphase_group_map();

        for ( const auto& iMeshIndex : mMeshIndices )
        {
            this->construct_follower_basis_using_basis_support( iMeshIndex );

            this->construct_cell_aggregates( iMeshIndex );

            this->visualize_cell_aggregates( iMeshIndex );

            this->construct_dependent_basis_to_root_map( iMeshIndex );

            this->override_t_matrices( iMeshIndex );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::compute_averaging_weights( moris_index aMeshIndex )
    {
        // get a reference to pointer to the required data for readability
        Cell< moris::Matrix< IndexMat > > const &         tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup     = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        const Bspline_Mesh_Info*                          tBsplineMeshInfo                  = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< Subphase_Group* > const &            tSPGs                             = tBsplineMeshInfo->mSubphaseGroups;
        uint                                              tNumEnrichedBF                    = tSubphaseGroupIndsInEnrichedBasis.size();

        // assign the correct size to the averaging weights
        mBasisData( aMeshIndex ).mAveragingWeights.resize( tNumEnrichedBF );

        // loop over the enriched basis functions and determine how many spgs they are active
        for ( uint iEnrBasisIndex = 0; iEnrBasisIndex < tNumEnrichedBF; iEnrBasisIndex++ )
        {
            // if it is a bad basis then process further, good basis do not get effected
            if ( 1 == mBasisData( aMeshIndex ).mFollowerBasis( iEnrBasisIndex ) )
            {
                // get the SPG indices for the enriched BF index
                moris::Matrix< IndexMat > const & tSPGIndices = tSubphaseGroupIndsInEnrichedBasis( iEnrBasisIndex );

                //  reserve enough space for the averaging weights for a specific enriched basis function
                mBasisData( aMeshIndex ).mAveragingWeights( iEnrBasisIndex ).reserve( tSPGIndices.numel() );

                // loop over the SPGs and get the volume of each
                for ( const auto& iSPGIndex : tSPGIndices )
                {
                    // get the subphase group based on the index
                    Subphase_Group* tSubphaseGroup = tSPGs( iSPGIndex );

                    // otherwise get the IG cells in that SPG and integrate over the volume of those
                    const moris::Cell< moris_index >& tIGCellsInGroup = tSubphaseGroup->get_ig_cell_indices_in_group();

                    // compute the volume of the cell
                    real tVolume = std::accumulate( tIGCellsInGroup.begin(), tIGCellsInGroup.end(), 0.0, [ & ]( real aVol, int index ) { 
                                    const mtk::Cell & tCell = mXTKModelPtr->mCutIntegrationMesh->get_mtk_cell(index) ;
                                    return aVol + tCell.compute_cell_measure(); } );

                    // add the volume as a first step
                    mBasisData( aMeshIndex ).mAveragingWeights( iEnrBasisIndex ).push_back( tVolume );
                }

                // create a reference for readability and ease of use,  this is a cell of volumes
                moris::Cell< real >& tAveragingWeightsInEnriched = mBasisData( aMeshIndex ).mAveragingWeights( iEnrBasisIndex );

                // normalize the averaging wights against their sum
                real tAccumulatedVolume = std::accumulate( tAveragingWeightsInEnriched.begin(), tAveragingWeightsInEnriched.end(), 0.0 );

                // normalize the values
                std::transform( tAveragingWeightsInEnriched.begin(), tAveragingWeightsInEnriched.end(), tAveragingWeightsInEnriched.begin(), [ &tAccumulatedVolume ]( real aVolume )    //
                        { return aVolume / tAccumulatedVolume; } );
            }
        }

        /*-----------------------------------------------------------------------------------------*/
        /* construct the transpose of the map above mAveragingWeightsSPGBased*/

        /// initialize a counter to keep track of SPGs that need an extension
        uint aSPGIndex = 0;

        // define a lambda function
        // TODO: add a condition to certify only the owned SPGs to optimize for the size
        auto tSPGNeedsEXtension = [ &aSPGIndex, &tBsplineMeshInfo ]( moris_index aRootSPGId ) {
            // this function checks for a given index the id of the SPG equals to the its root ID provided by mRootSPGIds
            // it also increments the SPG index for the nex iteration
            return aRootSPGId != tBsplineMeshInfo->get_id_for_spg_index( aSPGIndex++ );
        };

        // loop over the SPGs and determine if their root is the same as SPG index
        uint tExtenedSPGSize = std::count_if( mRootSPGIds.begin(), mRootSPGIds.end(), tSPGNeedsEXtension );

        // reserve enough space in the map for the SPGs that need extension
        mBasisData( aMeshIndex ).mAveragingWeightsSPGBased.reserve( tExtenedSPGSize );

        // loop over the owned SPGs and get the averaging weights for the enriched basis functions
        for ( uint iEnrBasisIndex = 0; iEnrBasisIndex < tNumEnrichedBF; iEnrBasisIndex++ )
        {
            // it it is a good basis then skip it
            if ( 0 == mBasisData( aMeshIndex ).mFollowerBasis( iEnrBasisIndex ) )
            {
                continue;
            }

            // if it is a bad basis then process further
            // get the SPG indices for the enriched BF index
            moris::Matrix< IndexMat > const & tSPGIndices                 = tSubphaseGroupIndsInEnrichedBasis( iEnrBasisIndex );
            moris::Cell< real > const &       tAveragingWeightsInEnriched = mBasisData( aMeshIndex ).mAveragingWeights( iEnrBasisIndex );

            MORIS_ASSERT( tSPGIndices.numel() == tAveragingWeightsInEnriched.size(), "The number of SPGs and the number of averaging weights are not the same" );

            // loop over the SPGs and get the volume of each
            for ( uint iSPGOrd = 0; iSPGOrd < tSPGIndices.numel(); iSPGOrd++ )
            {
                // get the SPG index
                moris_index iSPGIndex = tSPGIndices( iSPGOrd );

                // get all the enriched basis that exist on the SPG
                moris::Cell< moris_index > const & tEnrichedBasisInSPG = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                // resize the weights for the enriched basis in the SPG
                mBasisData( aMeshIndex ).mAveragingWeightsSPGBased[ iSPGIndex ].resize( tEnrichedBasisInSPG.size(), gNoIndex );

                // find the local index of the enriched basis in the SPG
                auto tEnrichedBasisOrdinalIterator = std::find( tEnrichedBasisInSPG.begin(), tEnrichedBasisInSPG.end(), iEnrBasisIndex );

                // get the local index of the enriched basis in the SPG
                moris_index tEnrichedBasisOrdinal = std::distance( tEnrichedBasisInSPG.begin(), tEnrichedBasisOrdinalIterator );

                // Otherwise get the IP cells in that SPG and integrate over the volume of those
                mBasisData( aMeshIndex ).mAveragingWeightsSPGBased[ iSPGIndex ]( tEnrichedBasisOrdinal ) = tAveragingWeightsInEnriched( iSPGOrd );
            }
        }
    }

    //----------------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_cells_using_volume( moris_index aMeshIndex )
    {
        // compute the threshold volume to decide good/bad basis
        real tVolThreshold = this->compute_threshold_volume( aMeshIndex );

        // get the spg to enr basis and its transpose( enr basis to spg map)
        // moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup     = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        Cell< moris::Matrix< IndexMat > > const & tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;

        // get the bspline mesh info to access the spgs
        const Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the spgs from b-spline mesh info
        moris::Cell< Subphase_Group* > const & tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the number of active B-spline elements ( this is un-enriched)
        uint tNumBspElems      = tBsplineMeshInfo->get_num_Bspline_cells();
        uint tNumEnrichedBasis = tSubphaseGroupIndsInEnrichedBasis.size();

        // Get the number of subphase groups (on the current proc) and resize
        uint tNumSubphaseGroups = tBsplineMeshInfo->get_num_SPGs();

        // allocate size for the basis info
        // for follower basis assume all are follower and then change the leader one to zero
        mBasisData( aMeshIndex ).mFollowerBasis.resize( tNumEnrichedBasis, 1 );
        mBasisData( aMeshIndex ).mFollowerToLeaderBasis.resize( tNumEnrichedBasis );
        mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights.resize( tNumEnrichedBasis );
        mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners.resize( tNumEnrichedBasis, 0 );

        // allocate root spgs index for the spgs, initialize with -1 to catch the errors later
        mRootSPGIds.resize( tNumSubphaseGroups, gNoID );
        mRootSPGOwner.resize( tNumSubphaseGroups, gNoIndex );

        // loop over the b-spline elements to get the SPGs living on them
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            // get the SPGs that live on the bspline element
            moris::Cell< moris_index > const & tSPGIndicesInBsplineCell = tBsplineMeshInfo->get_SPG_indices_in_bspline_cell( iBspElem );

            // if only one spg then B-spline element is not cut
            if ( 1 == tSPGIndicesInBsplineCell.size() )
            {
                // get the subphase group based on the index
                const Subphase_Group* tSubphaseGroup = tSPGs( tSPGIndicesInBsplineCell( 0 ) );

                if ( tSubphaseGroup->get_owner() == par_rank() )
                {
                    // assign as the root as itself because it is not cut
                    mRootSPGIds( tSPGIndicesInBsplineCell( 0 ) )   = tSubphaseGroup->get_id();
                    mRootSPGOwner( tSPGIndicesInBsplineCell( 0 ) ) = par_rank();
                    // go to the next b-spline element
                    continue;
                }
                continue;
            }

            // if the bspline element is cut find if the volume is below a certain threshold;
            for ( const auto& iSPGIndex : tSPGIndicesInBsplineCell )
            {
                // get the subphase group based on the index
                const Subphase_Group* tSubphaseGroup = tSPGs( iSPGIndex );

                // if SPG is not owned then skip
                if ( tSubphaseGroup->get_owner() != par_rank() )
                {
                    continue;
                }

                // otherwise get the IG cells in that SPG and integrate over the volume of those
                const moris::Cell< moris_index >& tIGCellsInGroup = tSubphaseGroup->get_ig_cell_indices_in_group();

                // compute the volume of the cell
                real tVolume = std::accumulate( tIGCellsInGroup.begin(), tIGCellsInGroup.end(), 0.0, [ & ]( real aVol, int index ) { 
                                    const mtk::Cell & tCell = mXTKModelPtr->mCutIntegrationMesh->get_mtk_cell(index) ;
                                    return aVol + tCell.compute_cell_measure(); } );

                // If volume is less than a specific threshold then add it to the list of
                if ( tVolume >= tVolThreshold )
                {
                    // if a volume threshold is met mark the spgs as its own root
                    mRootSPGIds( iSPGIndex )   = tSubphaseGroup->get_id();
                    mRootSPGOwner( iSPGIndex ) = par_rank();
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    real
    Basis_Processor::compute_threshold_volume( moris_index aMeshIndex )
    {
        // read off the volume fraction from the parameter list
        real tVolumeFraction = mParameterList->get< real >( "volume_fraction" );

        // get the all the lagrange cells of the that SPG index 0
        // will be used to get number of lagrange cells
        moris::Cell< moris_index > const & tUIPCIndices = mXTKModelPtr->mEnrichment->get_UIPC_indices_on_SPG( aMeshIndex, 0 );

        // get the enriched IP cell with index 0
        mtk::Cell& tFirstCell = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_mtk_cell( 0 );

        // compute the volume of the first cell
        moris::real tVolumeLagrangeCell = tFirstCell.compute_cell_measure();

        // compute the volume of the threshold bspline element
        return tVolumeFraction * tVolumeLagrangeCell * (real)tUIPCIndices.size();
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_cell_association( moris_index aMeshIndex )
    {
        // get maximum extension of the basis length
        moris_index mMaxExtention = 2;

        // get the bspline mesh info to access the spgs
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the spgs from b-spline mesh info
        moris::Cell< Subphase_Group* > const & tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the bspline cell
        moris::Cell< mtk::Cell* > const & tBsplineCells = tBsplineMeshInfo->mBsplineCells;

        // get the subphase connectivity
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseGroupNeighborhood = mXTKModelPtr->mCutIntegrationMesh->get_subphase_group_neighborhood( aMeshIndex );

        // initialize cells to keep possible candidate root cells and their distance
        moris::Cell< moris_index > tCandidateRootCell;
        moris::Cell< moris_index > tDistanceCell;

        // loop over the spgs and process the ones that do not have root spgs
        for ( size_t iSPGIndex = 0; iSPGIndex < mRootSPGIds.size(); iSPGIndex++ )
        {
            // this mean the SPG does not have not have a root cell
            tCandidateRootCell.resize( 0 );
            tDistanceCell.resize( 0 );

            // this means that the root spg is not assigned
            if ( -1 == mRootSPGIds( iSPGIndex ) )
            {
                // get the first degree neighbor of the SPG
                // TODO: only first degree neighbor is considered, need to consider higher degree neighbors
                moris::Cell< moris_index > tNeighbour;

                for ( moris_index iNeighborDegree = 1; iNeighborDegree <= mMaxExtention; iNeighborDegree++ )
                {
                    // reserve the space for the neighbors
                    tNeighbour.reserve( mSpatialDim * iNeighborDegree );
                    tSubphaseGroupNeighborhood->get_kth_degree_neighbors( iSPGIndex, iNeighborDegree, tNeighbour );

                    // between the neighbor get the ones that are is a root cell
                    std::copy_if( tNeighbour.begin(), tNeighbour.end(), std::back_inserter( tCandidateRootCell ), [ & ]( moris_index aNeighbourSPGIndex )    //
                            {
                                return mRootSPGIds( aNeighbourSPGIndex ) == aNeighbourSPGIndex;
                            } );

                    // assign as root cell the cell with smaller index
                    if ( tCandidateRootCell.size() )
                    {
                        // define a lambda function with sorting criteria
                        auto tSortBasedOnTheBsplineCellIndex = [ & ]( moris_index aSPGIndex1, moris_index aSPGIndex2 ) {
                            return tSPGs( aSPGIndex1 )->get_bspline_cell_index() < tSPGs( aSPGIndex2 )->get_bspline_cell_index();
                        };

                        // this sorting is necessary such that in the event of a tie in euclidean distance the cell with smallest B-spline index
                        std::sort( tCandidateRootCell.begin(), tCandidateRootCell.end(), tSortBasedOnTheBsplineCellIndex );

                        // loop over the candidate cells that have the find the closest
                        // TODO: could be improved by an elimination strategy instead of computing all distances
                        for ( const auto& iNeighbourSPGIndex : tCandidateRootCell )
                        {
                            // compute euclidean distance of the two cells
                            uint         tLevel = 0;
                            const luint* tIJK   = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tBsplineCells( tSPGs( iSPGIndex )->get_bspline_cell_index() ), tLevel );

                            // find the root spg of the neighbor SPG and its ijk
                            moris_index  tRootOfNeighbourSPG = mRootSPGIds( iNeighbourSPGIndex );
                            uint         tLevelRoot          = 0;
                            const luint* tIJKRoot            = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tBsplineCells( tSPGs( tRootOfNeighbourSPG )->get_bspline_cell_index() ), tLevelRoot );

                            // initialize the difference
                            moris_index tSumOfSquaredDifferences = 0;

                            // find the distance
                            for ( uint iDim = 0; iDim < mSpatialDim; iDim++ )
                            {
                                // define a difference to get the correct subtraction operator for the unsigned ints
                                moris_index tDiff = ( tIJKRoot[ iDim ] < tIJK[ iDim ] ) ? -( tIJK[ iDim ] - tIJKRoot[ iDim ] ) : ( tIJKRoot[ iDim ] - tIJK[ iDim ] );

                                // sum up the value
                                tSumOfSquaredDifferences += tDiff * tDiff;
                            }

                            tDistanceCell.push_back( tSumOfSquaredDifferences );
                        }

                        // find the one with smallest distance
                        auto tSelectedNeighborSPGIndex = std::min_element( tDistanceCell.begin(), tDistanceCell.end() );
                        auto tIndex                    = std::distance( tDistanceCell.begin(), tSelectedNeighborSPGIndex );
                        mRootSPGIds( iSPGIndex )       = mRootSPGIds( tCandidateRootCell( tIndex ) );

                        // a root has been found, so break the loop over neighbors and go the next SPG
                        break;
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_to_leader_basis_weights_indices( moris_index aMeshIndex )
    {
        // Get a reference or a pointer to the required data
        Bspline_Mesh_Info*                        tBsplineMeshInfo                  = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        Cell< moris::Matrix< IndexMat > > const & tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
        moris::Cell< Subphase_Group* > const &    tSPGs                             = tBsplineMeshInfo->mSubphaseGroups;
        moris::Cell< mtk::Cell* > const &         tBsplineCells                     = tBsplineMeshInfo->mBsplineCells;

        // loop over the basis
        for ( uint iBasisIndex = 0; iBasisIndex < mBasisData( aMeshIndex ).mFollowerBasis.size(); iBasisIndex++ )
        {
            // check if it is a bad basis
            if ( mBasisData( aMeshIndex ).mFollowerBasis( iBasisIndex ) == 1 )
            {
                // get the SPGs that this basis is active
                const moris::Matrix< IndexMat >& tSPGSinBasis = tSubphaseGroupIndsInEnrichedBasis( iBasisIndex );

                // loop over the SPGs that are in the support of the basis
                for ( uint iSPGOrd = 0; iSPGOrd < tSPGSinBasis.numel(); iSPGOrd++ )
                {
                    // spg index
                    moris_index iSPGIndex = tSPGSinBasis( iSPGOrd );

                    // get the averaging weight
                    real tAveragingWeight = mBasisData( aMeshIndex ).mAveragingWeights( iBasisIndex )( iSPGOrd );

                    // check if the SPG have a different root, then it needs extension
                    moris_index tRootSPGIndex = tBsplineMeshInfo->get_index_for_spg_id( mRootSPGIds( iSPGIndex ) );

                    // get the bspline cell of the root cell
                    moris_index tRootBSplineCellIndex = tSPGs( tRootSPGIndex )->get_bspline_cell_index();

                    // get the bspline cell of the root cell
                    moris_index tExtentionBSplineCellIndex = tSPGs( iSPGIndex )->get_bspline_cell_index();

                    // initialize the basis extension/projection data
                    moris::Cell< moris::Cell< const mtk::Vertex* > > tRootBsplineBasis;
                    moris::Cell< const mtk::Vertex* >                tExtendedBsplineBasis;
                    moris::Cell< Matrix< DDRMat > >                  tWeights;

                    // get the L2-projection matrix along from the root to the extended b-spline
                    mXTKModelPtr->mBackgroundMesh->get_L2_projection_matrix( aMeshIndex, tBsplineCells( tRootBSplineCellIndex ), tBsplineCells( tExtentionBSplineCellIndex ), tRootBsplineBasis, tExtendedBsplineBasis, tWeights );

                    // This two maps help to get the the enriched basis index from the
                    moris::Cell< moris_index > const & tBGBasisIndicesRoot = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisIndices( tRootSPGIndex );
                    moris::Cell< moris_index > const & tBGBasisLevelsRoot  = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisEnrLev( tRootSPGIndex );
                    IndexMap                           tBasisToLocalIndexMapRoot;

                    // need to convert the first one to map
                    convert_cell_to_map( tBGBasisIndicesRoot, tBasisToLocalIndexMapRoot );

                    // get the BG basis function index of the basis
                    moris::Cell< moris_index > const & tNonEnrBfIndForEnrBfInd = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mNonEnrBfIndForEnrBfInd;
                    moris_index                        tNonEnrBasisIndex       = tNonEnrBfIndForEnrBfInd( iBasisIndex );

                    // find the location in the L2 projection matrix
                    auto tIterator = std::find_if( tExtendedBsplineBasis.begin(), tExtendedBsplineBasis.end(),    //
                            [ &tNonEnrBasisIndex ]( const mtk::Vertex* aVertex ) { return aVertex->get_index() == tNonEnrBasisIndex; } );

                    // MORIS_ASSERT(tIterator ==  tExtendedBsplineBasis.end(), " The index does not exist in the list of current given basis");

                    // find the ordinal of the non-enriched basis index in the extended basis
                    uint tNonEnrBasisOrd = std::distance( tExtendedBsplineBasis.begin(), tIterator );

                    // get the map relating non-enriched bf index and level to enriched bf index
                    moris::Cell< Matrix< IndexMat > > const & tBasisEnrichmentIndices = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mBasisEnrichmentIndices;

                    // find the max size in enriched basis index
                    auto tCellMaxsize = std::max_element( tRootBsplineBasis.begin(), tRootBsplineBasis.end(),                      //
                            []( moris::Cell< const mtk::Vertex* >& aCellFirst, moris::Cell< const mtk::Vertex* >& aCellSecond )    //
                            { return aCellFirst.size() < aCellSecond.size(); } );

                    // reserve enough space for the enriched basis indices
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iBasisIndex ).reserve( tCellMaxsize->size() );
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iBasisIndex ).reserve( tCellMaxsize->size() );
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( iBasisIndex ).reserve( tCellMaxsize->size() );

                    // loop over the unenriched BSpline vertices
                    for ( uint iRootBGBasisOrd = 0; iRootBGBasisOrd < tRootBsplineBasis( tNonEnrBasisOrd ).size(); iRootBGBasisOrd++ )
                    {
                        // get the root basis index based on the ordinal
                        moris_index tBGBasisIndexRoot = tRootBsplineBasis( tNonEnrBasisOrd )( iRootBGBasisOrd )->get_index();

                        // find the location of the root basis that is being enriched
                        auto        tIter            = tBasisToLocalIndexMapRoot.find( tBGBasisIndexRoot );
                        moris_index tLocalBasisIndex = tIter->second;

                        // obtain the level of non-enriched root basis
                        moris_index tEnrLevel = tBGBasisLevelsRoot( tLocalBasisIndex );

                        // get the enriched basis index of the root basis
                        moris_index tEnrichedBasisIndex = tBasisEnrichmentIndices( tBGBasisIndexRoot )( tEnrLevel );

                        // determine the old basis that will be replaced by which enriched root basis and what weights
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iBasisIndex ).push_back( tEnrichedBasisIndex );
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iBasisIndex ).push_back( tWeights( tNonEnrBasisOrd )( iRootBGBasisOrd ) * tAveragingWeight );
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( iBasisIndex ).push_back( tSPGs( iSPGIndex )->get_owner() );
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::replace_t_matrices( moris_index aMeshIndex )
    {
        // get a reference or a pointer to required data
        Cell< Interpolation_Cell_Unzipped* > const & tEnrichedIPCells                  = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_interpolation_cells();
        Matrix< IndexMat > const &                   tBasisEnrichmentIndices           = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisIndexToId;
        Cell< moris::Matrix< IndexMat > > const &    tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
        Enriched_Interpolation_Mesh&                 tEnrInterpMesh                    = mXTKModelPtr->get_enriched_interp_mesh( aMeshIndex );

        // get the bspline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        // get the subphase groups
        moris::Cell< Subphase_Group* > const & tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // define an unordered set to keep track of the vertices that have been already processed
        std::unordered_set< moris_index > tListOfBGVertices;
        tListOfBGVertices.reserve( tSubphaseGroupIndsInEnrichedBasis.size() );

        // find the max size that a basis will be replaced by
        auto tCellMaxsize = std::max_element( mBasisData( aMeshIndex ).mFollowerToLeaderBasis.begin(), mBasisData( aMeshIndex ).mFollowerToLeaderBasis.end(),    //
                []( moris::Cell< moris_index >& aCellFirst, moris::Cell< moris_index >& aCellSecond )                                                            //
                { return aCellFirst.size() < aCellSecond.size(); } );

        // loop over the subphase groups to find the problematic ones
        for ( size_t iSPGIndex = 0; iSPGIndex < mRootSPGIds.size(); iSPGIndex++ )
        {
            if ( tSPGs( iSPGIndex )->get_owner() != par_rank() ) continue;

            // if it is an extended cell
            if ( tSPGs( iSPGIndex )->get_id() != mRootSPGIds( iSPGIndex ) )
            {
                // get the all the lagrange cells of the that SPG index
                moris::Cell< moris_index > const & tUIPCIndices = mXTKModelPtr->mEnrichment->get_UIPC_indices_on_SPG( aMeshIndex, iSPGIndex );

                // loop over the unzipped lagrange cells and evaluate t-matrices
                for ( const auto& iUIPCIndex : tUIPCIndices )
                {
                    // get base cell of the lagrange mesh
                    const Interpolation_Cell_Unzipped* tConstUIPC = tEnrichedIPCells( iUIPCIndex );

                    // Get the lagrange vertices of the UIPC to obtain the vertex interpolation object of each lagrange node
                    moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tXTKUIPVs = tConstUIPC->get_xtk_interpolation_vertices();

                    // loop over the unzipped interpolation vertices and get their t-matrix info to overwrite
                    for ( size_t iLocalVertIndex = 0; iLocalVertIndex < tXTKUIPVs.size(); iLocalVertIndex++ )
                    {
                        // if the vertex has been already processed, skip it
                        if ( tListOfBGVertices.count( tXTKUIPVs( iLocalVertIndex )->get_index() ) > 0 )
                        {
                            continue;
                        }

                        // get the vertex enrichment object
                        xtk::Vertex_Enrichment* tVertexEnrichment = tXTKUIPVs( iLocalVertIndex )->get_xtk_interpolation( aMeshIndex );

                        // add the vertex to the list of processed vertices
                        tListOfBGVertices.insert( tXTKUIPVs( iLocalVertIndex )->get_index() );

                        // access the basis indices and weights
                        moris::Matrix< moris::IndexMat > const & tBasisIndices = tVertexEnrichment->get_basis_function_indices();
                        moris::Matrix< moris::IndexMat > const & tBasisIds     = tVertexEnrichment->get_basis_function_ids();
                        moris::Matrix< moris::DDRMat >&          tBasisWeights = tVertexEnrichment->get_basis_function_weights();
                        moris::Matrix< IdMat >                   tBasisOwners  = tVertexEnrichment->get_owners();
                        IndexMap&                                tBasisMap     = tVertexEnrichment->get_basis_function_map();
                        tBasisMap.clear();

                        // probably need to do a resize call here
                        // overallocation basis Indices , assume that every basis index need to be replaced
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisIndices( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisIds( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        moris::Matrix< DDRMat >          tAgglomeratedBasisWeights( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisOwners( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        // tBasisMap.reserve( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ) ); // FIXME: only works if the Mini_Map is implemented using a std::unordered_map

                        // initialize a counter to count how many basis will be added and replaced
                        uint tBasisCounter = 0;

                        // loop over the basis indices and replace the ones that are not in the follower basis
                        for ( uint iBC = 0; iBC < tBasisIndices.numel(); iBC++ )
                        {
                            // get a reference to the basis index and weight
                            moris::moris_index tBasisIndex  = tBasisIndices( iBC );
                            real&              tBasisWeight = tBasisWeights( iBC );

                            // if it is a good basis keep it
                            if ( 0 == mBasisData( aMeshIndex ).mFollowerBasis( tBasisIndex ) )
                            {
                                // if it iis not in the map add it to map
                                if ( tBasisMap.find( tBasisIndex ) == tBasisMap.end() )
                                {
                                    tAgglomeratedBasisIndices( tBasisCounter ) = tBasisIndex;
                                    tAgglomeratedBasisIds( tBasisCounter )     = tBasisIds( iBC );
                                    tAgglomeratedBasisWeights( tBasisCounter ) = tBasisWeights( iBC );
                                    tAgglomeratedBasisOwners( tBasisCounter )  = tBasisOwners( iBC );
                                    tBasisMap[ tBasisIndex ]                   = tBasisCounter;

                                    tBasisCounter++;
                                }
                                else
                                {
                                    // if it is repetitive just add weights
                                    moris_index tLocalBasisIndex = tBasisMap[ tBasisIndex ];
                                    tAgglomeratedBasisWeights( tLocalBasisIndex ) += tBasisWeight;
                                }
                            }
                            else
                            {
                                // obtain the replacement basis
                                moris::Cell< moris_index > const & tFollowerBasis        = mBasisData( aMeshIndex ).mFollowerToLeaderBasis( tBasisIndex );
                                moris::Cell< real > const &        tFollowerBasisWeights = mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( tBasisIndex );
                                moris::Cell< real > const &        tFollowerBasisOwners  = mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( tBasisIndex );

                                // loop over and replace the t matrix
                                for ( uint iSlaveBasisOrd = 0; iSlaveBasisOrd < tFollowerBasis.size(); iSlaveBasisOrd++ )
                                {
                                    moris_index const & tFollowerBasisIndex = tFollowerBasis( iSlaveBasisOrd );

                                    //  if not found in the map
                                    if ( tBasisMap.find( tFollowerBasisIndex ) == tBasisMap.end() )
                                    {
                                        // assign the basis index and weight
                                        tAgglomeratedBasisIndices( tBasisCounter ) = tFollowerBasisIndex;
                                        tAgglomeratedBasisIds( tBasisCounter )     = tEnrInterpMesh.get_enr_basis_id_from_enr_basis_index( aMeshIndex, tFollowerBasisIndex );
                                        tAgglomeratedBasisWeights( tBasisCounter ) = tBasisWeight * tFollowerBasisWeights( iSlaveBasisOrd );
                                        tAgglomeratedBasisOwners( tBasisCounter )  = tFollowerBasisOwners( iSlaveBasisOrd );
                                        tBasisMap[ tFollowerBasisIndex ]           = tBasisCounter;

                                        // increment the counter
                                        tBasisCounter++;
                                    }
                                    // if found in the map
                                    else
                                    {
                                        // find the local index and add the weight
                                        moris_index tLocalBasisIndex = tBasisMap[ tFollowerBasisIndex ];
                                        tAgglomeratedBasisWeights( tLocalBasisIndex ) += tBasisWeight * tFollowerBasisWeights( iSlaveBasisOrd );
                                    }
                                }

                                // This is done for the purpose of writing the basis functions
                                // keeping old indices such that it can be assured they are 0 and indices are not mixed
                                if ( mParameterList->get< bool >( "write_basis_functions" ) )
                                {
                                    // assign a basis of weight zero if replaces with the bad basis
                                    tAgglomeratedBasisIndices( tBasisCounter ) = tBasisIndex;
                                    tAgglomeratedBasisIds( tBasisCounter )     = tBasisEnrichmentIndices( tBasisIndex );
                                    tAgglomeratedBasisWeights( tBasisCounter ) = 0.0;
                                    tBasisMap[ tBasisIndex ]                   = tBasisCounter;

                                    // increment the counter
                                    tBasisCounter++;
                                }
                            }
                        }

                        // resize the basis vectors with the correct size
                        tAgglomeratedBasisIndices.resize( tBasisCounter, 1 );
                        tAgglomeratedBasisWeights.resize( tBasisCounter, 1 );
                        tAgglomeratedBasisIds.resize( tBasisCounter, 1 );
                        tAgglomeratedBasisOwners.resize( tBasisCounter, 1 );

                        //  add basis information to the vertex enrichment
                        tVertexEnrichment->add_basis_information( tAgglomeratedBasisIndices, tAgglomeratedBasisIds );
                        tVertexEnrichment->add_basis_function_weights( tAgglomeratedBasisIndices, tAgglomeratedBasisWeights );
                        tVertexEnrichment->add_basis_function_owners( tAgglomeratedBasisIndices, tAgglomeratedBasisOwners );
                    }
                }
            }
        }
    }

    // ==================================================================================

    void
    Basis_Processor::construct_cell_aggregates( moris_index aMeshIndex )
    {
        // get the graph connectivity for the subphase groups
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseGroupNeighborhood = mXTKModelPtr->mCutIntegrationMesh->get_subphase_group_neighborhood( aMeshIndex );

        // get the bspline mesh info
        Bspline_Mesh_Info*                tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< mtk::Cell* > const & tBSpCells        = tBsplineMeshInfo->mBsplineCells;

        // get the subphase groups
        moris::Cell< Subphase_Group* >& tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // initialize candidate root cell and distance cell( candidate root cell is the list of potential root cells for the current SPG)
        // Distance cell is the list of distances from the current SPG to the candidate root cells
        moris::Cell< moris_index > tCandidateRootCell;
        moris::Cell< moris_index > tDistanceCell;

        // reserve space for the candidate root cell and distance cell
        tCandidateRootCell.reserve( mSpatialDim * 2 );
        tDistanceCell.reserve( mSpatialDim * 2 );

        // initialize loop counter that kills the run if the search for a root B-spline element takes too long
        uint tLoopCounter = 0;

        // until all SPGs are aggregated run the graph algorithm to generate aggregates
        while ( this->determine_stopping_criteria() )
        {
            // make sure the code does not get stuck in an infinite while-loop here
            // FIXME: a smarter stopping criteria would be good, this is only a temporary fix to the problem
            MORIS_ERROR(
                    tLoopCounter++ < 100,
                    "Basis_Processor::construct_cell_aggregates() - "
                    "Could not find root B-spline element to extend basis functions from within 100 search iterations." );

            // exchange information about nearest neighbors across processor boundaries (in parallel)
            this->perform_nearest_neighbor_exchange( aMeshIndex );

            // loop over the cut bspline cells ( corrsponding spgs) that are owned
            for ( const auto& iSPGIndex : tBsplineMeshInfo->mOwnedSubphaseGroupIndices )
            {
                // clear the data in the cell while preserving the capacity
                tCandidateRootCell.clear();
                tDistanceCell.clear();

                // this means that it is not assigned
                if ( mRootSPGIds( iSPGIndex ) == -1 )
                {
                    std::shared_ptr< moris::Cell< moris_index > > tNeighbour = tSubphaseGroupNeighborhood->mSubphaseToSubPhase( iSPGIndex );

                    // get the cells that are already touched
                    std::copy_if( tNeighbour->begin(), tNeighbour->end(), std::back_inserter( tCandidateRootCell ), [ & ]( moris_index aNeighborSPGIndex )    //
                            {
                                return mRootSPGIds( aNeighborSPGIndex ) != -1;
                            } );

                    // assign as root cell the cell with smaller index
                    if ( tCandidateRootCell.size() )
                    {
                        auto tSortBasedOnTheBsplineCellIndex = [ & ]( moris_index aSPGIndex1, moris_index aSPGIndex2 ) {
                            return tSPGs( aSPGIndex1 )->get_bspline_cell_index() < tSPGs( aSPGIndex2 )->get_bspline_cell_index();
                        };

                        // this sorting is necessary such that in the event of a tie in euclidean distance the cell with smallest B-spline index
                        std::sort( tCandidateRootCell.begin(), tCandidateRootCell.end(), tSortBasedOnTheBsplineCellIndex );

                        // find the euclidean distance between all candidate points
                        // TODO: could be improved by an elimination strategy instead of computing all distances
                        for ( const auto& iNeighbourSPGIndex : tCandidateRootCell )
                        {
                            // compute euclidean distance of the two cells
                            const luint* tIJK = mHMRHelper( aMeshIndex )->get_ijk_bspline_cell( tBSpCells( tSPGs( iSPGIndex )->get_bspline_cell_index() ) );

                            // find the root spg of the neighbor SPG and its ijk
                            bool tRootSPGExistsOnPartition = tBsplineMeshInfo->spg_exists_on_partition( mRootSPGIds( iNeighbourSPGIndex ) );

                            // initialize the root of the neighbor SPG
                            moris_index tRootOfNeighbourSPG = iNeighbourSPGIndex;

                            // define the root of the neighbor SPG, if it is on processor take the root otherwise take the neighbor itself
                            if ( tRootSPGExistsOnPartition )
                            {
                                tRootOfNeighbourSPG = tBsplineMeshInfo->get_index_for_spg_id( mRootSPGIds( iNeighbourSPGIndex ) );
                            }

                            // get the ijk of the root of the neighbor SPG
                            const luint* tIJKRoot = mHMRHelper( aMeshIndex )->get_ijk_bspline_cell( tBSpCells( tSPGs( tRootOfNeighbourSPG )->get_bspline_cell_index() ) );

                            // initialize the difference
                            moris_index tSumOfSquaredDifferences = 0;

                            // find the distance
                            for ( uint iDim = 0; iDim < mSpatialDim; iDim++ )
                            {
                                // define a difference to get the correct subtraction operator for the unsigned ints
                                moris_index tDiff = ( tIJKRoot[ iDim ] < tIJK[ iDim ] ) ? -( tIJK[ iDim ] - tIJKRoot[ iDim ] ) : ( tIJKRoot[ iDim ] - tIJK[ iDim ] );

                                // sum up the value
                                tSumOfSquaredDifferences += tDiff * tDiff;
                            }

                            // store the distance
                            tDistanceCell.push_back( tSumOfSquaredDifferences );
                        }

                        // find the corresponding root cell based on the distance
                        auto tSelectedNeighborSPGIndex = std::min_element( tDistanceCell.begin(), tDistanceCell.end() );
                        auto tIndex                    = std::distance( tDistanceCell.begin(), tSelectedNeighborSPGIndex );
                        mRootSPGIds( iSPGIndex )       = mRootSPGIds( tCandidateRootCell( tIndex ) );
                        mRootSPGOwner( iSPGIndex )     = mRootSPGOwner( tCandidateRootCell( tIndex ) );

                    }    // end if: there's a candidate root cell

                }        // end if: no root SPG has been determined yet

            }            // end for: each sub-phase group (i.e. unzipped B-spline elements)

        }                // end while: stopping criteria is not fulfilled

    }                    // end function: Basis_Processor::construct_cell_aggregates()

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::override_t_matrices( moris_index aMeshIndex )
    {
        // get the bspline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the subphase groups
        moris::Cell< Subphase_Group* > tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get a reference to the enriched interpolation cells
        Cell< Interpolation_Cell_Unzipped* > const & tEnrichedIPCells = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_interpolation_cells();

        // get the number of subphase groups
        uint tNumSubphaseGroups = mXTKModelPtr->mCutIntegrationMesh->get_num_subphase_groups( aMeshIndex );

        // loop over the SPGs and if an SPG has a different root cell than itself process further
        for ( moris_index iSPGIndex = 0; iSPGIndex < (int)tNumSubphaseGroups; iSPGIndex++ )
        {
            // if the root cell of the SPG is not itself
            if ( mRootSPGIds( iSPGIndex ) != iSPGIndex )
            {
                // get the bspline cell of the root cell
                moris_index tRootBSplineCellIndex = tSPGs( mRootSPGIds( iSPGIndex ) )->get_bspline_cell_index();

                // get the all the lagrange cells of the that SPG index
                moris::Cell< moris_index > const & tUIPCIndices = mXTKModelPtr->mEnrichment->get_UIPC_indices_on_SPG( aMeshIndex, iSPGIndex );

                // loop over the unzipped lagrange cells and evaluate t-matrices
                for ( const auto& iUIPCIndex : tUIPCIndices )
                {
                    // get base cell of the lagrange mesh
                    Interpolation_Cell_Unzipped* tUIPC = tEnrichedIPCells( iUIPCIndex );

                    moris::Cell< moris::Cell< mtk::Vertex* > > tBsplineBasis;
                    moris::Cell< Matrix< DDRMat > >            tWeights;

                    // get the indices, ids and weights from the hmr mesh from the
                    mXTKModelPtr->mBackgroundMesh->get_extended_t_matrix( aMeshIndex, tRootBSplineCellIndex, *tUIPC->get_base_cell(), tBsplineBasis, tWeights );

                    // get a const pointer to UIPC such that the following xtk specific function can be called
                    const Interpolation_Cell_Unzipped* tConstUIPC = tUIPC;

                    // Get the lagrange vertices of the UIPC to obtain the vertex interpolation object of each lagrange node
                    moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tXTKUIPVs = tConstUIPC->get_xtk_interpolation_vertices();

                    // loop over the unzipped interpolation vertices and get their t-matrix info to overwrite
                    for ( size_t iLocalVertIndex = 0; iLocalVertIndex < tXTKUIPVs.size(); iLocalVertIndex++ )
                    {
                        // get the vertex enrichment object
                        xtk::Vertex_Enrichment* tVertexEnrichment = tXTKUIPVs( iLocalVertIndex )->get_xtk_interpolation( aMeshIndex );

                        // access the basis indices
                        moris::Matrix< moris::IndexMat > const & tBasisIndices = tVertexEnrichment->get_basis_function_indices();
                        moris::Matrix< moris::DDRMat >&          tBasisWeights = tVertexEnrichment->get_basis_function_weights();
                        IndexMap&                                tBasisMap     = tVertexEnrichment->get_basis_function_map();

                        // probably need to do a resize call here
                        // overallocation basis Indices
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisIndices( tBasisIndices.numel() + tBsplineBasis( iLocalVertIndex ).size() );
                        moris::Matrix< DDRMat >          tAgglomeratedBasisWeights( tBasisIndices.numel() + tBsplineBasis( iLocalVertIndex ).size() );
                        tBasisMap.clear();

                        // initialize a counter to count the number of basis that will be replaced
                        uint iBasisCounter = 0;

                        // Loop over the basis indices and check if there is a dependent basis
                        for ( size_t iLocalBasisIndex = 0; iLocalBasisIndex < tBasisIndices.numel(); iLocalBasisIndex++ )
                        {
                            moris_index iOldBasisIndex = tBasisIndices( iLocalBasisIndex );

                            // check if this basis is good basis or it does not belong to the current
                            if ( 0 == mBasisData( aMeshIndex ).mFollowerBasis( iOldBasisIndex ) )
                            {
                                tAgglomeratedBasisIndices( iBasisCounter ) = iOldBasisIndex;
                                tAgglomeratedBasisWeights( iBasisCounter ) = tBasisWeights( iLocalBasisIndex );
                                iBasisCounter++;
                            }
                            // if this is a bad basis place all the
                            else
                            {
                                // get the root SPG index
                                moris_index tRootSPGIndex = mRootSPGIds( iSPGIndex );

                                // This association checks such that 1 basis does not get modified twice
                                if ( mBasisData( aMeshIndex ).mFollowerBasisOwningCell( iOldBasisIndex ) != tRootSPGIndex )
                                {
                                    tAgglomeratedBasisIndices( iBasisCounter ) = iOldBasisIndex;
                                    tAgglomeratedBasisWeights( iBasisCounter ) = tBasisWeights( iLocalBasisIndex );
                                    iBasisCounter++;

                                    continue;
                                }

                                // This two maps help to get the the enriched basis index from the
                                moris::Cell< moris_index > const & tBGBasisIndices = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisIndices( tRootSPGIndex );
                                moris::Cell< moris_index > const & tBGBasisLevels  = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisEnrLev( tRootSPGIndex );

                                // convert the bg basis indices to the map
                                IndexMap tBasisToLocalIndexMap;
                                convert_cell_to_map( tBGBasisIndices, tBasisToLocalIndexMap );

                                // get the enrichment indices that contains combines data of the above cells
                                moris::Cell< Matrix< IndexMat > > const & tBasisEnrichmentIndices = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mBasisEnrichmentIndices;

                                // loop over the unenriched BSpline vertices
                                for ( uint iBGBasisIndex = 0; iBGBasisIndex < tBsplineBasis( iLocalVertIndex ).size(); iBGBasisIndex++ )
                                {
                                    // get the bg basis index
                                    moris_index tBGBasisIndex = tBsplineBasis( iLocalVertIndex )( iBGBasisIndex )->get_index();

                                    // get the local index of the bg basis index
                                    auto        tIter            = tBasisToLocalIndexMap.find( tBGBasisIndex );
                                    moris_index tLocalBasisIndex = tIter->second;

                                    // get the enrichment level of the bg basis index
                                    moris_index tEnrLevel = tBGBasisLevels( tLocalBasisIndex );

                                    // get the enriched index
                                    moris_index tEnrichedBasisIndex = tBasisEnrichmentIndices( tBGBasisIndex )( tEnrLevel );

                                    // get the local index of the enriched basis index
                                    tAgglomeratedBasisIndices( iBasisCounter ) = tEnrichedBasisIndex;
                                    tAgglomeratedBasisWeights( iBasisCounter ) = tWeights( iLocalVertIndex )( iBGBasisIndex );

                                    // increment the counter
                                    iBasisCounter++;
                                }
                            }
                        }

                        // resize the basis indices and weights and the basis ids based on the actual number of basis
                        tAgglomeratedBasisIndices.resize( iBasisCounter, 1 );
                        tAgglomeratedBasisWeights.resize( iBasisCounter, 1 );
                        Matrix< IndexMat > tAgglomeratedBasisIds( iBasisCounter, 1 );
                        // tBasisMap.reserve( iBasisCounter ); // FIXME: only works if the Mini_Map is implemented using a std::unordered_map

                        // get the map to convert indices to ids
                        Matrix< IndexMat > const & tBasisEnrichmentIndices = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisIndexToId;

                        // loop over the basis and get ids and add that to the basis map
                        for ( uint iBC = 0; iBC < tAgglomeratedBasisIndices.numel(); iBC++ )
                        {
                            moris::moris_index tBasisIndex = tAgglomeratedBasisIndices( iBC );
                            tBasisMap[ tBasisIndex ]       = iBC;
                            tAgglomeratedBasisIds( iBC )   = tBasisEnrichmentIndices( tBasisIndex );
                        }

                        // add vertex enrichment data
                        tVertexEnrichment->add_basis_information( tAgglomeratedBasisIndices, tAgglomeratedBasisIds );
                        tVertexEnrichment->add_basis_function_weights( tAgglomeratedBasisIndices, tAgglomeratedBasisWeights );
                    }
                }
            }
        }
    }

    //----------------------------------------------------------------------------------------------

    void
    Basis_Processor::perform_nearest_neighbor_exchange( moris_index aMeshIndex )
    {
        if ( par_size() == 1 )
        {
            return;
        }

        // get the communication table and map
        Matrix< IdMat > const & tCommTable = mXTKModelPtr->mCutIntegrationMesh->get_communication_table();

        // convert it to a cell
        moris::Cell< moris_index > tCommCell;
        tCommCell.insert( tCommCell.size(), tCommTable.begin(), tCommTable.end() );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Prepare requests for non-owned entities */

        // initialize lists of information that identifies entities (on other procs)
        Cell< Cell< moris_index > >     tNotOwnedSpgsToProcs;    // SPG index (local to current proc, just used for construction of arrays)
        Cell< moris::Cell< moris_id > > tSendSPGIds;             // first SP IDs in SPGs in each of the SPGs

        // TODO: this function could be further optimized by preparing data of the SPGs that don't have a root SPG
        this->prepare_requests_for_not_owned_subphase_groups( aMeshIndex, tNotOwnedSpgsToProcs, tSendSPGIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Send and Receive requests about non-owned entities to and from other procs */

        // initialize arrays for receiving
        Cell< Cell< moris_id > > tReceivedSPGIds;

        // communicate information
        moris::communicate_cells( tCommCell, tSendSPGIds, tReceivedSPGIds );

        // clear memory not needed anymore
        tSendSPGIds.clear();
        shrink_to_fit_all( tSendSPGIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Find answers to the requests */

        // initialize lists of ID answers to other procs
        Cell< moris::Cell< moris_id > > tSendSubphaseGroupRootIds;
        Cell< moris::Cell< moris_id > > tSendSubphaseGroupRootOwners;

        // answer requests from other procs
        this->prepare_answers_for_owned_subphase_groups( aMeshIndex,
                tSendSubphaseGroupRootIds,
                tSendSubphaseGroupRootOwners,
                tReceivedSPGIds );

        // clear memory from requests (the answers to which have been found)
        tReceivedSPGIds.clear();
        shrink_to_fit_all( tReceivedSPGIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 3: Send and receive answers to and from other procs */

        // initialize arrays for receiving
        Cell< moris::Cell< moris_id > > tReceivedSubphaseGroupRootIds;
        Cell< moris::Cell< moris_id > > tReceivedSubphaseGroupRootOwners;

        // communicate answers
        moris::communicate_cells( tCommCell, tSendSubphaseGroupRootIds, tReceivedSubphaseGroupRootIds );
        moris::communicate_cells( tCommCell, tSendSubphaseGroupRootOwners, tReceivedSubphaseGroupRootOwners );

        // clear unused memory
        tSendSubphaseGroupRootIds.clear();
        shrink_to_fit_all( tSendSubphaseGroupRootIds );
        /* ---------------------------------------------------------------------------------------- */
        /* Step 4: Use answers to assign IDs to non-owned entities */

        this->handle_requested_subphase_groups_answers(
                tReceivedSubphaseGroupRootIds,
                tReceivedSubphaseGroupRootOwners,
                tNotOwnedSpgsToProcs );
    }

    //----------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_requests_for_not_owned_subphase_groups(
            moris_index const &              aMeshIndex,
            Cell< Cell< moris_index > >&     aNotOwnedSpgsToProcs,
            Cell< moris::Cell< moris_id > >& aSPGIDs )
    {
        // get the bspline mesh info
        const Bspline_Mesh_Info*               tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< Subphase_Group* > const & tSPGs            = tBsplineMeshInfo->mSubphaseGroups;

        // get the communication table and map
        Matrix< IdMat > const &                   tCommTable              = mXTKModelPtr->mCutIntegrationMesh->get_communication_table();
        uint                                      tCommTableSize          = tCommTable.numel();
        std::map< moris_id, moris_index > const & tProcIdToCommTableIndex = mXTKModelPtr->mCutIntegrationMesh->get_communication_map();

        // get the non-owned Subphases on the executing processor
        moris::Cell< moris_index > const & tNotOwnedSPGs    = tBsplineMeshInfo->mNotOwnedSubphaseGroupIndices;
        uint                               tNumNotOwnedSPGs = tNotOwnedSPGs.size();

        // initialize lists of identifying information
        aNotOwnedSpgsToProcs.resize( tCommTableSize );
        aSPGIDs.resize( tCommTableSize );

        // reserve enough space for all non-owned SPGs  ( over allocation, but better than under allocation )
        std::for_each( aNotOwnedSpgsToProcs.begin(), aNotOwnedSpgsToProcs.end(), [ tNumNotOwnedSPGs ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumNotOwnedSPGs ); } );
        std::for_each( aSPGIDs.begin(), aSPGIDs.end(), [ tNumNotOwnedSPGs ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumNotOwnedSPGs ); } );

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( const auto iNotOwnedSPGIndex : tNotOwnedSPGs )
        {
            // ... get their respective owners, and position in the comm table ...
            moris_index tOwnerProc = tSPGs( iNotOwnedSPGIndex )->get_owner();
            auto        tIter      = tProcIdToCommTableIndex.find( tOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Basis_Processor::prepare_requests_for_not_owned_subphase_groups() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i which is: %s",
                    tOwnerProc,
                    par_rank(),
                    ios::stringify_log( tCommTable ).c_str() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SPGs in the list of SPs to be requested from that owning proc
            aNotOwnedSpgsToProcs( tProcDataIndex ).push_back( iNotOwnedSPGIndex );
        }

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned SPGs to be sent to each processor processor
            uint tNumNotOwnedSpgsOnProc = aNotOwnedSpgsToProcs( iProc ).size();

            // allocate matrix
            aSPGIDs( iProc ).resize( tNumNotOwnedSpgsOnProc );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iSPG = 0; iSPG < tNumNotOwnedSpgsOnProc; iSPG++ )
            {
                // get the index of the subphase group on the executing proc
                moris_index tSpgIndex = aNotOwnedSpgsToProcs( iProc )( iSPG );

                // store the identifying information of the Subphase group in the output arrays
                aSPGIDs( iProc )( iSPG ) = tSPGs( tSpgIndex )->get_id();
            }
        }

        // size out unused memory
        shrink_to_fit_all( aNotOwnedSpgsToProcs );
        shrink_to_fit_all( aSPGIDs );
    }

    //----------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_answers_for_owned_subphase_groups(
            moris_index const &                            aMeshIndex,
            moris::Cell< moris::Cell< moris_id > >&        aSendSubphaseGroupRootIds,
            moris::Cell< moris::Cell< moris_id > >&        aSendSubphaseGroupRootOwners,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSPGIds )
    {
        // get the bspline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        // moris::Cell<Subphase_Group*> const & tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the communication table and map
        Matrix< IdMat > const & tCommTable     = mXTKModelPtr->mCutIntegrationMesh->get_communication_table();
        uint                    tCommTableSize = tCommTable.numel();
        // std::map< moris_id, moris_index >& tProcIdToCommTableIndex = mXTKModelPtr->mCutIntegrationMesh->->get_communication_map();

        // initialize array of answers
        aSendSubphaseGroupRootIds.resize( tCommTableSize );
        aSendSubphaseGroupRootOwners.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedSPGIds.size() == tCommTableSize,
                "Basis_Processor::prepare_answers_for_owned_subphase_groups() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of entity IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedSPGIds( iProc ).size();

            // size the list of answers / IDs accordingly
            aSendSubphaseGroupRootIds( iProc ).resize( tNumReceivedReqs );
            aSendSubphaseGroupRootOwners( iProc ).resize( tNumReceivedReqs );

            // iterate through SPs for which the IDs are requested
            for ( uint iSPG = 0; iSPG < tNumReceivedReqs; iSPG++ )
            {
                // get the ID and index of the received SP
                moris_id    tSubphaseId    = aReceivedSPGIds( iProc )( iSPG );
                moris_index tSubphaseIndex = tBsplineMeshInfo->get_index_for_spg_id( tSubphaseId );

                // place Subphase ID in return data
                aSendSubphaseGroupRootIds( iProc )( iSPG )    = mRootSPGIds( tSubphaseIndex );
                aSendSubphaseGroupRootOwners( iProc )( iSPG ) = mRootSPGOwner( tSubphaseIndex );

            }    // end for: communication for each entity with current processor

        }        // end for: communication list for each processor
    }

    //----------------------------------------------------------------------------------------------

    void
    Basis_Processor::handle_requested_subphase_groups_answers(
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSubphaseGroupRootIds,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSubphaseGroupRootOwners,
            moris::Cell< moris::Cell< moris_id > > const & aNotOwnedSpgsToProcs )
    {
        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedSubphaseGroupRootIds.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedSpgIds = aReceivedSubphaseGroupRootIds( iProc ).size();

            // check that the
            MORIS_ASSERT( tNumReceivedSpgIds = aNotOwnedSpgsToProcs( iProc ).size(),
                    "Integration_Mesh_Generator::handle_requested_subphase_group_ID_answers() - "
                    "Number of SPG ID requests and answers are not the same." );

            // assign IDs to each communicated entity
            for ( uint iSPG = 0; iSPG < tNumReceivedSpgIds; iSPG++ )
            {
                // get the current SPG index and ID from the data provided
                moris_index tSubphaseGroupIndex     = aNotOwnedSpgsToProcs( iProc )( iSPG );
                moris_id    tSubphaseGroupRootId    = aReceivedSubphaseGroupRootIds( iProc )( iSPG );
                moris_id    tSubphaseGroupRootOwner = aReceivedSubphaseGroupRootOwners( iProc )( iSPG );

                // store the received SPG ID
                mRootSPGIds( tSubphaseGroupIndex )   = tSubphaseGroupRootId;
                mRootSPGOwner( tSubphaseGroupIndex ) = tSubphaseGroupRootOwner;
            }
        }
    }

    //----------------------------------------------------------------------------------------------

    bool
    Basis_Processor::determine_stopping_criteria()
    {
        // contains negative index
        bool tContainsNoIndex = std::find( mRootSPGIds.begin(), mRootSPGIds.end(), gNoIndex ) != mRootSPGIds.end();

        if ( par_size() == 1 )
        {
            return tContainsNoIndex;
        }

        // determine if all processors have no index
        bool tContainsNoIndexAllProcs = all_lor( tContainsNoIndex );

        // return true if all processors have no index
        return tContainsNoIndexAllProcs;
    }

    //----------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_comm_table()
    {
        moris_id tSize = par_size();
        moris_id tRank = par_rank();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Create a cell that shows what are the list of processor that current processor needs to receive data  */

        // create an unordered set to count the unique processors that current processor needs to receive dat
        std::unordered_set< moris_index > tProcessorsToReceiveFromSet( mRootSPGOwner.begin(), mRootSPGOwner.end() );

        // in order to communicate the set above it needs to be in a continues memory block ( cell )
        moris::Cell< moris_id > tProcessorsToReceiveFrom;
        tProcessorsToReceiveFrom.reserve( tSize );

        // copy the set to the cell
        std::copy( tProcessorsToReceiveFromSet.begin(), tProcessorsToReceiveFromSet.end(), std::back_inserter( tProcessorsToReceiveFrom ) );

        // resize the remaining cell member with gNoIndex
        tProcessorsToReceiveFrom.resize( tSize, gNoIndex );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: choose a root processor that will receive the data from all other processors about the neighborhood and will establish two-way communion tables  */

        moris_index tRootProc = 0;

        // create a continues memory array that will contain all the relationship between processors
        moris::Cell< moris_id > tProcessorsToReceiveFromAllProcs( tSize * tSize );
        if ( tRootProc == par_rank() )
        {
            // gather everything on processor zero
            MPI_Gather( tProcessorsToReceiveFrom.memptr(), tSize, MPI_INT, tProcessorsToReceiveFromAllProcs.memptr(), tSize, MPI_INT, tRootProc, MPI_COMM_WORLD );
        }
        else
        {
            // send the data to the root processor
            MPI_Gather( tProcessorsToReceiveFrom.memptr(), tSize, MPI_INT, NULL, 0, MPI_INT, tRootProc, MPI_COMM_WORLD );
        }

        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: process the data that will need to be distrusted to individual processors as the comm table  */
        // create a continues memory array that will contain all the relationship between processors
        moris::Cell< moris::Cell< moris_id > > tCommTable( tSize );
        if ( tRootProc == par_rank() )
        {
            for ( auto& iCell : tCommTable )
            {
                iCell.resize( tSize, gNoIndex );
            }

            // iterate through the data and create the comm table
            for ( moris_index iProcIndex = 0; iProcIndex < tSize; iProcIndex++ )
            {
                for ( moris_index iProcToReceiveFromIndex = 0; iProcToReceiveFromIndex < tSize; iProcToReceiveFromIndex++ )
                {
                    moris_index iProcToReceiveFrom = tProcessorsToReceiveFromAllProcs( iProcIndex * tSize + iProcToReceiveFromIndex );

                    // if the processor is not gNoIndex then add it to the comm table
                    if ( iProcToReceiveFrom != gNoIndex )
                    {
                        // add it to the comm table
                        tCommTable( iProcIndex )( iProcToReceiveFrom ) = 1;

                        // add the reverse relationship to the comm table
                        tCommTable( iProcToReceiveFrom )( iProcIndex ) = 1;
                    }
                }
            }
        }

        /* ---------------------------------------------------------------------------------------- */
        /* Step 3: scatter the data to their corresponding  processor*/

        mCommTable.resize( tSize );

        if ( tRootProc == par_rank() )
        {
            // create a continues memory array that will contain all the relationship between processors
            moris::Cell< moris_id > tCommTableFlattened;
            tCommTableFlattened.resize( tSize * tSize );

            // iterate through the data and create the comm table
            for ( moris_index iProcIndex = 0; iProcIndex < tSize; iProcIndex++ )
            {
                for ( moris_index iProcToReceiveFromIndex = 0; iProcToReceiveFromIndex < tSize; iProcToReceiveFromIndex++ )
                {
                    moris_index iNeighbors = tCommTable( iProcIndex )( iProcToReceiveFromIndex );

                    // if the processor is not gNoIndex then add it to the comm table
                    if ( iNeighbors != gNoIndex )
                    {
                        // add it to the comm table
                        tCommTableFlattened( iProcIndex * tSize + iProcToReceiveFromIndex ) = iProcToReceiveFromIndex;
                    }
                    else
                    {
                        // add it to the comm table
                        tCommTableFlattened( iProcIndex * tSize + iProcToReceiveFromIndex ) = gNoIndex;
                    }
                }
            }

            MPI_Scatter( tCommTableFlattened.memptr(), tSize, MPI_INT, mCommTable.memptr(), tSize, MPI_INT, tRootProc, MPI_COMM_WORLD );
        }
        else
        {
            MPI_Scatter( NULL, tSize, MPI_INT, mCommTable.memptr(), tSize, MPI_INT, tRootProc, MPI_COMM_WORLD );
        }

        /* ---------------------------------------------------------------------------------------- */
        /* Step 4: remove negative values from the mCommTable*/
        // Remove negative values using erase-remove idiom

        auto tIsNegative = [ tRank ]( int aNum ) { return aNum < 0 || aNum == tRank; };
        mCommTable.data().erase( std::remove_if( mCommTable.begin(), mCommTable.end(), tIsNegative ), mCommTable.end() );
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_to_leader_basis_weights_indices_neighbors( moris_index aMeshIndex )
    {
        if ( par_size() == 1 )
        {
            return;
        }

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Prepare requests for non-owned entities */

        // initialize lists of information that identifies entities (on other procs)
        Cell< Cell< moris_index > >     tNotOwnedSpgsToProcs;    // SPG index (local to current proc, just used for construction of arrays)
        Cell< moris::Cell< moris_id > > tSendSPGIds;             // first SP IDs in SPGs in each of the SPGs

        // TODO: this function could be further optimized by preparing data of the SPGs that don't have a root SPG
        this->prepare_requests_for_not_owned_root_spg_projections( aMeshIndex, tNotOwnedSpgsToProcs, tSendSPGIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Send and Receive requests about non-owned entities to and from other procs */

        // initialize arrays for receiving
        Cell< Cell< moris_id > > tReceivedSPGIds;

        // communicate information
        moris::communicate_cells( mCommTable, tSendSPGIds, tReceivedSPGIds );

        // clear memory not needed anymore
        tSendSPGIds.clear();
        shrink_to_fit_all( tSendSPGIds );
        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Find answers to the requests */

        // initialize lists of ID answers to other procs
        Cell< moris::Cell< moris_id > > tSendSubphaseGroupRootBsplineElementIds;
        Cell< moris::Cell< moris_id > > tSendSubphaseGroupRootBsplineBasisIds;
        Cell< moris::Cell< moris_id > > tSendSubphaseGroupRootBsplineBasisOwners;

        // answer requests from other procs
        this->prepare_answers_for_owned_root_spg_projections( aMeshIndex,
                tSendSubphaseGroupRootBsplineElementIds,
                tSendSubphaseGroupRootBsplineBasisIds,
                tSendSubphaseGroupRootBsplineBasisOwners,
                tReceivedSPGIds );

        // clear memory from requests (the answers to which have been found)
        tReceivedSPGIds.clear();
        shrink_to_fit_all( tReceivedSPGIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 3: Send and receive answers to and from other procs */

        // initialize arrays for receiving
        Cell< moris::Cell< moris_id > > tReceivedSubphaseGroupRootBsplineElementIds;
        Cell< moris::Cell< moris_id > > tReceivedSubphaseGroupRootBsplineBasisIds;
        Cell< moris::Cell< moris_id > > tReceivedSubphaseGroupRootBsplineBasisOwners;

        // communicate answers
        moris::communicate_cells( mCommTable, tSendSubphaseGroupRootBsplineElementIds, tReceivedSubphaseGroupRootBsplineElementIds );
        moris::communicate_cells( mCommTable, tSendSubphaseGroupRootBsplineBasisIds, tReceivedSubphaseGroupRootBsplineBasisIds );
        moris::communicate_cells( mCommTable, tSendSubphaseGroupRootBsplineBasisOwners, tReceivedSubphaseGroupRootBsplineBasisOwners );

        // clear unused memory
        tSendSubphaseGroupRootBsplineElementIds.clear();
        tSendSubphaseGroupRootBsplineBasisIds.clear();
        tSendSubphaseGroupRootBsplineBasisOwners.clear();
        shrink_to_fit_all( tSendSubphaseGroupRootBsplineElementIds );
        shrink_to_fit_all( tSendSubphaseGroupRootBsplineBasisIds );
        shrink_to_fit_all( tSendSubphaseGroupRootBsplineBasisOwners );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 4: Use answers to assign IDs to non-owned entities */

        this->handle_requested_root_spg_projections( aMeshIndex,
                tReceivedSubphaseGroupRootBsplineElementIds,
                tReceivedSubphaseGroupRootBsplineBasisIds,
                tReceivedSubphaseGroupRootBsplineBasisOwners,
                tNotOwnedSpgsToProcs );
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_requests_for_not_owned_root_spg_projections( moris_index aMeshIndex,
            Cell< Cell< moris_index > >&                                              aNotOwnedSpgsToProcs,
            Cell< moris::Cell< moris_id > >&                                          aSPGIDs )
    {
        // get the bspline mesh info
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        // moris::Cell< Subphase_Group* > const & tSPGs            = tBsplineMeshInfo->mSubphaseGroups;

        // get the communication table and map
        uint                                        tCommTableSize = mCommTable.size();
        std::unordered_map< moris_id, moris_index > tProcIdToCommTableIndex;
        tProcIdToCommTableIndex.reserve( tCommTableSize );

        // create  a map from proc id to comm table index
        for ( moris_index iIndex = 0; iIndex < (moris_index)tCommTableSize; ++iIndex )
        {
            tProcIdToCommTableIndex[ mCommTable( iIndex ) ] = iIndex;
        }

        // get the non-owned Subphases on the executing processor
        moris::Cell< moris_index >& tOwnedSPGs    = tBsplineMeshInfo->mOwnedSubphaseGroupIndices;
        uint                        tNumOwnedSPGs = tOwnedSPGs.size();

        // initialize lists of identifying information
        aNotOwnedSpgsToProcs.resize( tCommTableSize );
        aSPGIDs.resize( tCommTableSize );

        // reserve enough space for all non-owned SPGs  ( over allocation, but better than under allocation )
        std::for_each( aNotOwnedSpgsToProcs.begin(), aNotOwnedSpgsToProcs.end(), [ tNumOwnedSPGs ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumOwnedSPGs ); } );
        std::for_each( aSPGIDs.begin(), aSPGIDs.end(), [ tNumOwnedSPGs ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumOwnedSPGs ); } );

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( const auto& iNotOwnedSPGIndex : tOwnedSPGs )
        {
            // ... for the non-owned spgs get the root and determine if the root is on the same processor
            // ... if not add it to the list
            moris_index tRootOwnerProc = mRootSPGOwner( iNotOwnedSPGIndex );

            // if the root is on the same proc, then skip it
            if ( tRootOwnerProc == par_rank() )
            {
                continue;
            }

            // otherwise add the root to the list of SPs to be requested from that owning proc
            auto tIter = tProcIdToCommTableIndex.find( tRootOwnerProc );
            MORIS_ASSERT(
                    tIter != tProcIdToCommTableIndex.end(),
                    "Integration_Mesh_Generator::prepare_requests_for_not_owned_subphase_group_IDs() - "
                    "Entity owner (Proc #%i) not found in communication table of current proc #%i when preparing",
                    tRootOwnerProc,
                    par_rank() );
            moris_index tProcDataIndex = tIter->second;

            // ... and finally add the non-owned SPGs in the list of SPs to be requested from that owning proc
            aNotOwnedSpgsToProcs( tProcDataIndex ).push_back( iNotOwnedSPGIndex );
        }

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of non-owned SPGs to be sent to each processor processor
            uint tNumNotOwnedSpgsOnProc = aNotOwnedSpgsToProcs( iProc ).size();

            // allocate matrix
            aSPGIDs( iProc ).resize( tNumNotOwnedSpgsOnProc );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iSPG = 0; iSPG < tNumNotOwnedSpgsOnProc; iSPG++ )
            {
                // get the index of the subphase group on the executing proc
                moris_index tSpgIndex = aNotOwnedSpgsToProcs( iProc )( iSPG );

                // store the identifying information of the Subphase group in the output arrays
                aSPGIDs( iProc )( iSPG ) = mRootSPGIds( tSpgIndex );
            }
        }

        // size out unused memory
        shrink_to_fit_all( aNotOwnedSpgsToProcs );
        shrink_to_fit_all( aSPGIDs );
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::handle_requested_root_spg_projections(
            moris_index                                    aMeshIndex,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSubphaseGroupRootBsplineElementIds,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSubphaseGroupRootBsplineBasisIds,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedSubphaseGroupRootBsplineBasisOwners,
            moris::Cell< moris::Cell< moris_id > > const & aNotOwnedSpgsToProcs )
    {

        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedSubphaseGroupRootBsplineElementIds.size(); iProc++ )
        {
            // get the number of requests and answers from the current proc
            uint tNumReceivedSpgIds = aReceivedSubphaseGroupRootBsplineElementIds( iProc ).size();

            // check that the
            MORIS_ASSERT( tNumReceivedSpgIds == aNotOwnedSpgsToProcs( iProc ).size(),
                    "Integration_Mesh_Generator::handle_requested_subphase_group_ID_answers() - "
                    "Number of SPG ID requests and answers are not the same." );

            uint tStartRange = 0;

            // assign IDs to each communicated entity
            for ( uint iSPG = 0; iSPG < tNumReceivedSpgIds; iSPG++ )
            {
                // get the current SPG index and ID from the data provided
                moris_index                     tSubphaseGroupIndex     = aNotOwnedSpgsToProcs( iProc )( iSPG );
                moris_id                        tRootBsplineId          = aReceivedSubphaseGroupRootBsplineElementIds( iProc )( iSPG );
                moris::Cell< moris_id > const & tRootBsplineBasisId     = aReceivedSubphaseGroupRootBsplineBasisIds( iProc );
                moris::Cell< moris_id > const & tRootBsplineBasisOwners = aReceivedSubphaseGroupRootBsplineBasisOwners( iProc );

                // construct the extension for the current SPG
                this->construct_follower_to_leader_basis_relationship_for_spg(
                        tSubphaseGroupIndex,
                        tRootBsplineId,
                        tRootBsplineBasisId,
                        tRootBsplineBasisOwners,
                        tStartRange,
                        aMeshIndex );

                tStartRange = tStartRange + mHMRHelper( aMeshIndex )->get_number_of_bases_per_element();
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_answers_for_owned_root_spg_projections( moris_index const & aMeshIndex,
            moris::Cell< moris::Cell< moris_id > >&                                      aSendSubphaseGroupRootBsplineElementIds,
            moris::Cell< moris::Cell< moris_id > >&                                      aSendSubphaseGroupRootBsplineBasisIds,
            moris::Cell< moris::Cell< moris_id > >&                                      aSendSubphaseGroupRootBsplineBasisOwners,
            moris::Cell< moris::Cell< moris_id > > const &                               aReceivedSPGIds )
    {
        // get the bspline mesh info
        Bspline_Mesh_Info*                     tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< mtk::Cell* > const &      tSPGBsplineCells = tBsplineMeshInfo->mBsplineCells;
        moris::Cell< Subphase_Group* > const & tSPGs            = tBsplineMeshInfo->mSubphaseGroups;
        uint                                   tNumberOfBasis   = mHMRHelper( aMeshIndex )->get_number_of_bases_per_element();

        // get the communication table and map
        uint tCommTableSize = mCommTable.size();
        // std::map< moris_id, moris_index >& tProcIdToCommTableIndex = mXTKModelPtr->mCutIntegrationMesh->->get_communication_map();

        // initialize array of answers
        aSendSubphaseGroupRootBsplineElementIds.resize( tCommTableSize );
        aSendSubphaseGroupRootBsplineBasisIds.resize( tCommTableSize );
        aSendSubphaseGroupRootBsplineBasisOwners.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedSPGIds.size() == tCommTableSize,
                "Basis_Processor::prepare_answers_for_owned_subphase_groups() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of entity IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedSPGIds( iProc ).size();

            // size the list of answers / IDs accordingly
            aSendSubphaseGroupRootBsplineElementIds( iProc ).resize( tNumReceivedReqs );
            aSendSubphaseGroupRootBsplineBasisIds( iProc ).reserve( tNumReceivedReqs * tNumberOfBasis );
            aSendSubphaseGroupRootBsplineBasisOwners( iProc ).reserve( tNumReceivedReqs * tNumberOfBasis );

            // iterate through SPs for which the IDs are requested
            for ( uint iSPG = 0; iSPG < tNumReceivedReqs; iSPG++ )
            {
                // get the ID and index of the received SP
                moris_id    tSubphaseId    = aReceivedSPGIds( iProc )( iSPG );
                moris_index tSubphaseIndex = tBsplineMeshInfo->get_index_for_spg_id( tSubphaseId );

                // get the bspline cell index based on the subphase index of the SPG
                moris_index tBsplineCellIndex = tSPGs( tSubphaseIndex )->get_bspline_cell_index();
                mtk::Cell*  tBsplineCell      = tSPGBsplineCells( tBsplineCellIndex );

                // HMR interface to get the global ID of the bspline cell
                luint tGlobalElementID = mHMRHelper( aMeshIndex )->get_global_domain_id_of_cell( tBsplineCell );

                // assign the first data that needs to be passed on to the requesting processor
                aSendSubphaseGroupRootBsplineElementIds( iProc )( iSPG ) = tGlobalElementID;

                // get the enriched basis of the bspline cell based on the enrichment level
                moris::Cell< moris_id > const & tEnrichedBasisIds    = mHMRHelper( aMeshIndex )->get_enriched_basis_id_of_cell( tBsplineCell, tSubphaseIndex );
                moris::Cell< moris_id > const & tEnrichedBasisOwners = mHMRHelper( aMeshIndex )->get_enriched_basis_owners_of_cell();

                // insert the data in the appropriate arrays
                aSendSubphaseGroupRootBsplineBasisIds( iProc )
                        .insert( aSendSubphaseGroupRootBsplineBasisIds( iProc ).size(),    //
                                tEnrichedBasisIds.begin(),
                                tEnrichedBasisIds.end() );

                // insert the data in the appropriate arrays
                aSendSubphaseGroupRootBsplineBasisOwners( iProc )
                        .insert( aSendSubphaseGroupRootBsplineBasisOwners( iProc ).size(),    //
                                tEnrichedBasisOwners.begin(),
                                tEnrichedBasisOwners.end() );

            }    // end for: communication for each entity with current processor

        }        // end for: communication list for each processor
    }


    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_to_leader_basis_relationship_for_spg(
            moris_index                     aSubphaseGroupIndex,
            moris_id                        aRootBsplineId,
            moris::Cell< moris_id > const & aRootBsplineBasisId,
            moris::Cell< moris_id > const & aRootBsplineBasisOwner,
            uint                            aStartRange,
            uint                            aMeshIndex )
    {
        // Get a reference or a pointer to the required data
        Bspline_Mesh_Info*                                tBsplineMeshInfo              = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        moris::Cell< Subphase_Group* > const &            tSPGs                         = tBsplineMeshInfo->mSubphaseGroups;
        moris::Cell< moris_index > const &                tNonEnrBfIndForEnrBfInd       = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mNonEnrBfIndForEnrBfInd;
        Enriched_Interpolation_Mesh&                      tEnrInterpMesh                = mXTKModelPtr->get_enriched_interp_mesh( aMeshIndex );

        std::unordered_map< moris_index, moris::Cell< real > >& tAveragingWeightsSPGBased = mBasisData( aMeshIndex ).mAveragingWeightsSPGBased;

        // get the enriched basis function that might be bad in the extened cell ( old ones )
        moris::Cell< moris_index > const & tEnrichedBFInSPG = tEnrichedBasisInSubphaseGroup( aSubphaseGroupIndex );

        // get the Bspline cell corresponding to the extened SPG index
        moris_index tExtenedBspCellIndex = tSPGs( aSubphaseGroupIndex )->get_bspline_cell_index();
        mtk::Cell*  tExtendedBspCell     = tBsplineMeshInfo->mBsplineCells( tExtenedBspCellIndex );

        // get the L2 projection matrix for the extended cell
        Matrix< DDRMat > const & tL2Projection = mHMRHelper( aMeshIndex )->get_l2_projection_matrix( tExtendedBspCell, aRootBsplineId );

        // get the list of BG Basis function Indices that exist on the extened SPG index
        moris::Cell< moris_index > const & tExtenedBGBspIndices = mHMRHelper( aMeshIndex )->get_bg_basis_indices_of_cell( tExtendedBspCell );

        // loop over the enriched basis functions in the SPG
        for ( uint iEnrichedBasisOrd = 0; iEnrichedBasisOrd < tEnrichedBFInSPG.size(); iEnrichedBasisOrd++ )
        {
            moris_index const & iEnrichedBasis = tEnrichedBFInSPG( iEnrichedBasisOrd );

            // if it is a good basis jump to the next one
            if ( mBasisData( aMeshIndex ).mFollowerBasis( iEnrichedBasis ) == 0 )
            {
                continue;
            }

            // get the averaging weights for the enriched basis
            real const & tAveragingWeights = tAveragingWeightsSPGBased[ aSubphaseGroupIndex ]( iEnrichedBasisOrd );

            //  else statement( to prevent further indentation it is written without else statement)
            //  if the enriched basis is follower, then first convert it to bg basis index
            moris_index tBGExtenedBasisIndex = tNonEnrBfIndForEnrBfInd( iEnrichedBasis );

            // then find the location(ordinal) of this basis in the ordered bg basis indices
            auto tIterator = std::find_if( tExtenedBGBspIndices.begin(), tExtenedBGBspIndices.end(),    //
                    [ &tBGExtenedBasisIndex ]( moris_index aBasisIndex ) { return aBasisIndex == tBGExtenedBasisIndex; } );

            // find the ordinal of the non-enriched basis index in the extended basis
            uint tBGExtenedBasisOrdinal = std::distance( tExtenedBGBspIndices.begin(), tIterator );

            // resize the extened to root basis relationship to proper size
            mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iEnrichedBasis ).reserve( mHMRHelper( aMeshIndex )->get_number_of_bases_per_element() );
            mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iEnrichedBasis ).reserve( mHMRHelper( aMeshIndex )->get_number_of_bases_per_element() );
            mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( iEnrichedBasis ).reserve( mHMRHelper( aMeshIndex )->get_number_of_bases_per_element() );

            // loop over the all the root basis and if the projection coefficient is not zero, then add it to the follower to leader basis relationship
            for ( uint iRootBasisOrd = 0; iRootBasisOrd < tL2Projection.n_cols(); iRootBasisOrd++ )
            {
                // get the projection coefficient
                real tProjectionCoefficient = tL2Projection( tBGExtenedBasisOrdinal, iRootBasisOrd );

                //  if the projection coefficient is not zero, then add it to the follower to leader basis relationship
                if ( std::abs( tProjectionCoefficient ) > MORIS_REAL_EPS )
                {
                    // get the root basis index from the list
                    moris_index tRootEnrichedBasisId = aRootBsplineBasisId( iRootBasisOrd + aStartRange );

                    moris_id tRootEnrichedBasisOwners = aRootBsplineBasisOwner( iRootBasisOrd + aStartRange );

                    // If the root basis does not exist in the current processor add it to the list and and generate a new index for it
                    // add this basis to the mesh if it does not exists on the current partition
                    if ( !tEnrInterpMesh.basis_exists_on_partition( aMeshIndex, tRootEnrichedBasisId ) )
                    {
                        // get the bulk-phase the basis interpolates into
                        moris_index tBulkPhase = tSPGs( aSubphaseGroupIndex )->get_bulk_phase();

                        // add basis ID to partition
                        tEnrInterpMesh.add_basis_function( aMeshIndex, tRootEnrichedBasisId, tRootEnrichedBasisOwners, tBulkPhase );
                    }

                    // find and store the basis index local to the executing processor
                    moris_index tRootEnrichedBasisIndex = tEnrInterpMesh.get_enr_basis_index_from_enr_basis_id( aMeshIndex, tRootEnrichedBasisId );
                    moris_id    tOwner                  = tEnrInterpMesh.get_entity_owner( tRootEnrichedBasisIndex, mtk::EntityRank::BSPLINE, aMeshIndex );

                    //  add the root basis index to the follower to leader basis relationship
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iEnrichedBasis ).push_back( tRootEnrichedBasisIndex );

                    // add the projection coefficient to the follower to leader basis relationship
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iEnrichedBasis ).push_back( tAveragingWeights * tProjectionCoefficient );
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( iEnrichedBasis ).push_back( tOwner );
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_basis_using_volume( moris_index aMeshIndex )
    {
        // Get a reference or a pointer to the required data
        Bspline_Mesh_Info*                                tBsplineMeshInfo              = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        moris::Cell< Subphase_Group* > const &            tSPGs                         = tBsplineMeshInfo->mSubphaseGroups;


        // loop over the SPGs and determine if they are their own SPG , if not look for the neighbor, if the neighbor is complete then use it
        for ( uint iSPGIndex = 0; iSPGIndex < tSPGs.size(); iSPGIndex++ )
        {
            const Subphase_Group* tSubphaseGroup = tSPGs( iSPGIndex );
            // then all the basis attached to this are good basis
            if ( tSubphaseGroup->get_id() == mRootSPGIds( iSPGIndex ) )
            {
                moris::Cell< moris_index > const & tEnrBFsInInSPG = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                for ( const auto iEnrBF : tEnrBFsInInSPG )
                {
                    mBasisData( aMeshIndex ).mFollowerBasis( iEnrBF ) = 0;
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------------

    void
    Basis_Processor::construct_follower_to_leader_basis_weights_indices_mine( moris_index aMeshIndex )
    {
        // Get a reference or a pointer to the required data
        Bspline_Mesh_Info*                                tBsplineMeshInfo              = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        // Cell< moris::Matrix< IndexMat > > const & tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
        moris::Cell< Subphase_Group* > const &                  tSPGs                     = tBsplineMeshInfo->mSubphaseGroups;
        moris::Cell< moris_index > const &                      tNonEnrBfIndForEnrBfInd   = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mNonEnrBfIndForEnrBfInd;
        std::unordered_map< moris_index, moris::Cell< real > >& tAveragingWeightsSPGBased = mBasisData( aMeshIndex ).mAveragingWeightsSPGBased;


        // loop over the SPGs and find out the ones that need extension such that the root processor is within processor domain
        for ( uint iSPGIndex = 0; iSPGIndex < tSPGs.size(); iSPGIndex++ )
        {
            const Subphase_Group* tSubphaseGroup = tSPGs( iSPGIndex );

            // if the SPG is not on the current processor skip it
            if ( tSubphaseGroup->get_owner() != par_rank() )
            {
                continue;
            }
            // this means the SPG does not need extension,  SPGs root is itself, thus skip it
            if ( tSubphaseGroup->get_id() == mRootSPGIds( iSPGIndex ) )
            {
                continue;
            }
            // else case, it requires extension, check if the root SPG for it is on the processor domain, if the root is another processor skip this step
            if ( !tBsplineMeshInfo->spg_exists_on_partition( mRootSPGIds( iSPGIndex ) ) )
            {
                continue;
            }
            else
            {
                // even if it is on the processor domain, check if the root is owned by the current processor, if not skip this step
                if ( tSPGs( tBsplineMeshInfo->get_index_for_spg_id( mRootSPGIds( iSPGIndex ) ) )->get_owner() != par_rank() )
                {
                    continue;
                }
            }

            // else case, spg needs an extension and the root is owned on the processor domain
            //  step 1: construct the L2 projection matrix based on the b-spline cell
            moris_index tExtendedBsplineCellIndex = tSubphaseGroup->get_bspline_cell_index();
            moris_index tRootSPGIndex             = tBsplineMeshInfo->get_index_for_spg_id( mRootSPGIds( iSPGIndex ) );
            moris_index tRootBsplineCellIndex     = tSPGs( tRootSPGIndex )->get_bspline_cell_index();

            // get the Bspline cells corresponding to those indices
            mtk::Cell* tExtendedBsplineCell = tBsplineMeshInfo->mBsplineCells( tExtendedBsplineCellIndex );
            mtk::Cell* tRootBsplineCell     = tBsplineMeshInfo->mBsplineCells( tRootBsplineCellIndex );

            Matrix< DDRMat > const & tL2projectionMatrix = mHMRHelper( aMeshIndex )->get_l2_projection_matrix( tExtendedBsplineCell, tRootBsplineCell );

            // step 2: get the enriched basis
            //  TODO: the first one is copied, can add a member data to avoid copying
            moris::Cell< moris_index >         tEnrichedRootBasisIndices = mHMRHelper( aMeshIndex )->get_enriched_basis_indicies_of_cell( tRootBsplineCell, tRootSPGIndex );
            moris::Cell< moris_index > const & tBGExtendedBasisIndices   = mHMRHelper( aMeshIndex )->get_bg_basis_indices_of_cell( tExtendedBsplineCell );

            // step 3: get the enriched basis connected to the extended cell
            moris::Cell< moris_index > const & tEnrBFinSPGExtended = tEnrichedBasisInSubphaseGroup( iSPGIndex );

            // step 4: loop over the enriched basis in the extened call and find out which ones need to be replaces
            for ( uint iEnrichedBasisOrd = 0; iEnrichedBasisOrd < tEnrBFinSPGExtended.size(); iEnrichedBasisOrd++ )
            {
                moris_index iEneBFExtended = tEnrBFinSPGExtended( iEnrichedBasisOrd );

                // if it is not a follower basis, skip it
                if ( mBasisData( aMeshIndex ).mFollowerBasis( iEneBFExtended ) == 0 )
                {
                    continue;
                }

                // get the averaging weights for the enriched basis
                real tAveragingWeights = tAveragingWeightsSPGBased[ iSPGIndex ]( iEnrichedBasisOrd );

                // else condition if it is a follower then find the BG basis and the ordinal of the basis in the BG basis
                moris_index tBGBasisExtended = tNonEnrBfIndForEnrBfInd( iEneBFExtended );

                // find the ordinal of the basis in the BG Basis
                auto tIterator = std::find( tBGExtendedBasisIndices.begin(), tBGExtendedBasisIndices.end(), tBGBasisExtended );    //

                // find the ordinal of the non-enriched basis index in the extended basis
                uint tBGExtendedBasisOrd = std::distance( tBGExtendedBasisIndices.begin(), tIterator );

                // loop over the L2 projection matrix columns for the specific ordinal and find the one that is non-zero
                for ( uint iEnrRootBasisOrd = 0; iEnrRootBasisOrd < tL2projectionMatrix.n_cols(); iEnrRootBasisOrd++ )
                {
                    // get the projection coefficient
                    real tProjectionCoefficient = tL2projectionMatrix( tBGExtendedBasisOrd, iEnrRootBasisOrd );

                    // if the entry is non-zero, then it is the follower basis
                    if ( std::abs( tProjectionCoefficient ) > MORIS_REAL_EPS )
                    {
                        // find the non-enriched basis index in the root cell
                        moris_index tEnrRootBasis = tEnrichedRootBasisIndices( iEnrRootBasisOrd );

                        moris_id tOwner = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_entity_owner( tEnrRootBasis, mtk::EntityRank::BSPLINE, aMeshIndex );

                        // add the root basis index to the follower to leader basis relationship
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iEneBFExtended ).push_back( tEnrRootBasis );

                        // add the projection coefficient to the follower to leader basis relationship
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iEneBFExtended ).push_back( tAveragingWeights * tProjectionCoefficient );
                        mBasisData( aMeshIndex ).mFollowerToLeaderBasisOwners( iEneBFExtended ).push_back( tOwner );
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------------

    void
    Basis_Processor::communicate_shared_basis( moris_index aMeshIndex )
    {
        if ( par_size() == 1 )
        {
            return;
        }

        // get the communication table and map
        Matrix< IdMat > const & tCommTable = mXTKModelPtr->mCutIntegrationMesh->get_communication_table();

        // convert it to a cell
        moris::Cell< moris_index > tCommCell;
        tCommCell.insert( tCommCell.size(), tCommTable.begin(), tCommTable.end() );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 0: Prepare requests for non-owned entities */

        // initialize lists of information that identifies entities (on other procs)
        Cell< Cell< moris_index > >     tBasisIndexToProcs;    // SPG index (local to current proc, just used for construction of arrays)
        Cell< moris::Cell< moris_id > > tSendBasisIds;         // first SP IDs in SPGs in each of the SPGs

        // TODO: this function could be further optimized by preparing data of the SPGs that don't have a root SPG
        this->prepare_requests_for_shared_follower_basis( aMeshIndex, tBasisIndexToProcs, tSendBasisIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 1: Send and Receive requests about non-owned entities to and from other procs */

        // initialize arrays for receiving
        Cell< Cell< moris_id > > tReceivedBasisIds;

        // communicate information
        moris::communicate_cells( tCommCell, tSendBasisIds, tReceivedBasisIds );

        //  clear memory not needed anymore
        tSendBasisIds.clear();
        shrink_to_fit_all( tSendBasisIds );
        /* ---------------------------------------------------------------------------------------- */
        /* Step 2: Find answers to the requests */

        // initialize lists of ID answers to other procs
        Cell< moris::Cell< real > >     tSendFollowerToLeaderBasisWeights;
        Cell< moris::Cell< moris_id > > tSendFollowerToLeaderBasisIds;
        Cell< moris::Cell< moris_id > > tSendFollowerToLeaderBasisOwners;
        Cell< moris::Cell< moris_id > > tSendFollowerToLeaderOffset;

        // answer requests from other procs
        this->prepare_answers_for_follower_shared_basis( aMeshIndex,
                tSendFollowerToLeaderBasisIds,
                tSendFollowerToLeaderBasisOwners,
                tSendFollowerToLeaderBasisWeights,
                tSendFollowerToLeaderOffset,
                tReceivedBasisIds );


        // clear memory from requests (the answers to which have been found)
        tReceivedBasisIds.clear();
        shrink_to_fit_all( tReceivedBasisIds );

        /* ---------------------------------------------------------------------------------------- */
        /* Step 3: Send and receive answers to and from other procs */

        // initialize arrays for receiving
        Cell< moris::Cell< moris_id > > tReceivedFollowerToLeaderBasisIds;
        Cell< moris::Cell< moris_id > > tReceivedFollowerToLeaderBasisOwners;
        Cell< moris::Cell< real > >     tReceivedFollowerToLeaderBasisWeights;
        Cell< moris::Cell< moris_id > > tReceivedFollowerToLeaderOffset;

        // communicate answers
        moris::communicate_cells( tCommCell, tSendFollowerToLeaderBasisIds, tReceivedFollowerToLeaderBasisIds );
        moris::communicate_cells( tCommCell, tSendFollowerToLeaderBasisOwners, tReceivedFollowerToLeaderBasisOwners );
        moris::communicate_cells( tCommCell, tSendFollowerToLeaderBasisWeights, tReceivedFollowerToLeaderBasisWeights );
        moris::communicate_cells( tCommCell, tSendFollowerToLeaderOffset, tReceivedFollowerToLeaderOffset );

        // clear unused memory
        tSendFollowerToLeaderBasisWeights.clear();
        tSendFollowerToLeaderBasisIds.clear();
        tSendFollowerToLeaderBasisOwners.clear();
        tSendFollowerToLeaderOffset.clear();

        // clear unused capacity
        shrink_to_fit_all( tSendFollowerToLeaderBasisWeights );
        tSendFollowerToLeaderBasisIds.shrink_to_fit();
        tSendFollowerToLeaderBasisOwners.shrink_to_fit();
        tSendFollowerToLeaderOffset.shrink_to_fit();

        /* ---------------------------------------------------------------------------------------- */
        /* Step 4: Use answers to assign IDs to non-owned entities */

        this->handle_requested_shared_follower_basis( aMeshIndex,
                tReceivedFollowerToLeaderBasisIds,
                tReceivedFollowerToLeaderBasisOwners,
                tReceivedFollowerToLeaderBasisWeights,
                tReceivedFollowerToLeaderOffset,
                tBasisIndexToProcs );
    }


    //------------------------------------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_requests_for_shared_follower_basis( moris_index aMeshIndex,
            Cell< Cell< moris_index > >&                                     aBasisIndexToProcs,
            Cell< moris::Cell< moris_id > >&                                 aSendBasisIds )
    {
        // get the bspline mesh info
        Bspline_Mesh_Info*                        tBsplineMeshInfo                  = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< Subphase_Group* > const &    tSPGs                             = tBsplineMeshInfo->mSubphaseGroups;
        Cell< moris::Matrix< IndexMat > > const & tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;

        // get the communication table and map
        uint                                        tCommTableSize = mCommTable.size();
        std::unordered_map< moris_id, moris_index > tProcIdToCommTableIndex;
        tProcIdToCommTableIndex.reserve( tCommTableSize );

        // create  a map from proc id to comm table index
        for ( moris_index iIndex = 0; iIndex < (moris_index)tCommTableSize; ++iIndex )
        {
            tProcIdToCommTableIndex[ mCommTable( iIndex ) ] = iIndex;
        }

        // count number of follower basis functions in the basis data object
        uint tNumFollowerBasis = std::count_if( mBasisData( aMeshIndex ).mFollowerBasis.begin(),    //
                mBasisData( aMeshIndex ).mFollowerBasis.end(),                                      //
                []( moris_index aBasisFlag ) { return aBasisFlag == 1; } );                         // 1 indicates it is  follower basis

        // initialize lists of identifying information
        aBasisIndexToProcs.resize( tCommTableSize );
        aSendBasisIds.resize( tCommTableSize );

        // reserve enough space for all non-owned SPGs  ( over allocation, but better than under allocation )
        std::for_each( aBasisIndexToProcs.begin(), aBasisIndexToProcs.end(), [ tNumFollowerBasis ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumFollowerBasis ); } );
        std::for_each( aSendBasisIds.begin(), aSendBasisIds.end(), [ tNumFollowerBasis ]( Cell< moris_index >& aCell ) { aCell.reserve( tNumFollowerBasis ); } );

        // go through SPGs that executing proc knows about, but doesn't own, ...
        for ( uint iBasisIndex = 0; iBasisIndex < mBasisData( aMeshIndex ).mFollowerBasis.size(); iBasisIndex++ )
        {
            // check if the basis function is not a follower basis function, then skip
            if ( mBasisData( aMeshIndex ).mFollowerBasis( iBasisIndex ) == 0 )
            {
                continue;
            }

            // else condition, if it is find out the SPGs that this basis is interpolating by getting the SPGs that this basis is in
            moris::Matrix< IndexMat > const & tSubphaseGroups = tSubphaseGroupIndsInEnrichedBasis( iBasisIndex );

            for ( moris_index iSPGIndex = 0; iSPGIndex < (moris_index)tSubphaseGroups.numel(); iSPGIndex++ )
            {
                // get the SPG index
                moris_index tSPGIndex = tSubphaseGroups( iSPGIndex );

                // get the SPG
                Subphase_Group* tSPG = tSPGs( tSPGIndex );

                // get the owner of the SPG
                moris_index tSPGOwner = tSPG->get_owner();

                // if the owner is the executing proc, then skip
                if ( tSPGOwner == par_rank() )
                {
                    continue;
                }

                // otherwise add the basis to the list of basis to be requested from that owning proc
                auto tIter = tProcIdToCommTableIndex.find( tSPGOwner );
                MORIS_ASSERT(
                        tIter != tProcIdToCommTableIndex.end(),
                        "Integration_Mesh_Generator::prepare_requests_for_not_owned_subphase_group_IDs() - "
                        "Entity owner (Proc #%i) not found in communication table of current proc #%i when preparing",
                        tSPGOwner,
                        par_rank() );
                moris_index tCommTableIndex = tIter->second;

                //  add the basis index to the list of basis to be requested from that owning proc
                aBasisIndexToProcs( tCommTableIndex ).push_back( iBasisIndex );
            }
        }

        // assemble identifying information for every processor communicated with
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // make the indices unique
            unique( aBasisIndexToProcs( iProc ) );

            // get the number of non-owned SPGs to be sent to each processor processor
            uint tNumSharedBasis = aBasisIndexToProcs( iProc ).size();

            // allocate matrix
            aSendBasisIds( iProc ).resize( tNumSharedBasis );

            // go through the Subphase groups for which IDs will be requested by the other processor
            for ( uint iBasisOrd = 0; iBasisOrd < tNumSharedBasis; iBasisOrd++ )
            {
                // get the index of the subphase group on the executing proc
                moris_index tBasisIndex = aBasisIndexToProcs( iProc )( iBasisOrd );

                // get the basis Id
                moris_index tBasisId = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enr_basis_id_from_enr_basis_index( aMeshIndex, tBasisIndex );

                // store the identifying information of the Subphase group in the output arrays
                aSendBasisIds( iProc )( iBasisOrd ) = tBasisId;
            }
        }

        // size out unused memory
        shrink_to_fit_all( aBasisIndexToProcs );
        shrink_to_fit_all( aSendBasisIds );
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::prepare_answers_for_follower_shared_basis( moris_index const & aMeshIndex,
            moris::Cell< moris::Cell< moris_id > >&                                 aSendFollowerToLeaderBasisIds,
            moris::Cell< moris::Cell< moris_id > >&                                 aSendFollowerToLeaderBasisOwners,
            moris::Cell< moris::Cell< real > >&                                     aSendFollowerToLeaderBasisWeights,
            moris::Cell< moris::Cell< moris_id > >&                                 aSendFollowerToLeaderBasisOffset,
            moris::Cell< moris::Cell< moris_id > > const &                          aReceivedBasisIds )
    {
        // get a pointer or reference to required data for readability
        xtk::Enriched_Interpolation_Mesh* tEnrichedInterpMesh    = mXTKModelPtr->mEnrichedInterpMesh( 0 );
        Basis_Data&                       tBasisData             = mBasisData( aMeshIndex );
        Matrix< IndexMat > const &        tLocalToGlobalBasisMap = tEnrichedInterpMesh->get_enriched_coefficient_local_to_global_map( aMeshIndex );
        Matrix< IdMat > const &           tCommTable             = mXTKModelPtr->mCutIntegrationMesh->get_communication_table();
        uint                              tCommTableSize         = tCommTable.numel();

        // initialize array of answers
        aSendFollowerToLeaderBasisIds.resize( tCommTableSize );
        aSendFollowerToLeaderBasisOwners.resize( tCommTableSize );
        aSendFollowerToLeaderBasisWeights.resize( tCommTableSize );
        aSendFollowerToLeaderBasisOffset.resize( tCommTableSize );

        // check that the received data is complete
        MORIS_ASSERT(
                aReceivedBasisIds.size() == tCommTableSize,
                "Basis_Processor::prepare_answers_for_owned_subphase_groups() - Received information incomplete." );

        // go through the list of processors in the array of ID requests
        for ( uint iProc = 0; iProc < tCommTableSize; iProc++ )
        {
            // get the number of entity IDs requested from the current proc position
            uint tNumReceivedReqs = aReceivedBasisIds( iProc ).size();

            // size the list of answers / IDs accordingly
            aSendFollowerToLeaderBasisIds( iProc ).reserve( tNumReceivedReqs );
            aSendFollowerToLeaderBasisOwners( iProc ).reserve( tNumReceivedReqs * 4 );
            aSendFollowerToLeaderBasisWeights( iProc ).reserve( tNumReceivedReqs * 4 );
            aSendFollowerToLeaderBasisOffset( iProc ).reserve( tNumReceivedReqs * 4 );

            aSendFollowerToLeaderBasisOffset( iProc ).push_back( 0 );

            uint tNumBasis = 0;
            // iterate through SPs for which the IDs are requested
            for ( uint iBasisOrd = 0; iBasisOrd < tNumReceivedReqs; iBasisOrd++ )
            {
                // get the ID and index of the received SP
                moris_id tBasisId    = aReceivedBasisIds( iProc )( iBasisOrd );
                moris_id tBasisIndex = tEnrichedInterpMesh->get_enr_basis_index_from_enr_basis_id( aMeshIndex, tBasisId );

                // insert at the end of the vector of IDs
                aSendFollowerToLeaderBasisIds( iProc )
                        .insert( aSendFollowerToLeaderBasisIds( iProc ).size(),    //
                                tBasisData.mFollowerToLeaderBasis( tBasisIndex ).begin(),
                                tBasisData.mFollowerToLeaderBasis( tBasisIndex ).end() );

                aSendFollowerToLeaderBasisOwners( iProc )
                        .insert( aSendFollowerToLeaderBasisOwners( iProc ).size(),    //
                                tBasisData.mFollowerToLeaderBasisOwners( tBasisIndex ).begin(),
                                tBasisData.mFollowerToLeaderBasisOwners( tBasisIndex ).end() );

                aSendFollowerToLeaderBasisWeights( iProc )
                        .insert( aSendFollowerToLeaderBasisWeights( iProc ).size(),    //
                                tBasisData.mFollowerToLeaderBasisWeights( tBasisIndex ).begin(),
                                tBasisData.mFollowerToLeaderBasisWeights( tBasisIndex ).end() );

                tNumBasis += tBasisData.mFollowerToLeaderBasis( tBasisIndex ).size();

                aSendFollowerToLeaderBasisOffset( iProc ).push_back( tNumBasis );

            }    // end for: communication for each entity with current processor

            // transform the indices to ids using the map
            std::transform( aSendFollowerToLeaderBasisIds( iProc ).begin(),    //
                    aSendFollowerToLeaderBasisIds( iProc ).end(),              //
                    aSendFollowerToLeaderBasisIds( iProc ).begin(),            //
                    [ &tLocalToGlobalBasisMap ]( moris_index const & aIndex ) { return tLocalToGlobalBasisMap( aIndex ); } );
        }                                                                      // end for: communication list for each processor
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::handle_requested_shared_follower_basis(
            moris_index                                    aMeshIndex,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedFollowerToLeaderBasisIds,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedFollowerToLeaderBasisOwners,
            moris::Cell< moris::Cell< real > > const &     aReceivedFollowerToLeaderBasisWeights,
            moris::Cell< moris::Cell< moris_id > > const & aReceivedFollowerToLeaderBasisOffset,
            moris::Cell< moris::Cell< moris_id > > const & tBasisIndexToProcs )
    {
        // get a pointer or reference to required data for readability
        Basis_Data&                  tBasisData     = mBasisData( aMeshIndex );
        Enriched_Interpolation_Mesh& tEnrInterpMesh = mXTKModelPtr->get_enriched_interp_mesh( aMeshIndex );

        // process answers from each proc communicated with
        for ( uint iProc = 0; iProc < aReceivedFollowerToLeaderBasisIds.size(); iProc++ )
        {
            // check the number of requests and answers from the current proc
            MORIS_ASSERT( aReceivedFollowerToLeaderBasisOffset( iProc ).size() - 1 == tBasisIndexToProcs( iProc ).size(),
                    "Integration_Mesh_Generator::handle_requested_subphase_group_ID_answers() - "
                    "Number of SPG ID requests and answers are not the same." );

            // assign IDs to each communicated entity
            for ( uint iBasisOrd = 0; iBasisOrd < tBasisIndexToProcs( iProc ).size(); iBasisOrd++ )
            {
                // get the current SPG index and ID from the data provided
                moris_index tBasisIndex = tBasisIndexToProcs( iProc )( iBasisOrd );

                // determine the range of basis IDs to be added to the processor domain
                uint tStartRange = aReceivedFollowerToLeaderBasisOffset( iProc )( iBasisOrd );
                uint tEndRange   = aReceivedFollowerToLeaderBasisOffset( iProc )( iBasisOrd + 1 );

                for ( uint iBase = tStartRange; iBase < tEndRange; iBase++ )
                {
                    // get the basis ID that needs to be added to the processor domain
                    moris_index tBasisId = aReceivedFollowerToLeaderBasisIds( iProc )( iBase );

                    // If the root basis does not exist in the current processor add it to the list and and generate a new index for it
                    // add this basis to the mesh if it does not exists on the current partition
                    if ( !tEnrInterpMesh.basis_exists_on_partition( aMeshIndex, tBasisId ) )
                    {
                        // NOTE: the bulk phase information will never be used later so it is set to zero
                        moris_index tBulkPhase = 0;

                        // add basis ID to partition
                        tEnrInterpMesh.add_basis_function( aMeshIndex, tBasisId, aReceivedFollowerToLeaderBasisOwners( iProc )( iBase ), tBulkPhase );
                    }

                    // find and store the basis index local to the executing processor
                    moris_index tLeaderBasisIndex = tEnrInterpMesh.get_enr_basis_index_from_enr_basis_id( aMeshIndex, tBasisId );

                    tBasisData.mFollowerToLeaderBasis( tBasisIndex ).push_back( tLeaderBasisIndex );
                    tBasisData.mFollowerToLeaderBasisOwners( tBasisIndex ).push_back( aReceivedFollowerToLeaderBasisOwners( iProc )( iBase ) );
                    tBasisData.mFollowerToLeaderBasisWeights( tBasisIndex ).push_back( aReceivedFollowerToLeaderBasisWeights( iProc )( iBase ) );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------

    void
    Basis_Processor::update_comm_table()
    {
        // update the communication table for the enriched interpolation mesh
        mXTKModelPtr->mEnrichedInterpMesh( 0 )->update_communication_table( mCommTable );
    }
}    // namespace xtk
