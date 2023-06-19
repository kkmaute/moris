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
    }

    // ----------------------------------------------------------------------------------

    Basis_Processor::~Basis_Processor()
    {
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
        for ( const auto& iMeshIndex : mMeshIndices )
        {
            // creat a cell to keep track of the basis that has been grouped
            moris::Cell< uint > tBasisHasBeenUsed( mXTKModelPtr->mEnrichment->mEnrichmentData( iMeshIndex ).mEnrLvlOfEnrBf.size(), 0 );

            moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup = mXTKModelPtr->mEnrichment->mEnrichmentData( iMeshIndex ).mEnrichedBasisInSubphaseGroup;

            Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( iMeshIndex );
            // Get the number of subphase groups (on the current proc) and reisze
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

                // If volume is less than a specific treshhold then add it to the list of
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
        mRootSPGIndex.resize( tBsplineMeshInfo->get_num_SPGs(), -1 );

        // loop over the B-spline elements to determine the cut elements
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            moris::Cell< moris_index > const & tSPGIndicesInBsplineCell = tBsplineMeshInfo->get_SPG_indices_in_bspline_cell( iBspElem );

            // B-spline element is not cut thus does not need to be grouped
            if ( 1 == tSPGIndicesInBsplineCell.size() )
            {
                // assign as the root as itself because it is not cut
                mRootSPGIndex( tSPGIndicesInBsplineCell( 0 ) ) = tSPGIndicesInBsplineCell( 0 );

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
            if ( iSPGIndex != mRootSPGIndex( iSPGIndex ) )
            {
                // get the enriched basis indices present in the SPGs
                moris::Cell< moris_index > const & tEnrichedBasis = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                // loop over the enriched basis and if there a is a follower then assign the root as the owner, might be overwritten later
                for ( const auto& iEnrichedBase : tEnrichedBasis )
                {
                    // if it a bad basis then this is true
                    if ( mBasisData( aMeshIndex ).mFollowerBasis( iEnrichedBase ) == 1 )
                    {
                        mBasisData( aMeshIndex ).mFollowerBasisOwningCell( iEnrichedBase ) = mRootSPGIndex( iSPGIndex );
                    }
                }
            }
        }
    }

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
        moris_index tFieldIndex = mXTKModelPtr->mEnrichedIntegMesh( 0 )->create_field( tCellFields( 0 ), EntityRank::ELEMENT, MORIS_INDEX_MAX );

        moris::Matrix< moris::DDRMat > tCellIdField( 1, tNumIGCells, -1.0 );

        // loop over the SPGs and assign dpg id to each cell index
        for ( uint iSPG = 0; iSPG < tNumSubphaseGroups; iSPG++ )
        {
            const moris::Cell< moris_index >& tIGCellIndices = tBsplineMeshInfo->mSubphaseGroups( iSPG )->get_ig_cell_indices_in_group();

            std::for_each( tIGCellIndices.begin(), tIGCellIndices.end(), [ &tCellIdField, this, iSPG ]( moris_index aIGCellIndex )    //
                    { tCellIdField( aIGCellIndex ) = mRootSPGIndex( iSPG ); } );
        }

        // add the field data to the mesh
        mXTKModelPtr->mEnrichedIntegMesh( 0 )->add_field_data( tFieldIndex, EntityRank::ELEMENT, tCellIdField, MORIS_INDEX_MAX );

        // write this field to the exodus file if there is a cell that is not grouped
        auto it = std::find( mRootSPGIndex.begin(), mRootSPGIndex.end(), -1 );
        if ( it != mRootSPGIndex.end() )
        {
            // output the xtk mesh for debugging
            mXTKModelPtr->mEnrichedIntegMesh( 0 )->write_mesh( mParameterList );

            MORIS_ASSERT( false, "The SPG with the index %d is not grouped, refer to XTK IG mesh for more info", *it );
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::perform_basis_extention()
    {
        Tracer tTracer( "XTK", "Basis Processor", "Basis Extension" );

        // generate the subsphase to enriched basis map
        mXTKModelPtr->mEnrichment->construct_enriched_basis_in_subphase_group_map();

        for ( const auto& iMeshIndex : mMeshIndices )
        {
            this->construct_follower_basis_using_volume( iMeshIndex );

            this->compute_averaging_weights( iMeshIndex );

            this->construct_cell_association( iMeshIndex );

            this->visualize_cell_aggregates( iMeshIndex );

            this->construct_follower_to_leader_basis_weights_indices( iMeshIndex );

            this->replace_t_matrices( iMeshIndex );
        }
    }

    void
    Basis_Processor::perform_cell_agglomeration()
    {
        Tracer tTracer( "XTK", "Basis Processor", "Cell Agglomeration" );

        // generate the subsphase to enriched basis map
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


    void
    Basis_Processor::compute_averaging_weights( moris_index aMeshIndex )
    {
        // get a reference to pointer to the required data
        Cell< moris::Matrix< IndexMat > > const & tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;
        Bspline_Mesh_Info*                        tBsplineMeshInfo                  = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );
        moris::Cell< Subphase_Group* >&           tSPGs                             = tBsplineMeshInfo->mSubphaseGroups;
        uint                                      tNumEnrichedBF                    = tSubphaseGroupIndsInEnrichedBasis.size();


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
    }

    void
    Basis_Processor::construct_follower_basis_using_volume( moris_index aMeshIndex )
    {
        // compute the threshold volume to decide good/bad basis
        real tVolThreshold = this->compute_threshold_volume( aMeshIndex );

        // get the spg to enr basis and its transpose( enr basis to spg map)
        moris::Cell< moris::Cell< moris_index > > const & tEnrichedBasisInSubphaseGroup     = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisInSubphaseGroup;
        Cell< moris::Matrix< IndexMat > > const &         tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;

        // get the bspline mesh info to access the spgs
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the spgs from b-spline mesh info
        moris::Cell< Subphase_Group* > tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the number of active B-spline elements ( this is un-enriched)
        uint tNumBspElems      = tBsplineMeshInfo->get_num_Bspline_cells();
        uint tNumEnrichedBasis = tSubphaseGroupIndsInEnrichedBasis.size();
        // Get the number of subphase groups (on the current proc) and reisze
        uint tNumSubphaseGroups = tBsplineMeshInfo->get_num_SPGs();

        // allocate size for the basis info
        // for follower basis assume all are follower and then change the leader one to zero
        mBasisData( aMeshIndex ).mFollowerBasis.resize( tNumEnrichedBasis, 1 );
        mBasisData( aMeshIndex ).mFollowerToLeaderBasis.resize( tNumEnrichedBasis );
        mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights.resize( tNumEnrichedBasis );

        // allocate root spgs index for the spgs, initialize with -1 to catch the errors later
        mRootSPGIndex.resize( tNumSubphaseGroups, -1 );

        // loop over the b-spline elements to get the SPGs living on them
        for ( uint iBspElem = 0; iBspElem < tNumBspElems; iBspElem++ )
        {
            // get the SPGs that live on the bspline element
            moris::Cell< moris_index > const & tSPGIndicesInBsplineCell = tBsplineMeshInfo->get_SPG_indices_in_bspline_cell( iBspElem );

            // if only one spg then B-spline element is not cut
            if ( 1 == tSPGIndicesInBsplineCell.size() )
            {
                // assign as the root as itself because it is not cut
                mRootSPGIndex( tSPGIndicesInBsplineCell( 0 ) ) = tSPGIndicesInBsplineCell( 0 );

                // access the associated basis of the SPG, all marked as not follower(i.e. leader)
                // get the basis that are attached to this SPG
                moris::Cell< moris_index > const & tEnrichedBFsInSPGs = tEnrichedBasisInSubphaseGroup( tSPGIndicesInBsplineCell( 0 ) );

                // loop over the basis and mark them as leader
                for ( const auto& iEnrichedBF : tEnrichedBFsInSPGs )
                {
                    // set it to false that basis is not a follower
                    mBasisData( aMeshIndex ).mFollowerBasis( iEnrichedBF ) = 0;
                }

                // go to the next b-spline element
                continue;
            }

            // if the bspline element is cut find if the volume is below a certain treshhold;
            for ( const auto& iSPGIndex : tSPGIndicesInBsplineCell )
            {
                // get the IG cells and integarte over them
                // get the subphase group based on the index
                Subphase_Group* tSubphaseGroup = tBsplineMeshInfo->mSubphaseGroups( iSPGIndex );

                // otherwise get the IG cells in that SPG and integrate over the volume of those
                const moris::Cell< moris_index >& tIGCellsInGroup = tSubphaseGroup->get_ig_cell_indices_in_group();

                // compute the volume of the cell
                real tVolume = std::accumulate( tIGCellsInGroup.begin(), tIGCellsInGroup.end(), 0.0, [ & ]( real aVol, int index ) { 
                                    const mtk::Cell & tCell = mXTKModelPtr->mCutIntegrationMesh->get_mtk_cell(index) ;
                                    return aVol + tCell.compute_cell_measure(); } );

                // If volume is less than a specific treshhold then add it to the list of
                if ( tVolume >= tVolThreshold )
                {
                    // if a volume threshold is met mark the spgs as its own root
                    mRootSPGIndex( iSPGIndex ) = iSPGIndex;

                    // get the enriched basis that are interpolating into the spg
                    moris::Cell< moris_index > const & tEnrichedBFsInSPGs = tEnrichedBasisInSubphaseGroup( iSPGIndex );

                    // mark these basis as not follower basis
                    for ( const auto& iEnrichedBF : tEnrichedBFsInSPGs )
                    {
                        // set it to false that basis is not a follower
                        mBasisData( aMeshIndex ).mFollowerBasis( iEnrichedBF ) = 0;
                    }
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
        // get the bspline mesh info to access the spgs
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the spgs from b-spline mesh info
        moris::Cell< Subphase_Group* >& tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the subphase connectivity
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tSubphaseGroupNeighborhood = mXTKModelPtr->mCutIntegrationMesh->get_subphase_group_neighborhood( aMeshIndex );

        // initialize cells to keep possible candidate root cells and their distance
        moris::Cell< moris_index > tCandidateRootCell;
        moris::Cell< moris_index > tDistanceCell;

        // loop over the spgs and process the ones that do not have root spgs
        for ( size_t iSPGIndex = 0; iSPGIndex < mRootSPGIndex.size(); iSPGIndex++ )
        {
            // this mean the SPG does not have not have a root cell
            tCandidateRootCell.resize( 0 );
            tDistanceCell.resize( 0 );

            // this means that the root spg is not assigned
            if ( -1 == mRootSPGIndex( iSPGIndex ) )
            {
                // get the first degree neighbour of the SPG
                // TODO: only first degree neighbour is considered, need to consider higher degree neighbours
                std::shared_ptr< moris::Cell< moris_index > > tNeighbour = tSubphaseGroupNeighborhood->mSubphaseToSubPhase( iSPGIndex );

                // between the neighbour get the ones that are is a root cell
                std::copy_if( tNeighbour->begin(), tNeighbour->end(), std::back_inserter( tCandidateRootCell ), [ & ]( moris_index aNeighbourSPGIndex )    //
                        {
                            return mRootSPGIndex( aNeighbourSPGIndex ) == aNeighbourSPGIndex;
                        } );

                // assign as root cell the cell with smaller index
                if ( tCandidateRootCell.size() )
                {
                    // define a lambda function with sorting criteria
                    auto tSortBasedOntheBsplineCellIndex = [ & ]( moris_index aSPGIndex1, moris_index aSPGIndex2 ) {
                        return tSPGs( aSPGIndex1 )->get_bspline_cell_index() < tSPGs( aSPGIndex2 )->get_bspline_cell_index();
                    };

                    // this sorting is necessary such that in the event of a tie in euclaidam distance the cell with smallest B-spline index
                    std::sort( tCandidateRootCell.begin(), tCandidateRootCell.end(), tSortBasedOntheBsplineCellIndex );

                    // loop over the candidate cells that have the find the closest
                    // TODO: could be improved by an elimination strategy instead of computing all distances
                    for ( const auto& iNeighbourSPGIndex : tCandidateRootCell )
                    {
                        // compute euclidean distance of the two cells
                        uint         tLevel = 0;
                        const luint* tIJK   = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tSPGs( iSPGIndex )->get_bspline_cell_index(), tLevel );

                        // find the root spg of the neighbour SPG and its ijk
                        moris_index  tRootOfNeighbourSPG = mRootSPGIndex( iNeighbourSPGIndex );
                        uint         tLevelRoot          = 0;
                        const luint* tIJKRoot            = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tSPGs( tRootOfNeighbourSPG )->get_bspline_cell_index(), tLevelRoot );

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
                    auto tSelectedNeigbourSPGIndex = std::min_element( tDistanceCell.begin(), tDistanceCell.end() );
                    auto tIndex                    = std::distance( tDistanceCell.begin(), tSelectedNeigbourSPGIndex );
                    mRootSPGIndex( iSPGIndex )     = mRootSPGIndex( tCandidateRootCell( tIndex ) );
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
        moris::Cell< Subphase_Group* >&           tSPGs                             = tBsplineMeshInfo->mSubphaseGroups;

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
                    moris_index tRootSPGIndex = mRootSPGIndex( iSPGIndex );

                    // get the bspline cell of the root cell
                    moris_index tRootBSplineCellIndex = tSPGs( tRootSPGIndex )->get_bspline_cell_index();

                    // get the bspline cell of the root cell
                    moris_index tExtentionBSplineCellIndex = tSPGs( iSPGIndex )->get_bspline_cell_index();

                    // initialize the basis extension/projection data
                    moris::Cell< moris::Cell< mtk::Vertex* > > tRootBsplineBasis;
                    moris::Cell< mtk::Vertex* >                tExtendedBsplineBasis;
                    moris::Cell< Matrix< DDRMat > >            tWeights;

                    // get the L2-projection matrix along from the root to the extended b-spline
                    mXTKModelPtr->mBackgroundMesh->get_L2_projection_matrix( aMeshIndex, tRootBSplineCellIndex, tExtentionBSplineCellIndex, tRootBsplineBasis, tExtendedBsplineBasis, tWeights );

                    // This two maps help to get the the enriched basis index from the
                    moris::Cell< moris_index > const &             tBGBasisIndicesRoot = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisIndices( tRootSPGIndex );
                    moris::Cell< moris_index > const &             tBGBasisLevelsRoot  = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupBGBasisEnrLev( tRootSPGIndex );
                    std::unordered_map< moris_index, moris_index > tBasisToLocalIndexMapRoot;

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


                    //! comment
                    moris::Cell< moris_index > tRootBGBSplineIndices( tRootBsplineBasis( tNonEnrBasisOrd ).size() );

                    std::transform( tRootBsplineBasis( tNonEnrBasisOrd ).begin(), tRootBsplineBasis( tNonEnrBasisOrd ).end(), tRootBGBSplineIndices.begin(), []( mtk::Vertex* aVertex ) { return aVertex->get_index(); } );

                    // find the max size in enriched basis index
                    auto tCellMaxsize = std::max_element( tRootBsplineBasis.begin(), tRootBsplineBasis.end(),          //
                            []( moris::Cell< mtk::Vertex* >& aCellFirst, moris::Cell< mtk::Vertex* >& aCellSecond )    //
                            { return aCellFirst.size() < aCellSecond.size(); } );

                    // reserve enough space for the enriched basis indices
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasis( iBasisIndex ).reserve( tCellMaxsize->size() );
                    mBasisData( aMeshIndex ).mFollowerToLeaderBasisWeights( iBasisIndex ).reserve( tCellMaxsize->size() );

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
                    }
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Basis_Processor::replace_t_matrices( moris_index aMeshIndex )
    {
        // get a reference or a pointer to reuqired data
        Cell< Interpolation_Cell_Unzipped* > const & tEnrichedIPCells                  = mXTKModelPtr->mEnrichedInterpMesh( 0 )->get_enriched_interpolation_cells();
        Matrix< IndexMat > const &                   tBasisEnrichmentIndices           = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mEnrichedBasisIndexToId;
        Cell< moris::Matrix< IndexMat > > const &    tSubphaseGroupIndsInEnrichedBasis = mXTKModelPtr->mEnrichment->mEnrichmentData( aMeshIndex ).mSubphaseGroupIndsInEnrichedBasis;

        // define an unordered set to keep track of the vertices that have been already processed
        std::unordered_set< moris_index > tListOfBGVertices;
        tListOfBGVertices.reserve( tSubphaseGroupIndsInEnrichedBasis.size() );

        // find the max size that a basis will be replaced by
        auto tCellMaxsize = std::max_element( mBasisData( aMeshIndex ).mFollowerToLeaderBasis.begin(), mBasisData( aMeshIndex ).mFollowerToLeaderBasis.end(),    //
                []( moris::Cell< moris_index >& aCellFirst, moris::Cell< moris_index >& aCellSecond )                                                            //
                { return aCellFirst.size() < aCellSecond.size(); } );

        // loop over the subphase groups to find the problematic ones
        for ( size_t iSPGIndex = 0; iSPGIndex < mRootSPGIndex.size(); iSPGIndex++ )
        {
            // if it is an extended cell
            if ( (moris_index)iSPGIndex != mRootSPGIndex( iSPGIndex ) )
            {
                // get the all the lagrange cells of the that SPG index
                moris::Cell< moris_index > const & tUIPCIndices = mXTKModelPtr->mEnrichment->get_UIPC_indices_on_SPG( aMeshIndex, iSPGIndex );

                // loop over the unzipped lagrange cells and evaluate t-matrices
                for ( const auto& iUIPCIndex : tUIPCIndices )
                {
                    // get base cell of the lagrange mesh
                    const Interpolation_Cell_Unzipped* tConsttUIPC = tEnrichedIPCells( iUIPCIndex );

                    // Get the lagrange vertices of the UIPC to obtain the vertex interpolation object of each lagrange node
                    moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tXTKUIPVs = tConsttUIPC->get_xtk_interpolation_vertices();

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
                        moris::Matrix< moris::IndexMat > const &                      tBasisIndices = tVertexEnrichment->get_basis_indices();
                        moris::Matrix< moris::IndexMat > const &                      tBasisIds     = tVertexEnrichment->get_basis_ids();
                        moris::Matrix< moris::DDRMat >&                               tBasisWeights = tVertexEnrichment->get_basis_weights();
                        std::unordered_map< moris::moris_index, moris::moris_index >& tBasisMap     = tVertexEnrichment->get_basis_map();
                        tBasisMap.clear();

                        // probably need to do a resize call here
                        // overallocation basis Indices , assume that every basis index need to be replaced
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisIndices( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        moris::Matrix< moris::IndexMat > tAgglomeratedBasisIds( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        moris::Matrix< DDRMat >          tAgglomeratedBasisWeights( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ), 1 );
                        tBasisMap.reserve( tCellMaxsize->size() * ( tBasisIndices.numel() + 1 ) );

                        // initialize a counter to count how many basis will be added and replaced
                        uint tBasisCounter = 0;

                        // loop over the basis indices and replace the ones that are not in the follower basis
                        for ( uint iBC = 0; iBC < tBasisIndices.numel(); iBC++ )
                        {
                            // get a reference to the basis index and weight
                            moris::moris_index const & tBasisIndex  = tBasisIndices( iBC );
                            real&                      tBasisWeight = tBasisWeights( iBC );

                            // if it is a good basis keep it
                            if ( mBasisData( aMeshIndex ).mFollowerBasis( tBasisIndex ) == 0 )
                            {
                                // if it iis not in the map add it to map
                                if ( tBasisMap.find( tBasisIndex ) == tBasisMap.end() )
                                {
                                    tAgglomeratedBasisIndices( tBasisCounter ) = tBasisIndex;
                                    tAgglomeratedBasisIds( tBasisCounter )     = tBasisIds( iBC );
                                    tAgglomeratedBasisWeights( tBasisCounter ) = tBasisWeights( iBC );
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

                                // loop over and replace the t matrix
                                for ( uint iSlaveBasisOrd = 0; iSlaveBasisOrd < tFollowerBasis.size(); iSlaveBasisOrd++ )
                                {
                                    moris_index const & tFollowerBasisIndex = tFollowerBasis( iSlaveBasisOrd );

                                    //  if not found in the map
                                    if ( tBasisMap.find( tFollowerBasisIndex ) == tBasisMap.end() )
                                    {
                                        // assign the basis index and weight
                                        tAgglomeratedBasisIndices( tBasisCounter ) = tFollowerBasisIndex;
                                        tAgglomeratedBasisIds( tBasisCounter )     = tBasisEnrichmentIndices( tFollowerBasisIndex );
                                        tAgglomeratedBasisWeights( tBasisCounter ) = tBasisWeight * tFollowerBasisWeights( iSlaveBasisOrd );
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

                                // This is done for the puprose of writing the basis functions
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

                        // add basis information to the vertex enrichment
                        tVertexEnrichment->add_basis_information( tAgglomeratedBasisIndices, tAgglomeratedBasisIds );
                        tVertexEnrichment->add_basis_weights( tAgglomeratedBasisIndices, tAgglomeratedBasisWeights );
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
        Bspline_Mesh_Info* tBsplineMeshInfo = mXTKModelPtr->mEnrichment->mBsplineMeshInfos( aMeshIndex );

        // get the subphase groups
        moris::Cell< Subphase_Group* > tSPGs = tBsplineMeshInfo->mSubphaseGroups;

        // get the number of active B-spline elements and subphase groups
        uint tNumSubphaseGroups = tBsplineMeshInfo->get_num_SPGs();

        // initialize a counter to count number of spg that have been touched
        size_t tNumTouchedSPGs = std::count_if( mRootSPGIndex.begin(), mRootSPGIndex.end(), []( moris_index const & aRootSPGIndex ) {
            return aRootSPGIndex >= 0;
        } );


        // initialize candidate root cell and distance cell( candidate root cell is the list of potential root cells for the current SPG)
        // Distance cell is the list of distances from the current SPG to the candidate root cells
        moris::Cell< moris_index > tCandidateRootCell;
        moris::Cell< moris_index > tDistanceCell;

        // until all SPGs are aggregated run the graph algorithm to generate aggregates
        while ( tNumTouchedSPGs != tNumSubphaseGroups )
        {
            // loop over the cut bspline cells ( corrsopsing spgs)
            for ( size_t iSPGIndex = 0; iSPGIndex < (size_t)tNumSubphaseGroups; iSPGIndex++ )
            {
                tCandidateRootCell.resize( 0 );
                tDistanceCell.resize( 0 );

                // this means that it is not assigned
                if ( -1 == mRootSPGIndex( iSPGIndex ) )
                {
                    std::shared_ptr< moris::Cell< moris_index > > tNeighbour = tSubphaseGroupNeighborhood->mSubphaseToSubPhase( iSPGIndex );

                    // get the cells that are already touched
                    std::copy_if( tNeighbour->begin(), tNeighbour->end(), std::back_inserter( tCandidateRootCell ), [ & ]( moris_index aNiegbourSPGIndex )    //
                            {
                                return mRootSPGIndex( aNiegbourSPGIndex ) != -1;
                            } );

                    // assign as root cell the cell with smaller index
                    if ( tCandidateRootCell.size() )
                    {
                        auto tSortBasedOntheBsplineCellIndex = [ & ]( moris_index aSPGIndex1, moris_index aSPGIndex2 ) {
                            return tSPGs( aSPGIndex1 )->get_bspline_cell_index() < tSPGs( aSPGIndex2 )->get_bspline_cell_index();
                        };

                        // this sorting is necessary such that in the event of a tie in euclidean distance the cell with smallest B-spline index
                        std::sort( tCandidateRootCell.begin(), tCandidateRootCell.end(), tSortBasedOntheBsplineCellIndex );

                        // find the euclidean distance between all candicate points
                        // TODO: could be improved by an elimination strategy instead of computing all distances
                        for ( const auto& iNeighbourSPGIndex : tCandidateRootCell )
                        {
                            // compute euclidean distance of the two cells
                            uint         tLevel = 0;
                            const luint* tIJK   = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tSPGs( iSPGIndex )->get_bspline_cell_index(), tLevel );

                            // find the root spg of the neighbour SPG and its ijk
                            moris_index  tRootOfNeighbourSPG = mRootSPGIndex( iNeighbourSPGIndex );
                            uint         tLevelRoot          = 0;
                            const luint* tIJKRoot            = mXTKModelPtr->mBackgroundMesh->get_bspline_element_ijk_level( aMeshIndex, tSPGs( tRootOfNeighbourSPG )->get_bspline_cell_index(), tLevelRoot );

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
                        auto tSelectedNeigbourSPGIndex = std::min_element( tDistanceCell.begin(), tDistanceCell.end() );
                        auto tIndex                    = std::distance( tDistanceCell.begin(), tSelectedNeigbourSPGIndex );
                        mRootSPGIndex( iSPGIndex )     = mRootSPGIndex( tCandidateRootCell( tIndex ) );
                        tNumTouchedSPGs++;
                    }
                }
            }
        }
    }

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
            if ( mRootSPGIndex( iSPGIndex ) != iSPGIndex )
            {
                // get the bspline cell of the root cell
                moris_index tRootBSplineCellIndex = tSPGs( mRootSPGIndex( iSPGIndex ) )->get_bspline_cell_index();

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

                    // get a const pointer to UPIC such that the following xtk specific function can be called
                    const Interpolation_Cell_Unzipped* tConsttUIPC = tUIPC;

                    // Get the lagrange vertices of the UIPC to obtain the vertex interpolation object of each lagrange node
                    moris::Cell< xtk::Interpolation_Vertex_Unzipped* > const & tXTKUIPVs = tConsttUIPC->get_xtk_interpolation_vertices();

                    // loop over the unzipped interpolation vertices and get their t-matrix info to overwrite
                    for ( size_t iLocalVertIndex = 0; iLocalVertIndex < tXTKUIPVs.size(); iLocalVertIndex++ )
                    {
                        // get the vertex enrichment object
                        xtk::Vertex_Enrichment* tVertexEnrichment = tXTKUIPVs( iLocalVertIndex )->get_xtk_interpolation( aMeshIndex );

                        // access the basis indices
                        moris::Matrix< moris::IndexMat > const &                      tBasisIndices = tVertexEnrichment->get_basis_indices();
                        moris::Matrix< moris::DDRMat >&                               tBasisWeights = tVertexEnrichment->get_basis_weights();
                        std::unordered_map< moris::moris_index, moris::moris_index >& tBasisMap     = tVertexEnrichment->get_basis_map();

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
                                moris_index tRootSPGIndex = mRootSPGIndex( iSPGIndex );

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
                                std::unordered_map< moris_index, moris_index > tBasisToLocalIndexMap;
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
                        tBasisMap.reserve( iBasisCounter );

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
                        tVertexEnrichment->add_basis_weights( tAgglomeratedBasisIndices, tAgglomeratedBasisWeights );
                    }
                }
            }
        }
    }
}    // namespace xtk
