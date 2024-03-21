/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_MTK_Contact_Mesh_Editor.cpp
 *
 */

#include "cl_MTK_Contact_Mesh_Editor.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_MappingResult.hpp"
#include "cl_MTK_PointPairs.hpp"
#include "cl_MTK_Nonconformal_Side_Cluster.hpp"
#include "cl_Matrix_Arma_Dynamic.hpp"
#include "cl_Vector.hpp"
#include "fn_assert.hpp"
#include "fn_join_horiz.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include <algorithm>
#include <iterator>
#include <set>
#include <string>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <map>

namespace moris::mtk
{
    std::unordered_map< moris_index, std::map< std::pair< moris_index, moris_index >, Contact_Mesh_Editor::ResultIndices > >
    Contact_Mesh_Editor::extract_cluster_pairing( MappingResult const &aMappingResult ) const
    {
        // Map from every cluster on the leader side to every cluster on which a mapping was successful on the follower side.
        // The second map stores all indices in the aMappingResult that belong to this cluster pair.
        // This stage does not yet group the indices by the cell pairs.
        // Note: Even if a cluster has not been succesfully mapped, we need to store an empty map for this cluster!
        // This is necessary to add "no-op" nonconformal side sets to the IGMesh to get the displacement values at those nodes as well.
        uint const tNumClusters = mSideSets( aMappingResult.mSourceMeshIndex )->get_num_clusters_on_set();

        std::unordered_map< moris_index, std::map< std::pair< moris_index, moris_index >, Vector< moris_index > > > tClusterPairs( tNumClusters );
        for ( moris_index i = 0; i < static_cast< moris_index >( aMappingResult.mSourceClusterIndex.size() ); ++i )
        {
            moris_index const &tSourceClusterIndex = aMappingResult.mSourceClusterIndex( i );
            if ( aMappingResult.mTargetClusterIndex( i ) == -1 )    // the mapping on this cluster (in this result index) was not successful
            {
                if ( tClusterPairs.find( tSourceClusterIndex ) == tClusterPairs.end() )    // the cluster has not been mapped yet (prevents overwriting of the map)
                {
                    tClusterPairs[ tSourceClusterIndex ] = {};    // empty map for clusters that have not been mapped
                }
            }
            else
            {
                moris_index const &tTargetClusterIndex = aMappingResult.mTargetClusterIndex( i );
                moris_index const &tTargetSideSetIndex = aMappingResult.mTargetSideSetIndices( i );

                // to locate the correct cluster pair, we need to provide the cluster index as well as the side set on which the cluster is located
                auto const tTargetLocator = std::make_pair( tTargetClusterIndex, tTargetSideSetIndex );
                tClusterPairs[ tSourceClusterIndex ][ tTargetLocator ].push_back( i );
            }
        }
        return tClusterPairs;
    }

    std::map< Contact_Mesh_Editor::CellPair, Contact_Mesh_Editor::ResultIndices >
    Contact_Mesh_Editor::extract_cell_pairing( MappingResult const &aMappingResult, Vector< moris_index > aResultIndices )
    {
        std::map< Contact_Mesh_Editor::CellPair, Contact_Mesh_Editor::ResultIndices > tCellPairs;
        for ( auto const &tResultIndex : aResultIndices )
        {
            moris_index const &tSourceCellIndex = aMappingResult.mSourceCellIndex( tResultIndex );
            moris_index const &tTargetCellIndex = aMappingResult.mTargetCellIndices( tResultIndex );
            auto const         tCellPair        = std::make_pair( tSourceCellIndex, tTargetCellIndex );
            tCellPairs[ tCellPair ].push_back( tResultIndex );
        }
        return tCellPairs;
    }

    std::map< Contact_Mesh_Editor::SetPair, Vector< Nonconformal_Side_Cluster > >
    Contact_Mesh_Editor::convert_mapping_result_to_nonconformal_side_clusters( MappingResult const &aMappingResult ) const
    {
        moris_index const tSourceMeshIndex = aMappingResult.mSourceMeshIndex;

        // extract the cluster pairs and the corresponding cell pairs from the mapping result this will act as a basis to create the nonconformal side clusters
        auto tClusterPairs = extract_cluster_pairing( aMappingResult );
        // Since we know the number of unique cluster-pairs, we can reserve the correct amount of memory for the nonconformal side clusters
        std::map< SetPair, Vector< Nonconformal_Side_Cluster > > tNonconformalSideClusters;

        // loop over every cluster of the source side and create a nonconformal side cluster for each mapped target cluster
        for ( const auto &[ tSourceClusterIndex, tTargetClusterToResultMap ] : tClusterPairs )
        {
            // every cluster contains the list of points that got mapped from one cell to another. The first integration point pairs object corresponds to the first
            // cell in the cluster and so on. If the did not get mapped, the lists will be empty. However, we still need to create at least one empty cluster object
            // for each cluster to be able to request displacement values at the nodes of those clusters during the remapping.
            if ( tTargetClusterToResultMap.empty() )
            {
                // Simply use an arbitrary cluster from another set... since this cluster will never be actually evaluated, it does not matter which one we use.
                Cluster const *tLeaderCluster = mSideSets( tSourceMeshIndex )->get_clusters_by_index( tSourceClusterIndex );
                uint const     tDummySetIndex = tSourceMeshIndex == 0 ? mSideSets.size() - 1 : 0;
                Cluster const *tDummyCluster  = mSideSets( tDummySetIndex )->get_clusters_on_set()( 0 );
                SetPair const  tSetPair       = std::make_pair( tSourceMeshIndex, tDummySetIndex );
                tNonconformalSideClusters[ tSetPair ].emplace_back(
                        tLeaderCluster,
                        tDummyCluster,
                        Vector< IntegrationPointPairs >{},
                        Vector< NodalPointPairs >{} );
                continue;    // skip the rest of the loop
            }
            for ( const auto &[ tTargetClusterLocator, tResultIndices ] : tTargetClusterToResultMap )
            {
                Vector< IntegrationPointPairs > tIntegrationPointPairs;    // bundles of integration points that were mapped from the leader side to the follower side cells
                Vector< NodalPointPairs >       tNodePointPairs;           // mapped points of all nodes on the leader side to the follower side
                // for this cluster-cluster pairing, get the pairing of their cells (using the indices, the access to the corresponding indices
                // in the mapping result is guaranteed to be correct)
                auto tCellPairing = extract_cell_pairing( aMappingResult, tResultIndices );
                for ( auto const &[ tCellPair, tCellResults ] : tCellPairing )
                {
                    populate_integration_and_nodal_point_pairs( aMappingResult, tIntegrationPointPairs, tNodePointPairs, tCellResults );
                }
                auto const &[ tTargetClusterIndex, tTargetMeshIndex ] = tTargetClusterLocator;
                SetPair const  tSetPair                               = std::make_pair( tSourceMeshIndex, tTargetMeshIndex );
                Cluster const *tLeaderCluster                         = mSideSets( tSourceMeshIndex )->get_clusters_by_index( tSourceClusterIndex );
                Cluster const *tFollowerCluster                       = mSideSets( tTargetMeshIndex )->get_clusters_by_index( tTargetClusterIndex );

                // append a new nonconformal side cluster
                tNonconformalSideClusters[ tSetPair ].emplace_back( tLeaderCluster, tFollowerCluster, tIntegrationPointPairs, tNodePointPairs );
            }
        }
        return tNonconformalSideClusters;
    }
    void Contact_Mesh_Editor::populate_integration_and_nodal_point_pairs(
            MappingResult const                      &aMappingResult,
            Vector< IntegrationPointPairs >          &aIntegrationPointPairs,
            Vector< NodalPointPairs >                &aNodePointPairs,
            Contact_Mesh_Editor::ResultIndices const &aCellResults ) const
    {
        Vector< moris_index > tNodalResultColumns;
        Vector< moris_index > tIntegrationPointResultColumns;
        for ( uint iResult = 0; iResult < aCellResults.size(); iResult++ )
        {
            size_t const tMappingResultColumn = aCellResults( iResult );
            if ( is_integration_point_result_index( tMappingResultColumn ) )
            {
                tIntegrationPointResultColumns.push_back( tMappingResultColumn );
            }
            else
            {
                tNodalResultColumns.push_back( tMappingResultColumn );
            }
        }

        if ( tIntegrationPointResultColumns.size() > 0 )
        {
            aIntegrationPointPairs.push_back( create_integration_point_pairs_from_results( tIntegrationPointResultColumns, aMappingResult ) );
        }
        if ( tNodalResultColumns.size() > 0 )
        {
            aNodePointPairs.push_back( create_nodal_point_pairs_from_results( tNodalResultColumns, aMappingResult ) );
        }
    }

    IntegrationPointPairs Contact_Mesh_Editor::create_integration_point_pairs_from_results( Vector< moris_index > aResultIndices, MappingResult aMappingResult ) const
    {
        size_t const tNumResults = aResultIndices.size();

        Matrix< DDRMat > tQWeights = mIntegrator.get_weights();
        Matrix< DDRMat > tQPoints  = mIntegrator.get_points();
        MORIS_ASSERT( tQPoints.n_rows() == 2 && sum( tQPoints.get_row( 1 ) ) < MORIS_REAL_EPS, "Currently, only 1D parametric coordinates with a constant time dimension are supported!" );

        // the leader and follower cell indices are the same for each result in this method!
        moris_index const tLeaderCellIndex   = aMappingResult.mSourceCellIndex( aResultIndices( 0 ) );
        moris_index const tFollowerCellIndex = aMappingResult.mTargetCellIndices( aResultIndices( 0 ) );
        Vector< real >    tWeights( tNumResults );
        Vector< real >    tDistances( tNumResults );
        Matrix< DDRMat >  tNormals( aMappingResult.mNormals.n_rows(), tNumResults );
        Matrix< DDRMat >  tReferenceNormals( aMappingResult.mReferenceNormals.n_rows(), tNumResults );
        Matrix< DDRMat >  tLeaderParametricCoords( tQPoints.n_rows(), tNumResults );
        Matrix< DDRMat >  tFollowerParametricCoords( tQPoints.n_rows(), tNumResults, -1.0 );

        // loop over each index in the mapping result for this pair of cells
        for ( uint iIndex = 0; iIndex < tNumResults; ++iIndex )
        {
            /* To get the correct column in the mapping result, we have to index into the list of result columns.
             * E.g. for a cell-cell pair, the mapping result columns are [ 3, 6, 8, 9 ] which means that the
             * mapping results of the first integration point (that got mapped successfully on the other cell) is stored in column 3, the second in column 6 and so on. */
            size_t const tMappingResultColumn = aResultIndices( iIndex );

            /* Depending on the index in the mapping result, the correct column in the list of integration points (tQPoints) has to be determined.
             * Since the integration points are repeated for each cell, the index is the remainder of the division by the number of integration points.
             * E.g we have 4 integration points per cell and 3 cells. The results in the mapping will refer to the integration points in columns
             * [ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1 ,2 ,3 ] and so on. */
            size_t const tIntegrationPointIndex = get_integration_point_index( tMappingResultColumn );
            tWeights( iIndex )                  = tQWeights( tIntegrationPointIndex );
            tDistances( iIndex )                = aMappingResult.mSignedDistance( tMappingResultColumn );
            tLeaderParametricCoords.set_column( iIndex, tQPoints.get_column( tIntegrationPointIndex ) );
            MORIS_ASSERT(
                    tLeaderCellIndex == aMappingResult.mSourceCellIndex( tMappingResultColumn ),
                    "Leader cell index does not match! This means that the pre-filtering of result indices did not work correctly!" );
            MORIS_ASSERT(
                    tFollowerCellIndex == aMappingResult.mTargetCellIndices( tMappingResultColumn ),
                    "Follower cell index does not match! This means that the pre-filtering of result indices did not work correctly!" );

            // Since the integration points have a constant time dimension, we can leave a -1.0 to the parametric coordinates.
            auto tCoordinate = aMappingResult.mTargetParametricCoordinate.get_column( tMappingResultColumn );
            for ( uint iCoord = 0; iCoord < tCoordinate.n_rows; ++iCoord )
            {
                tFollowerParametricCoords( iCoord, iIndex ) = tCoordinate( iCoord );
            }
            tNormals.set_column( iIndex, aMappingResult.mNormals.get_column( tMappingResultColumn ) );
            tReferenceNormals.set_column( iIndex, aMappingResult.mReferenceNormals.get_column( tMappingResultColumn ) );
        }

        return {
            tLeaderCellIndex,
            tLeaderParametricCoords,
            tFollowerCellIndex,
            tFollowerParametricCoords,
            tWeights,
            tDistances,
            tNormals,
            tReferenceNormals
        };
    }


    NodalPointPairs Contact_Mesh_Editor::create_nodal_point_pairs_from_results( Vector< moris_index > aResultIndices, MappingResult aMappingResult ) const
    {
        size_t const tNumResults = aResultIndices.size();

        // the leader and follower cluster and cell indices are the same for each result in this method!
        moris_index const tLeaderClusterIndex = aMappingResult.mSourceClusterIndex( aResultIndices( 0 ) );
        moris_index const tLeaderCellIndex    = aMappingResult.mSourceCellIndex( aResultIndices( 0 ) );
        moris_index const tFollowerCellIndex  = aMappingResult.mTargetCellIndices( aResultIndices( 0 ) );

        Vector< real >        tDistances( tNumResults );
        Vector< moris_index > tLeaderNodeIndices( tNumResults );
        Matrix< DDRMat >      tNormals( aMappingResult.mNormals.n_rows(), tNumResults );
        Matrix< DDRMat >      tReferenceNormals( aMappingResult.mReferenceNormals.n_rows(), tNumResults );
        Matrix< DDRMat >      tFollowerParametricCoords( mIntegrator.get_points().n_rows(), tNumResults, -1.0 );

        auto const            &tSideSet         = mSideSets( aMappingResult.mSourceMeshIndex );
        auto const            &tCluster         = dynamic_cast< mtk::Side_Cluster const                    *>( tSideSet->get_clusters_by_index( tLeaderClusterIndex ) );
        Matrix< IdMat >        tSideOrdinals    = tCluster->get_cell_side_ordinals();
        Vector< const Cell * > tCells           = tCluster->get_primary_cells_in_cluster();
        auto                   tHasCorrectIndex = [ &tLeaderCellIndex ]( const Cell *const &aCell ) { return aCell->get_index() == tLeaderCellIndex; };
        auto const            &tCell            = std::find_if( tCells.begin(), tCells.end(), tHasCorrectIndex );
        MORIS_ASSERT( tCell != tCells.end(), "Contact_Mesh_Editor::create_nodal_point_pairs_from_results: Could not find cell with index %d in cluster %d!", tLeaderCellIndex, tLeaderClusterIndex );
        size_t      tCellIndex = std::distance( tCells.begin(), tCell );
        auto const &tVertices  = tCells( tCellIndex )->get_vertices_on_side_ordinal( tSideOrdinals( tCellIndex ) );

        // loop over each index in the mapping result for this pair of cells
        for ( uint iIndex = 0; iIndex < tNumResults; ++iIndex )
        {
            /* To get the correct column in the mapping result, we have to index into the list of result columns.
             * E.g. for a cell-cell pair, the mapping result columns are [ 3, 6, 8, 9 ] which means that the
             * mapping results of the first integration point (that got mapped successfully on the other cell) is stored in column 3, the second in column 6 and so on. */
            size_t const tMappingResultColumn = aResultIndices( iIndex );

            /* The next line provides the local index of the node that has been mapped in the result. E.g. if we have a 2 node line cell, the two nodal points
             * will correspond to the node at the parametric coordinate (-1.0) and (1.0) at the local indices 0 and 1 respectively.
             * It is up to this function to determine the global index of the node in the mesh! */
            size_t const tNodeIndex      = get_node_coordinate_index( tMappingResultColumn );
            tLeaderNodeIndices( iIndex ) = tVertices( tNodeIndex )->get_index();

            // Since the integration points have a constant time dimension, we can leave a -1.0 to the parametric coordinates.
            auto tCoordinate = aMappingResult.mTargetParametricCoordinate.get_column( tMappingResultColumn );
            for ( uint iCoord = 0; iCoord < tCoordinate.n_rows; ++iCoord )
            {
                tFollowerParametricCoords( iCoord, iIndex ) = tCoordinate( iCoord );
            }
            tDistances( iIndex ) = aMappingResult.mSignedDistance( tMappingResultColumn );
            tNormals.set_column( iIndex, aMappingResult.mNormals.get_column( tMappingResultColumn ) );
            tReferenceNormals.set_column( iIndex, aMappingResult.mReferenceNormals.get_column( tMappingResultColumn ) );
        }

        return {
            tLeaderCellIndex,
            tLeaderNodeIndices,
            tFollowerCellIndex,
            tFollowerParametricCoords,
            tDistances,
            tNormals,
            tReferenceNormals
        };
    }


    Vector< MappingResult > Contact_Mesh_Editor::perform_mapping( Matrix< DDRMat > aPointsToMap ) const
    {
        // get all possible source side sets that have been specified in the candidate pairings
        std::set< moris_index > tSourceSideSets;
        std::transform(
                mCandidatePairs.begin(),
                mCandidatePairs.end(),
                std::inserter( tSourceSideSets, tSourceSideSets.end() ),
                []( auto const &aPair ) { return aPair.first; } );

        //        Json                    tMappingResultsJson; // TODO @ff: Remove! Only for debugging!
        Vector< MappingResult > tMappingResults;
        tMappingResults.reserve( tSourceSideSets.size() );

        for ( auto const &tSourceSideSet : tSourceSideSets )
        {
            MappingResult tResult = mPointMapper.map( tSourceSideSet, aPointsToMap, mMaxNegativeRayLength, mMaxPositiveRayLength );
            tMappingResults.push_back( tResult );
            //            tMappingResultsJson.put_child( mSideSets( tSourceSideSet )->get_set_name(), tResult.to_json() ); // TODO @ff: Remove! Only for debugging!
        }
        // TODO @ff: Remove! Only for debugging!
        //        uint const  tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
        //        std::string tFileName  = "mapping_result_" + std::to_string( tIteration ) + ".json";
        //        write_json( tFileName, tMappingResultsJson );

        return tMappingResults;
    }

    void Contact_Mesh_Editor::update_nonconformal_side_sets() const
    {
        // removes all nonconformal sidesets and clusters from the IGMesh
        mIGMesh->reset_nonconformal_side_set();


        // each MappingResult will contain the mapping results for one source (leader) side set to many target (follower) side sets
        Vector< MappingResult > const tMappingResults = perform_mapping( get_points_to_map() );

        // we convert each mapping result to a list of nonconformal side clusters, grouped by the source- and target side sets
        std::map< SetPair, Vector< Nonconformal_Side_Cluster > > tConvertedResults;
        for ( auto const &tMappingResult : tMappingResults )
        {
            tConvertedResults.merge( convert_mapping_result_to_nonconformal_side_clusters( tMappingResult ) );
        }

        // get the total number of nonconformal side clusters to reserve the correct amount of memory
        size_t const tNumNonconformalSideClusters = std::accumulate(
                tConvertedResults.begin(),
                tConvertedResults.end(),
                0,
                []( size_t aSum, auto const &aPair ) { return aSum + aPair.second.size(); } );

        // std::cout << "Number of nonconformal side clusters: " << tNumNonconformalSideClusters << std::endl;

        // the number of nonconformal side clusters that get stored in the IGMesh has to be reserved beforehand.
        // If a reallocation of memory would be necessary, the pointers to the nonconformal side clusters would be invalidated, which would
        // lead to wrong references being stored in the nonconformal side sets.
        mIGMesh->reserve_nonconformal_side_clusters( tNumNonconformalSideClusters );
        // each group of unique source- and target side set-pairs will be stored in a nonconformal side set
        for ( auto const &[ tSetPair, tNonconformalSideClusters ] : tConvertedResults )
        {
            mIGMesh->add_nonconformal_side_set(
                    this->get_nonconformal_side_set_name( tSetPair ),
                    tNonconformalSideClusters,
                    mSideSets( tSetPair.first )->get_set_colors()    // TODO: Is this correct? Just using the color of the first set...
            );
        }
    }

    std::string Contact_Mesh_Editor::get_nonconformal_side_set_name( Contact_Mesh_Editor::SetPair const &tSetPair ) const
    {
        // as the bulk phase information is not available at this point of the code,
        // we have to use the side set names instead to extract the correct bulk phase information.
        std::string const &tSourceMeshName = mSideSets( tSetPair.first )->get_set_name();
        std::string const &tTargetMeshName = mSideSets( tSetPair.second )->get_set_name();

        auto const &[ tSourceB0, tSourceB1 ] = get_leaderphase_from_set_name( tSourceMeshName );
        auto const &[ tTargetB0, tTargetB1 ] = get_leaderphase_from_set_name( tTargetMeshName );

        return "ncss|iside_b0_" + tSourceB0 + "_b1_" + tSourceB1 + "|iside_b0_" + tTargetB0 + "_b1_" + tTargetB1;
    }

    std::pair< std::string, std::string > Contact_Mesh_Editor::get_leaderphase_from_set_name( std::string const &aSideSetName )
    {
        // the format of the side set names is "iside_b0_<leader_phase>_b1_<follower_phase>"
        size_t const       tB0Pos   = aSideSetName.find( "_b0_" );
        size_t const       tB1Pos   = aSideSetName.find( "_b1_" );
        std::string const &tB0Phase = aSideSetName.substr( tB0Pos + 4, tB1Pos - tB0Pos - 4 );
        std::string const &tB1Phase = aSideSetName.substr( tB1Pos + 4 );
        return { tB0Phase, tB1Phase };
    }

    void Contact_Mesh_Editor::update_displacements( std::unordered_map< moris_index, Vector< real > > const &aNodalDisplacements )
    {
        mPointMapper.update_displacements( aNodalDisplacements );
    }

    Vector< Side_Set const * > Contact_Mesh_Editor::get_side_sets() const
    {
        return mSideSets;
    }

    Matrix< DDRMat > Contact_Mesh_Editor::get_nodal_parametric_coordinates() const
    {
        auto const &tCluster      = dynamic_cast< mtk::Side_Cluster const      *>( mSideSets( 0 )->get_clusters_on_set()( 0 ) );
        auto const &tCell         = tCluster->get_primary_cells_in_cluster()( 0 );
        auto const &tSideOrdinals = tCluster->get_cell_side_ordinals();
        auto const &tVertices     = tCell->get_vertices_on_side_ordinal( tSideOrdinals( 0 ) );
        MORIS_ASSERT( tVertices.size() == 2, "Currently, only Line elements are supported in the nonconformal mapping!" );
        return { { -1.0, 1.0 } };    // has to be generalized for different cell types
    }

    Matrix< DDRMat > Contact_Mesh_Editor::get_points_to_map() const
    {
        Matrix< DDRMat > const tNodalParametricCoordinates = get_nodal_parametric_coordinates();
        Matrix< DDRMat > const tIntegrationPoints          = mIntegrator.get_points().get_row( 0 );    // only use the spatial parametric coordinates, not the time dimension
        return join_horiz( tNodalParametricCoordinates, tIntegrationPoints );
    }


    bool Contact_Mesh_Editor::is_integration_point_result_index( moris_index aMappingResultColumnIndex ) const
    {
        size_t const nNodalPoints       = get_nodal_parametric_coordinates().n_cols();
        size_t const nIntegrationPoints = mIntegrator.get_number_of_points();
        size_t const nLocalIndex        = aMappingResultColumnIndex % ( nNodalPoints + nIntegrationPoints );
        return nLocalIndex >= nNodalPoints;
    }

    moris_index Contact_Mesh_Editor::get_integration_point_index( moris_index aMappingResultColumnIndex ) const
    {
        size_t const nNodalPoints = get_nodal_parametric_coordinates().n_cols();
        return aMappingResultColumnIndex % ( nNodalPoints + mIntegrator.get_number_of_points() ) - nNodalPoints;
    }

    moris_index Contact_Mesh_Editor::get_node_coordinate_index( moris_index aMappingResultColumnIndex ) const
    {
        size_t const nNodalPoints = get_nodal_parametric_coordinates().n_cols();
        return aMappingResultColumnIndex % ( nNodalPoints + mIntegrator.get_number_of_points() );
    }


}    // namespace moris::mtk
