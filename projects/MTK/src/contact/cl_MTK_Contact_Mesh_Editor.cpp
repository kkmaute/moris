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
#include "cl_MTK_MappingResult.hpp"
#include "cl_MTK_IntegrationPointPairs.hpp"
#include "cl_MTK_Nonconformal_Side_Cluster.hpp"

namespace moris::mtk
{
    std::map< Contact_Mesh_Editor::ClusterPair, std::map< Contact_Mesh_Editor::CellPair, Contact_Mesh_Editor::ResultIndices > >
    Contact_Mesh_Editor::extract_cluster_and_cell_pairing( MappingResult const &aResult )
    {

        std::map< ClusterPair, std::map< CellPair, ResultIndices > > tClusterPairs;

        for ( moris_index i = 0; i < static_cast< moris_index >( aResult.mSourceClusterIndex.size() ); ++i )
        {
            if ( aResult.mTargetClusterIndex( i ) == -1 )
            {
                continue;    // skip unsuccessful mappings
            }
            moris_index const &tSourceClusterIndex = aResult.mSourceClusterIndex( i );
            moris_index const &tTargetClusterIndex = aResult.mTargetClusterIndex( i );
            moris_index const &tTargetSideSetIndex = aResult.mTargetSideSetIndices( i );
            moris_index const &tSourceCellIndex    = aResult.mSourceCellIndex( i );
            moris_index const &tTargetCellIndex    = aResult.mTargetCellIndices( i );

            auto const tClusterPair = std::make_tuple( tSourceClusterIndex, tTargetClusterIndex, tTargetSideSetIndex );
            auto const tCellPair    = std::make_pair( tSourceCellIndex, tTargetCellIndex );
            tClusterPairs[ tClusterPair ][ tCellPair ].push_back( i );
        }
        return tClusterPairs;
    }

    moris::Cell< Nonconformal_Side_Cluster >
    Contact_Mesh_Editor::convert_mapping_result_to_nonconformal_side_clusters( moris_index aSourceSideSetIndex, MappingResult aResult )
    {
        Matrix< DDRMat > tQWeights;
        Matrix< DDRMat > tQPoints;
        mIntegrator.get_weights( tQWeights );
        mIntegrator.get_points( tQPoints );
        tQPoints = tQPoints.get_row( 0 ).eval();    // for some reason, the dimension of the matrix is one too high

        size_t const tNumIntegrationPoints = mIntegrator.get_number_of_points();

        // extract the cluster pairs and the corresponding cell pairs from the mapping result this will act as a basis to create the nonconformal side clusters
        auto tClusterPairs = extract_cluster_and_cell_pairing( aResult );

        // Since we know the number of unique cluster-pairs, we can reserve the correct amount of memory for the nonconformal side clusters
        moris::Cell< Nonconformal_Side_Cluster > tNonconformalSideClusters;
        tNonconformalSideClusters.reserve( tClusterPairs.size() );

        // loop over each unique cluster pair which itself has multiple cell-cell pairs
        for ( const auto &[ tClusterPair, tCellMap ] : tClusterPairs )
        {
            // the integration point pairs will store the bundles of integration points that were mapped from the follower side to the leader side cells
            moris::Cell< IntegrationPointPairs > tIntegrationPointPairs;

            // the cell map contains all pairs of source and target cells that were mapped onto each other
            // the values of this map are the list of indices to get access to the correct entries of the mapping result
            for ( auto const &[ tCellMap, tResultColumns ] : tCellMap )
            {
                moris::Cell< real > tWeights( tResultColumns.size() );
                Matrix< DDRMat >    tNormals( aResult.mNormal.n_rows(), tResultColumns.size() );
                Matrix< DDRMat >    tFollowerParametricCoords( tQPoints.n_rows(), tResultColumns.size() );
                Matrix< DDRMat >    tLeaderParametricCoords( tQPoints.n_rows(), tResultColumns.size() );

                // loop over each index in the mapping result for this pair of cells
                for ( uint iIndex = 0; iIndex < tResultColumns.size(); ++iIndex )
                {
                    /* Depending on the index in the mapping result, the correct column in the list of integration points (tQPoints) has to be determined.
                     * Since the integration points are repeated for each cell, the index is the remainder of the division by the number of integration points.
                     * E.g we have 4 integration points per cell and 3 cells. The results in the mapping will refer to the integration points in columns
                     * [ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1 ,2 ,3 ] and so on. */
                    size_t tIntegrationPointColumn = iIndex % tNumIntegrationPoints;
                    tWeights( iIndex )             = tQWeights( tIntegrationPointColumn );
                    tFollowerParametricCoords.set_column( iIndex, tQPoints.get_column( tIntegrationPointColumn ) );

                    /*To get the correct column in the mapping result, we have to index into the list of result columns.
                     * E.g. for this cell-cell pair, the mapping result columns are [ 3, 6, 8, 9 ] which means that the
                     * mapping results of the first integration point is stored in column 3, the second in column 6 and so on. */
                    size_t const tMappingResultColumn = tResultColumns( iIndex );
                    tLeaderParametricCoords.set_column( iIndex, aResult.mTargetParametricCoordinate.get_column( tMappingResultColumn ) );
                    tNormals.set_column( iIndex, aResult.mNormal.get_column( tMappingResultColumn ) );
                }

                auto const &[ iSourceCell, iTargetCell ] = tCellMap;
                tIntegrationPointPairs.emplace_back(
                        iSourceCell,
                        tFollowerParametricCoords,
                        iTargetCell,
                        tLeaderParametricCoords,
                        tWeights,
                        tNormals );
            }

            auto const &[ tSourceCluster, tTargetCluster, tTargetMesh ] = tClusterPair;

            Cluster const *tFollowerCluster = mSideSets( aSourceSideSetIndex )->get_clusters_by_index( tSourceCluster );
            Cluster const *tLeaderCluster   = mSideSets( tTargetMesh )->get_clusters_by_index( tTargetCluster );
            //
            tNonconformalSideClusters.emplace_back( tFollowerCluster, tLeaderCluster, tIntegrationPointPairs );
        }

        return tNonconformalSideClusters;
    }


    void Contact_Mesh_Editor::update_nonconformal_side_sets()
    {
        Matrix< DDRMat > tQPoints;
        mIntegrator.get_points( tQPoints );
        tQPoints = tQPoints.get_row( 0 ).eval();

        Matrix< DDRMat > tQWeights;
        mIntegrator.get_weights( tQWeights );

        // get all possible source side sets that have been specified in the candidate pairings
        std::set< moris_index > tSourceSideSets;
        std::transform(
                mCandidatePairs.begin(),
                mCandidatePairs.end(),
                std::inserter( tSourceSideSets, tSourceSideSets.end() ),
                []( auto const &aPair ) { return aPair.first; } );

        // clear the nonconformal side clusters in the IG mesh database
        mIGMesh->mNonconformalSideClusters.clear();
        Json tMappingResults;
        for ( auto const &tSourceSideSet : tSourceSideSets )
        {
            MappingResult tResult = mPointMapper.map( tSourceSideSet, tQPoints );
            tMappingResults.put_child( mSideSets( tSourceSideSet )->get_set_name(), tResult.to_json() );

            auto tNonconformalSideClusters = convert_mapping_result_to_nonconformal_side_clusters( tSourceSideSet, tResult );
            mIGMesh->mNonconformalSideClusters.append( tNonconformalSideClusters );
        }
        write_json( "mapping_result.json", tMappingResults );
    }

    void Contact_Mesh_Editor::update_displacements( Matrix< DDRMat > &aDisplacements ) {}

}    // namespace moris::mtk
