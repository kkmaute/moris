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
#include <algorithm>
#include <iterator>

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

    std::map< Contact_Mesh_Editor::SetPair, Vector< Nonconformal_Side_Cluster > >
    Contact_Mesh_Editor::convert_mapping_result_to_nonconformal_side_clusters( MappingResult const &aMappingResult ) const
    {
        Matrix< DDRMat > tQWeights             = mIntegrator.get_weights();
        Matrix< DDRMat > tQPoints              = mIntegrator.get_points().get_row( 0 );    // for some reason, the dimension of the matrix is one too high
        size_t const     tNumIntegrationPoints = mIntegrator.get_number_of_points();

        // extract the cluster pairs and the corresponding cell pairs from the mapping result this will act as a basis to create the nonconformal side clusters
        auto tClusterPairs = extract_cluster_and_cell_pairing( aMappingResult );

        // Since we know the number of unique cluster-pairs, we can reserve the correct amount of memory for the nonconformal side clusters
        std::map< SetPair, Vector< Nonconformal_Side_Cluster > > tNonconformalSideClusters;

        // loop over each unique cluster pair which itself has multiple cell-cell pairs
        for ( const auto &[ tClusterPair, tCellMaps ] : tClusterPairs )
        {
            // the integration point pairs will store the bundles of integration points that were mapped from the follower side to the leader side cells
            Vector< IntegrationPointPairs > tIntegrationPointPairs;

            // the cell map contains all pairs of source and target cells that were mapped onto each other
            // the values of this map are the list of indices to get access to the correct entries of the mapping result
            for ( auto const &[ tCellMap, tResultColumns ] : tCellMaps )
            {
                Vector< real >   tWeights( tResultColumns.size() );
                Matrix< DDRMat > tNormals( aMappingResult.mNormal.n_rows(), tResultColumns.size() );
                Matrix< DDRMat > tFollowerParametricCoords( tQPoints.n_rows(), tResultColumns.size() );
                Matrix< DDRMat > tLeaderParametricCoords( tQPoints.n_rows(), tResultColumns.size() );

                // loop over each index in the mapping result for this pair of cells
                for ( uint iIndex = 0; iIndex < tResultColumns.size(); ++iIndex )
                {
                    /* Depending on the index in the mapping result, the correct column in the list of integration points (tQPoints) has to be determined.
                     * Since the integration points are repeated for each cell, the index is the remainder of the division by the number of integration points.
                     * E.g we have 4 integration points per cell and 3 cells. The results in the mapping will refer to the integration points in columns
                     * [ 0, 1, 2, 3, 0, 1, 2, 3, 0, 1 ,2 ,3 ] and so on. */
                    size_t const tIntegrationPointColumn = iIndex % tNumIntegrationPoints;
                    tWeights( iIndex )                   = tQWeights( tIntegrationPointColumn );
                    tFollowerParametricCoords.set_column( iIndex, tQPoints.get_column( tIntegrationPointColumn ) );

                    /*To get the correct column in the mapping result, we have to index into the list of result columns.
                     * E.g. for this cell-cell pair, the mapping result columns are [ 3, 6, 8, 9 ] which means that the
                     * mapping results of the first integration point is stored in column 3, the second in column 6 and so on. */
                    size_t const tMappingResultColumn = tResultColumns( iIndex );
                    tLeaderParametricCoords.set_column( iIndex, aMappingResult.mTargetParametricCoordinate.get_column( tMappingResultColumn ) );
                    tNormals.set_column( iIndex, aMappingResult.mNormal.get_column( tMappingResultColumn ) );
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

            auto const &[ tSourceClusterIndex, tTargetClusterIndex, tTargetMeshIndex ] = tClusterPair;

            moris_index const tSourceMeshIndex = aMappingResult.mSourceMeshIndex;

            SetPair const tSetPair = std::make_pair( tSourceMeshIndex, tTargetMeshIndex );

            Cluster const *tFollowerCluster = mSideSets( tSourceMeshIndex )->get_clusters_by_index( tSourceClusterIndex );
            Cluster const *tLeaderCluster   = mSideSets( tTargetMeshIndex )->get_clusters_by_index( tTargetClusterIndex );

            tNonconformalSideClusters[ tSetPair ].emplace_back( tFollowerCluster, tLeaderCluster, tIntegrationPointPairs );
        }

        return tNonconformalSideClusters;
    }


    Vector< MappingResult > Contact_Mesh_Editor::perform_mapping() const
    {
        // get all possible source side sets that have been specified in the candidate pairings
        std::set< moris_index > tSourceSideSets;
        std::transform(
                mCandidatePairs.begin(),
                mCandidatePairs.end(),
                std::inserter( tSourceSideSets, tSourceSideSets.end() ),
                []( auto const &aPair ) { return aPair.first; } );

        Json                    tMappingResultsJson;
        Vector< MappingResult > tMappingResults;
        tMappingResults.reserve( tSourceSideSets.size() );

        for ( auto const &tSourceSideSet : tSourceSideSets )
        {
            Matrix< DDRMat > const tQuadPoints = mIntegrator.get_points().get_row( 0 );    // TODO: for some reason, the parametric coordinates are 2d instead of 1d
            MappingResult          tResult     = mPointMapper.map( tSourceSideSet, tQuadPoints );
            tMappingResults.push_back( tResult );
            tMappingResultsJson.put_child( mSideSets( tSourceSideSet )->get_set_name(), tResult.to_json() );
        }

        write_json( "mapping_result.json", tMappingResultsJson );

        return tMappingResults;
    }

    void Contact_Mesh_Editor::update_nonconformal_side_sets() const
    {
        // removes all nonconformal sidesets and clusters from the IGMesh
        mIGMesh->reset_nonconformal_side_set();

        // each mapping result will contain the mapping results for one source side set to many target side sets
        Vector< MappingResult > const tMappingResults = perform_mapping();

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

        std::cout << "Number of nonconformal side clusters: " << tNumNonconformalSideClusters << std::endl;

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

        return "ncss_" + tSourceB0 + "_" + tSourceB1 + "_to_" + tTargetB0 + "_" + tTargetB1;
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

    void Contact_Mesh_Editor::update_displacements( Matrix< DDRMat > &aDisplacements ) {}

}    // namespace moris::mtk
