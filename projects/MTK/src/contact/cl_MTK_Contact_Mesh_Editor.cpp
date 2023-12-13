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

namespace moris::mtk
{
    void Contact_Mesh_Editor::update_nonconformal_side_sets()
    {
        Matrix< DDRMat > tQPoints;
        mIntegrationRule.get_points( tQPoints );
        tQPoints = tQPoints.get_row( 0 ).eval();

        Matrix< DDRMat > tQWeights;
        mIntegrationRule.get_weights( tQWeights );


        MappingResult tResult = mPointMapper.map( 0, tQPoints );
        write_json( "mapping_result.json", tResult.to_json() );


        /**
         * TODO: Can we make it cleaner?
         * This (rather/quite) complex data structure is used to store the correct indices of the mapping results for each unique pair of clusters (temporarily).
         * Each unique pair of clusters is identified by a tuple of the source cluster index, the target cluster index and the target side set index.
         * As value, a map is stored that maps the source-target cell pairs to the list of indices to get access to the correct entries of the mapping result in a later step.
         * As we do not know in advance how many integration points will be mapped onto each of the target cells, we first have to structure the plain indices into the
         * result map. Afterwards, we can create matrices in the correct size and fill them with the correct values according to the collected indices.
         *
         * (SourceCluster, TargetCluster, TargetMesh) -- n -> (SourceCell, TargetCell) -- n -> (Index in MappingResult, e.g. to access integration points)
         */
        std::map<
                std::tuple< moris_index, moris_index, moris_index >,    // (SourceCluster, TargetCluster, TargetMesh)
                std::map<
                        std::pair< moris_index, moris_index >,    // (SourceCell, TargetCell)
                        moris::Cell< moris_index >                // List of columns in the mapping result that contain coordinates for this pair of cells
                        > >
                tClusterPairs;

        for ( moris_index i = 0; i < static_cast< moris_index >( tResult.mSourceClusterIndex.size() ); ++i )
        {
            if ( tResult.mTargetClusterIndex( i ) == -1 )
            {
                continue;    // skip unsuccessful mappings
            }
            moris_index const &tSourceClusterIndex = tResult.mSourceClusterIndex( i );
            moris_index const &tTargetClusterIndex = tResult.mTargetClusterIndex( i );
            moris_index const &tTargetSideSetIndex = tResult.mTargetSideSetIndices( i );
            moris_index const &tSourceCellIndex    = tResult.mSourceCellIndex( i );
            moris_index const &tTargetCellIndex    = tResult.mTargetCellIndices( i );

            auto const tClusterPair = std::make_tuple( tSourceClusterIndex, tTargetClusterIndex, tTargetSideSetIndex );
            auto const tCellPair    = std::make_pair( tSourceCellIndex, tTargetCellIndex );
            tClusterPairs[ tClusterPair ][ tCellPair ].push_back( i );
        }


        moris::Cell< Nonconformal_Side_Cluster > tNonconformalSideClusters;
        tNonconformalSideClusters.reserve( tClusterPairs.size() );
        size_t const tNumIntegrationPoints = mIntegrationRule.get_number_of_points();
        for ( const auto &[ tClusterPair, tCellMap ] : tClusterPairs )
        {
            auto const &[ tSourceCluster, tTargetCluster, tTargetMesh ] = tClusterPair;
            moris::Cell< IntegrationPointPair > tIntegrationPointPairs;

            // the cell map contains all pairs of source and target cells that were mapped onto each other
            // the values of this map are the list of indices to get access to the correct entries of the mapping result
            for ( auto const &[ tCellMap, tResultColumns ] : tCellMap )
            {
                moris::Cell< real > tWeights( tResultColumns.size() );
                Matrix< DDRMat >    tNormals( tResult.mNormal.n_rows(), tResultColumns.size() );
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

                    //
                    size_t tMappingResultColumn = tResultColumns( iIndex );
                    tLeaderParametricCoords.set_column( iIndex, tResult.mTargetParametricCoordinate.get_column( tMappingResultColumn ) );
                    tNormals.set_column( iIndex, tResult.mNormal.get_column( tMappingResultColumn ) );
                }

                auto const &[ iSourceCell, iTargetCell ] = tCellMap;
                tIntegrationPointPairs.emplace_back( iSourceCell, tFollowerParametricCoords, iTargetCell, tLeaderParametricCoords, tWeights, tNormals );
            }

            Cluster const *tFollowerCluster = mSideSets( 0 )->get_clusters_by_index( tSourceCluster );    // TODO: Fix SourceMeshIndex
            Cluster const *tLeaderCluster   = mSideSets( tTargetMesh )->get_clusters_by_index( tTargetCluster );
            //
            tNonconformalSideClusters.emplace_back( tFollowerCluster, tLeaderCluster, tIntegrationPointPairs );
        }


        //        MappingResult tResult_1 = mPointMapper.map( 1, tQPoints );
        //        write_json( "mapping_result_1.json", tResult_1.to_json() );
    }

    void Contact_Mesh_Editor::update_displacements( Matrix< DDRMat > &aDisplacements ) {}

}    // namespace moris::mtk
