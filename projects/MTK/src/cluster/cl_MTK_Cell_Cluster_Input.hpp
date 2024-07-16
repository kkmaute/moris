/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Cluster_Input.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_

#include "cl_Vector.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include <unordered_map>
#include <iomanip>    // std::setw

namespace moris
{

    namespace mtk
    {
        ///////////////////////////////////
        // STRUC FOR CELL CLUSTER INPUT  //
        ///////////////////////////////////

        struct Cell_Cluster_Input
        {
            Cell_Cluster_Input()
            {
            }

            void
            add_cluster_data( mtk::Cell             *aInterpCell,
                    moris::Matrix< IndexMat > const *aPrimaryCellIds,
                    moris::Matrix< IndexMat > const *aVoidCellIds,
                    moris::Matrix< IndexMat > const *aVerticesInCluster,
                    moris::Matrix< DDRMat > const   *aLocalCoordsRelativeToInterpCell )
            {
                MORIS_ASSERT( aVerticesInCluster->numel() == aLocalCoordsRelativeToInterpCell->n_rows(), "Number of vertices in the cluster must match the number of rows in local coord mat" );

                moris::uint tIndex = mInterpolationCells.size();

                // add to map
                MORIS_ASSERT( mInterpCellIndexToIndex.find( aInterpCell->get_id() ) == mInterpCellIndexToIndex.end(), "Trying to add cluster on an interpolation cell which already has a cluster." );
                mInterpCellIndexToIndex[ aInterpCell->get_id() ] = tIndex;

                mInterpolationCells.push_back( aInterpCell );
                mPrimaryCellIds.push_back( aPrimaryCellIds );
                mVoidCellsIds.push_back( aVoidCellIds );
                mVerticesIdsInCluster.push_back( aVerticesInCluster );
                mLocalCoordsRelativeToInterpCell.push_back( aLocalCoordsRelativeToInterpCell );
            }

            moris_index
            get_cluster_index( moris_id aInterpCellId )
            {
                auto tIter = mInterpCellIndexToIndex.find( aInterpCellId );

                if ( tIter == mInterpCellIndexToIndex.end() )
                {
                    return MORIS_INDEX_MAX;
                }

                else
                {
                    return tIter->second;
                }
            }

            moris::uint
            get_num_cell_clusters()
            {
                return mInterpolationCells.size();
            }

            mtk::Cell *
            get_interp_cell( moris_index aClusterIndex )
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mInterpolationCells( aClusterIndex );
            }

            moris::Matrix< IndexMat > const *
            get_primary_cell_ids( moris_index aClusterIndex )
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mPrimaryCellIds( aClusterIndex );
            }

            moris::Matrix< IndexMat > const *
            get_void_cell_ids( moris_index aClusterIndex )
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mVoidCellsIds( aClusterIndex );
            }

            moris::Matrix< IndexMat > const *
            get_vertex_in_cluster_ids( moris_index aClusterIndex )
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mVerticesIdsInCluster( aClusterIndex );
            }

            moris::Matrix< DDRMat > const *
            get_vertex_local_coords_wrt_interpolation_cell( moris_index aClusterIndex )
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mLocalCoordsRelativeToInterpCell( aClusterIndex );
            }

            void
            print()
            {
                std::cout << std::left << std::setw( 10 ) << "Interpolation Cells with Cluster Information\n";

                std::cout << "Num Interp Cells W Clusters: " << mInterpolationCells.size() << std::endl;
                std::cout << "Num Primary Integ Cells:     " << mPrimaryCellIds.size() << std::endl;
                std::cout << "Num Void Integ Cells:        " << mVoidCellsIds.size() << std::endl;
                for ( moris::uint i = 0; i < mInterpolationCells.size(); i++ )
                {
                    std::cout << std::left << "Interpolation Cell Id: " << std::setw( 6 ) << mInterpolationCells( i )->get_id() << std::endl;

                    // print primary information
                    std::cout << std::left << "     Primary Integration Cell Ids: ";
                    for ( moris::uint j = 0; j < mPrimaryCellIds( i )->numel(); j++ )
                    {
                        std::cout << std::right << std::setw( 6 ) << ( *mPrimaryCellIds( i ) )( j ) << " ";
                    }

                    // print void information
                    std::cout << std::left << "\n     Void Integration Cell Ids:    ";
                    for ( moris::uint j = 0; j < mVoidCellsIds( i )->numel(); j++ )
                    {
                        std::cout << std::right << std::setw( 6 ) << ( *mVoidCellsIds( i ) )( j ) << " ";
                    }

                    std::cout << "\n      Vertex Parametric Coordinates:" << std::endl;
                    for ( moris::uint k = 0; k < mLocalCoordsRelativeToInterpCell( i )->n_rows(); k++ )
                    {
                        std::cout << "            Vert Id: " << std::right << std::setw( 5 ) << ( *mVerticesIdsInCluster( i ) )( k ) << " | ";
                        for ( moris::uint j = 0; j < mLocalCoordsRelativeToInterpCell( i )->n_cols(); j++ )
                        {
                            std::cout << std::scientific << std::right << std::setw( 14 ) << ( *mLocalCoordsRelativeToInterpCell( i ) )( k, j ) << "   ";
                        }
                        std::cout << std::endl;
                    }
                }
            }

          private:
            Vector< mtk::Cell * >                       mInterpolationCells;
            Vector< moris::Matrix< IndexMat > const * > mPrimaryCellIds;
            Vector< moris::Matrix< IndexMat > const * > mVoidCellsIds;
            Vector< moris::Matrix< IndexMat > const * > mVerticesIdsInCluster;
            Vector< moris::Matrix< DDRMat > const * >   mLocalCoordsRelativeToInterpCell;

            std::unordered_map< moris_index, moris_index > mInterpCellIndexToIndex;
        };

    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_CELL_CLUSTER_INPUT_HPP_ */
