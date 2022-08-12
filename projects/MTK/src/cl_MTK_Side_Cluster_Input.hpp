/*
 * cl_MTK_Side_Cluster_Input.hpp
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_INPUT_HPP_

#include "cl_Cell.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_Matrix.hpp"
#include <unordered_map>

namespace moris
{
    namespace mtk
    {
        class Side_Set_Cluster_Data
        {
          public:
            Side_Set_Cluster_Data() {}

            void
            add_cluster_data(
                    bool                             aTrivial,
                    mtk::Cell                       *aInterpCell,
                    moris::Matrix< IndexMat > const *aCellIdsAndSideOrds,
                    moris::Matrix< IndexMat > const *aVerticesInCluster,
                    moris::Matrix< DDRMat > const   *aLocalCoordsRelativeToInterpCell )
            {
                moris::uint tIndex = mInterpolationCells.size();

                // add to map
                MORIS_ASSERT( mInterpCellIndexToIndex.find( aInterpCell->get_id() ) == mInterpCellIndexToIndex.end(),
                        "Trying to add cluster on an interpolation cell which already has a cluster." );

                mInterpCellIndexToIndex[ aInterpCell->get_id() ] = tIndex;

                if ( aTrivial )
                {
                    mTrivialFlag.push_back( 1 );
                }
                else
                {
                    mTrivialFlag.push_back( 0 );
                }

                mInterpolationCells.push_back( aInterpCell );
                mCellIdsAndSideOrdinals.push_back( aCellIdsAndSideOrds );
                mVerticesIdsInCluster.push_back( aVerticesInCluster );
                mLocalCoordsRelativeToInterpCell.push_back( aLocalCoordsRelativeToInterpCell );
            }

            moris::uint
            get_num_cell_clusters() const
            {
                return mInterpolationCells.size();
            }

            bool
            is_trivial( moris_index aClusterIndex ) const
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );

                moris::uint tTrivFlag = mTrivialFlag( aClusterIndex );

                if ( tTrivFlag > 0 )
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            mtk::Cell const *
            get_interp_cell( moris_index aClusterIndex ) const
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(),
                        " Cell cluster index out of bounds" );

                return mInterpolationCells( aClusterIndex );
            }

            moris::Matrix< IndexMat > const *
            get_integration_cell_ids_and_side_ords( moris_index aClusterIndex ) const
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mCellIdsAndSideOrdinals( aClusterIndex );
            }

            moris::Matrix< IndexMat > const *
            get_vertex_in_cluster_ids( moris_index aClusterIndex ) const
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mVerticesIdsInCluster( aClusterIndex );
            }

            moris::Matrix< DDRMat > const *
            get_vertex_local_coords_wrt_interpolation_cell( moris_index aClusterIndex ) const
            {
                MORIS_ASSERT( aClusterIndex < (moris_index)get_num_cell_clusters(), " Cell cluster index out of bounds" );
                return mLocalCoordsRelativeToInterpCell( aClusterIndex );
            }

          private:
            moris::Cell< moris::uint >                       mTrivialFlag;
            moris::Cell< mtk::Cell * >                       mInterpolationCells;
            moris::Cell< moris::Matrix< IndexMat > const * > mCellIdsAndSideOrdinals;
            moris::Cell< moris::Matrix< IndexMat > const * > mVerticesIdsInCluster;
            moris::Cell< moris::Matrix< DDRMat > const * >   mLocalCoordsRelativeToInterpCell;
            std::unordered_map< moris_index, moris_index >   mInterpCellIndexToIndex;
        };

        class Side_Cluster_Input
        {
            moris::Cell< Side_Set_Cluster_Data > mSideClusters; /*all side cluster data for a given side set*/

            std::unordered_map< std::string, moris_index > mSideSetLabelToOrd;

            moris::Cell< std::string > mSideSetLabels;

          public:
            Side_Cluster_Input(){};

            /*!
             * Registers a side set label, also returns the ordinal of this side set
             */
            moris_index
            add_side_set_label( std::string aSideLabel )
            {
                MORIS_ASSERT( mSideSetLabelToOrd.find( aSideLabel ) == mSideSetLabelToOrd.end(),
                        "Trying to add a side set label which has already been registered" );

                moris::moris_index tIndex        = mSideClusters.size();
                mSideSetLabelToOrd[ aSideLabel ] = tIndex;

                mSideClusters.resize( tIndex + 1 );

                mSideSetLabels.push_back( aSideLabel );

                return tIndex;
            }

            /*!
             * Returns the side label ordinal in this input data for a given side label. Retunrs MORIS_INDEX_MAX if not in list
             */
            moris_index
            get_side_label_ordinal( std::string aSideLabel )
            {
                auto tIter = mSideSetLabelToOrd.find( aSideLabel );

                if ( tIter == mSideSetLabelToOrd.end() )
                {
                    return MORIS_INDEX_MAX;
                }
                else
                {
                    return tIter->second;
                }
            }

            void
            add_cluster_data(
                    bool                             aTrivial,
                    moris::uint                      aSideSetOrd,
                    mtk::Cell                       *aInterpCell,
                    moris::Matrix< IndexMat > const *aCellIdsAndSideOrds,
                    moris::Matrix< IndexMat > const *aVerticesInCluster,
                    moris::Matrix< DDRMat > const   *aLocalCoordsRelativeToInterpCell )
            {
                MORIS_ASSERT( aVerticesInCluster->numel() == aLocalCoordsRelativeToInterpCell->n_rows(),
                        "Number of vertices in the cluster must match the number of rows in local coord mat" );

                MORIS_ASSERT( aSideSetOrd < mSideSetLabelToOrd.size(),
                        "Side set ordinal out of bounds." );

                mSideClusters( aSideSetOrd ).add_cluster_data(    //
                        aTrivial,
                        aInterpCell,
                        aCellIdsAndSideOrds,
                        aVerticesInCluster,
                        aLocalCoordsRelativeToInterpCell );
            }

            Side_Set_Cluster_Data const &
            get_cluster_data( moris::uint aSideSetOrd ) const
            {
                MORIS_ASSERT( aSideSetOrd < mSideSetLabelToOrd.size(), "Side set ordinal out of bounds." );

                return mSideClusters( aSideSetOrd );
            }

            void
            print()
            {
                std::cout << std::left << std::setw( 10 ) << "Number of Side Sets with Side Cluster Information:" << mSideClusters.size() << "\n";

                for ( moris::uint i = 0; i < mSideClusters.size(); i++ )
                {
                    std::cout << "Vertex Set Name:" << mSideSetLabels( i ) << " | Number of Interp Cells: " << mSideClusters( i ).get_num_cell_clusters() << std::endl;
                }
            }
        };

    }    // namespace mtk
}    // namespace moris
#endif /* PROJECTS_MTK_SRC_CL_MTK_SIDE_CLUSTER_INPUT_HPP_ */
