/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Double_Side_Cluster_Input.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_INPUT_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_INPUT_HPP_

#include "cl_MTK_Side_Cluster_Input.hpp"

namespace moris
{
    namespace mtk
    {

        class Double_Side_Cluster_Input
        {

          public:
            moris::Cell< Side_Set_Cluster_Data >                         mLeftSideClusters;  /*all side cluster data for a given side set*/
            moris::Cell< Side_Set_Cluster_Data >                         mRightSideClusters; /*all side cluster data for a given side set*/
            moris::Cell< moris::Cell< moris::Matrix< moris::IdMat >* > > mVertexPairing;     /*vertex pairing for a given cluster*/

            moris::Cell< std::string >                     mDoubleSideSetLabel;
            std::unordered_map< std::string, moris_index > mSideSetLabelToOrd;

            moris_index
            add_double_side_set_label( std::string aSideLabel )
            {
                MORIS_ASSERT( mSideSetLabelToOrd.find( aSideLabel ) == mSideSetLabelToOrd.end(), "Trying to add a side set label which has already been registered" );

                moris::moris_index tIndex = mDoubleSideSetLabel.size();
                mDoubleSideSetLabel.push_back( aSideLabel );
                mSideSetLabelToOrd[ aSideLabel ] = tIndex;

                mLeftSideClusters.resize( tIndex + 1 );
                mRightSideClusters.resize( tIndex + 1 );
                mVertexPairing.resize( tIndex + 1 );

                return tIndex;
            }

            /*!
             * Returns the side label ordinal in this input data for a given side label. Retunrs MORIS_INDEX_MAX if not in list
             */
            moris_index
            get_double_side_label_ordinal( std::string aSideLabel )
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

            moris::Cell< std::string > const &
            get_double_side_set_labels()
            {
                return mDoubleSideSetLabel;
            }

            /*!
             * add left and right side to input for stk mesh. Note: this input list is long but is necessary to ensure
             * we add left and right at the same time.
             */
            void
            add_cluster_data( moris::uint      aSideSetOrd,
                    bool                       aLeftTrivial,
                    mtk::Cell*                 aLeftInterpCell,
                    moris::Matrix< IndexMat >* aLeftCellIdsAndSideOrds,
                    moris::Matrix< IndexMat >* aLeftVerticesInCluster,
                    moris::Matrix< DDRMat >*   aLeftLocalCoordsRelativeToInterpCell,
                    bool                       aRightTrivial,
                    mtk::Cell*                 aRightInterpCell,
                    moris::Matrix< IndexMat >* aRightCellIdsAndSideOrds,
                    moris::Matrix< IndexMat >* aRightVerticesInCluster,
                    moris::Matrix< DDRMat >*   aRightLocalCoordsRelativeToInterpCell,
                    moris::Matrix< IdMat >*    aVertexPairing )
            {
                MORIS_ASSERT( aLeftVerticesInCluster->numel() == aLeftLocalCoordsRelativeToInterpCell->n_rows(), "Number of vertices in the cluster must match the number of rows in local coord mat" );
                MORIS_ASSERT( aLeftVerticesInCluster->numel() == aLeftLocalCoordsRelativeToInterpCell->n_rows(), "Number of vertices in the cluster must match the number of rows in local coord mat" );

                MORIS_ASSERT( aSideSetOrd < mSideSetLabelToOrd.size(), "Side set ordinal out of bounds." );
                MORIS_ASSERT( aVertexPairing->n_cols() == 2, "There needs to be two columns (one for each vertex pairing) Col 0 - Left vert, Col 1 Right Vert" );

                mLeftSideClusters( aSideSetOrd ).add_cluster_data( aLeftTrivial, aLeftInterpCell, aLeftCellIdsAndSideOrds, aLeftVerticesInCluster, aLeftLocalCoordsRelativeToInterpCell );

                mRightSideClusters( aSideSetOrd ).add_cluster_data( aRightTrivial, aRightInterpCell, aRightCellIdsAndSideOrds, aRightVerticesInCluster, aRightLocalCoordsRelativeToInterpCell );

                mVertexPairing( aSideSetOrd ).push_back( aVertexPairing );
            }
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_DOUBLE_SIDE_CLUSTER_INPUT_HPP_ */

