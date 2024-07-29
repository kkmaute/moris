/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Sets_Info.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_

#include "cl_MTK_Node_Sets_Info.hpp"
#include "cl_MTK_Side_Sets_Info.hpp"
#include "cl_MTK_Block_Sets_Info.hpp"

namespace moris
{

    namespace mtk
    {
        ///////////////////////////////
        // STRUC FOR SETS CONTAINER  //
        ///////////////////////////////

        struct MtkSetsInfo
        {
            Vector< MtkNodeSetInfo* >  NodeSetsInfo;
            Vector< MtkSideSetInfo* >  SideSetsInfo;
            Vector< MtkBlockSetInfo* > BlockSetsInfo;

            MtkSetsInfo()
                    : NodeSetsInfo( 0, nullptr )
                    , SideSetsInfo( 0, nullptr )
                    , BlockSetsInfo( 0, nullptr )
            {
            }

            //------------------------------------------------
            // Add sets to data structure
            //------------------------------------------------

            void
            add_node_set( MtkNodeSetInfo* aNodeSet )
            {
                NodeSetsInfo.push_back( aNodeSet );
            }

            void
            add_side_set( MtkSideSetInfo* aSideSet )
            {
                SideSetsInfo.push_back( aSideSet );
            }

            void
            add_block_set( MtkBlockSetInfo* aBlockSet )
            {
                BlockSetsInfo.push_back( aBlockSet );
            }

            //------------------------------------------------
            // Node set access
            //------------------------------------------------
            uint
            get_num_node_sets() const
            {
                return NodeSetsInfo.size();
            }

            Vector< MtkNodeSetInfo* > const &
            get_node_sets() const
            {
                return NodeSetsInfo;
            }

            MtkNodeSetInfo*
            get_node_set( uint aNodeSetIndex ) const
            {
                return NodeSetsInfo( aNodeSetIndex );
            }

            //------------------------------------------------
            // Side set access
            //------------------------------------------------
            uint
            get_num_side_sets() const
            {
                return SideSetsInfo.size();
            }

            Vector< MtkSideSetInfo* > const &
            get_side_sets() const
            {
                return SideSetsInfo;
            }

            MtkSideSetInfo*
            get_side_set( uint aSideSetIndex ) const
            {
                return SideSetsInfo( aSideSetIndex );
            }

            //------------------------------------------------
            // Block set access
            //------------------------------------------------
            uint
            get_num_block_sets() const
            {
                return BlockSetsInfo.size();
            }

            Vector< MtkBlockSetInfo* > const &
            get_block_sets() const
            {
                return BlockSetsInfo;
            }

            MtkBlockSetInfo*
            get_block_set( uint aBlockSetIndex ) const
            {
                return BlockSetsInfo( aBlockSetIndex );
            }

            void
            print()
            {
                moris::print( NodeSetsInfo, "Vertex set information" );
                moris::print( SideSetsInfo, "Side set information" );
                moris::print( BlockSetsInfo, "Block set Information" );
            }
        };
    }    // namespace mtk
}    // namespace moris

#endif /* PROJECTS_MTK_SRC_CL_MTK_SETS_INFO_HPP_ */
