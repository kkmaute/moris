/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Mapper_Node.hpp
 *
 */

#ifndef PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_NODE_HPP_
#define PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_NODE_HPP_

#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "fn_norm.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_Cell.hpp"

namespace moris
{
    namespace mapper
    {
//------------------------------------------------------------------------------

        class Node
        {
            // ref to vertex on MTK mesh
            const mtk::Vertex * mVertex;

            //! multi purpose flag
            bool mFlag;

            uint mNeighborCounter = 0;
            moris::Cell< Node * > mNeighbors;

            real mDistance;
            Matrix< DDRMat >   mCoords;
            Matrix< DDRMat >   mWeights;
            Matrix< IndexMat > mNodeIndices;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Node( const mtk::Vertex * aVertex ) :
                mVertex( aVertex ),
                mCoords( aVertex->get_coords() )
            {
            }

//------------------------------------------------------------------------------

            ~Node()
            {
                mNeighbors.clear();
            };

//------------------------------------------------------------------------------

            void
            flag()
            {
                mFlag = true;
            }

//------------------------------------------------------------------------------

            void
            unflag()
            {
                mFlag = false;
            }

//------------------------------------------------------------------------------

            bool
            is_flagged() const
            {
                return mFlag;
            }

//------------------------------------------------------------------------------

            /**
             * Return level of node. Returns zero if not HMR
             */
            uint
            get_level() const
            {
                return mVertex->get_level();
            }

//------------------------------------------------------------------------------

            moris_index
            get_index() const
            {
                return mVertex->get_index();
            }

//------------------------------------------------------------------------------

            moris_id
            get_id() const
            {
                return mVertex->get_id();
            }

//------------------------------------------------------------------------------

            void
            init_neighbor_container( const uint aNumberOfNeighbors )
            {
                mNeighbors.resize( aNumberOfNeighbors, nullptr );
            }

//------------------------------------------------------------------------------

            void
            insert_neighbor( Node * aNeighbor )
            {
                mNeighbors( mNeighborCounter++ ) = aNeighbor;
            }

//------------------------------------------------------------------------------

            uint
            get_number_of_neighbors() const
            {
                return mNeighborCounter;
            }

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_coords()
            {
                return mCoords;
            }

//------------------------------------------------------------------------------

            Matrix< DDRMat > &
            get_weights()
            {
                return mWeights;
            }

//------------------------------------------------------------------------------

            const Matrix< DDRMat > &
            get_weights() const
            {
                return mWeights;
            }

//------------------------------------------------------------------------------

            Matrix< IndexMat >  &
            get_node_indices()
            {
                return mNodeIndices;
            }

//------------------------------------------------------------------------------

            const Matrix< IndexMat >  &
            get_node_indices() const
            {
                return mNodeIndices;
            }

//------------------------------------------------------------------------------

            real
            get_distance( const Matrix< DDRMat > & aCoords )
            {
                mDistance = norm ( aCoords - mCoords );
                return mDistance;
            }

//------------------------------------------------------------------------------

            real
            get_distance()
            {
                return mDistance;
            }

//------------------------------------------------------------------------------

            void
            get_nodes_in_proximity(
                    const Matrix< DDRMat > & aCoords,
                    const             real & aDistance,
                            Cell< Node * > & aNodes )
            {
                for( Node * tNeighbor : mNeighbors )
                {
                    // test if neighbor is flagged
                    if( ! tNeighbor->is_flagged() )
                    {
                        if( tNeighbor->get_distance( aCoords ) <= aDistance )
                        {
                            tNeighbor->flag();
                            aNodes.push_back( tNeighbor );
                            tNeighbor->get_nodes_in_proximity(
                                    aCoords,
                                    aDistance,
                                    aNodes );
                        }
                    }
                }
            }

        };

//------------------------------------------------------------------------------
    }
}

#endif /* PROJECTS_MTK_MAP_SRC_CL_MTK_MAPPER_NODE_HPP_ */

