/*
 * fn_HMR_Background_Element_Edges_3D.hpp
 *
 *  Created on: September 27, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGED_3D_HPP_
#define SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGES_3D_HPP_

#include "typedefs.hpp"
#include "../../../HMR/src/cl_HMR_Background_Element_Base.hpp"
#include "../../../HMR/src/cl_HMR_Background_Element.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        template <>
        void
        Background_Element< 3, 8, 26, 6, 12 >::create_edges()
        {
            // step 1: copy existing edges from neighbors 0 to 5

            // test if neighbor 0 exists and is on same level
            if( mNeighbors[ 0 ] != NULL )
            {
                if( mNeighbors[ 0 ]->get_level() == mLevel )
                {
                    // test edge 0 deos not exists
                    if( mEdges[ 0 ] == NULL )
                    {
                       // copy edge 2 from neighbor
                       this->insert_edge( mNeighbors[ 0 ]->get_edge( 2 ), 0 );
                    }
                    // test edge 0 deos not exists
                    if( mEdges[ 4 ] == NULL )
                    {
                       // copy edge 7 from neighbor
                       this->insert_edge( mNeighbors[ 0 ]->get_edge( 7 ), 4 );
                    }
                    // test edge 0 deos not exists
                    if( mEdges[ 5 ] == NULL )
                    {
                       // copy edge 6 from neighbor
                       this->insert_edge( mNeighbors[ 0 ]->get_edge( 6 ), 5 );
                    }
                    // test edge 0 deos not exists
                    if( mEdges[ 8 ] == NULL )
                    {
                       // copy edge 10 from neighbor
                       this->insert_edge( mNeighbors[ 0 ]->get_edge( 10 ), 8 );
                    }
                }

            }

            // test if neighbor 1 exists and is on same level
            if( mNeighbors[ 1 ] != NULL )
            {
                if( mNeighbors[ 1 ]->get_level() == mLevel )
                {
                    // test edge 1 deos not exists
                    if( mEdges[ 1 ] == NULL )
                    {
                       // copy edge 3 from neighbor
                       this->insert_edge( mNeighbors[ 1 ]->get_edge( 3 ), 1 );
                    }
                    // test edge 1 deos not exists
                    if( mEdges[ 5 ] == NULL )
                    {
                       // copy edge 4 from neighbor
                       this->insert_edge( mNeighbors[ 1 ]->get_edge( 4 ), 5 );
                    }
                    // test edge 1 deos not exists
                    if( mEdges[ 6 ] == NULL )
                    {
                       // copy edge 7 from neighbor
                       this->insert_edge( mNeighbors[ 1 ]->get_edge( 7 ), 6 );
                    }
                    // test edge 1 deos not exists
                    if( mEdges[ 9 ] == NULL )
                    {
                       // copy edge 11 from neighbor
                       this->insert_edge( mNeighbors[ 1 ]->get_edge( 11 ), 9 );
                    }
                }

            }

            // test if neighbor 2 exists and is on same level
            if( mNeighbors[ 2 ] != NULL )
            {
                if( mNeighbors[ 2 ]->get_level() == mLevel )
                {
                    // test edge 2 deos not exists
                    if( mEdges[ 2 ] == NULL )
                    {
                       // copy edge 0 from neighbor
                       this->insert_edge( mNeighbors[ 2 ]->get_edge( 0 ), 2 );
                    }
                    // test edge 2 deos not exists
                    if( mEdges[ 6 ] == NULL )
                    {
                       // copy edge 5 from neighbor
                       this->insert_edge( mNeighbors[ 2 ]->get_edge( 5 ), 6 );
                    }
                    // test edge 2 deos not exists
                    if( mEdges[ 7 ] == NULL )
                    {
                       // copy edge 4 from neighbor
                       this->insert_edge( mNeighbors[ 2 ]->get_edge( 4 ), 7 );
                    }
                    // test edge 2 deos not exists
                    if( mEdges[ 10 ] == NULL )
                    {
                       // copy edge 8 from neighbor
                       this->insert_edge( mNeighbors[ 2 ]->get_edge( 8 ), 10 );
                    }
                }

            }

            // test if neighbor 3 exists and is on same level
            if( mNeighbors[ 3 ] != NULL )
            {
                if( mNeighbors[ 3 ]->get_level() == mLevel )
                {
                    // test edge 3 deos not exists
                    if( mEdges[ 3 ] == NULL )
                    {
                       // copy edge 1 from neighbor
                       this->insert_edge( mNeighbors[ 3 ]->get_edge( 1 ), 3 );
                    }
                    // test edge 3 deos not exists
                    if( mEdges[ 4 ] == NULL )
                    {
                       // copy edge 5 from neighbor
                       this->insert_edge( mNeighbors[ 3 ]->get_edge( 5 ), 4 );
                    }
                    // test edge 3 deos not exists
                    if( mEdges[ 7 ] == NULL )
                    {
                       // copy edge 6 from neighbor
                       this->insert_edge( mNeighbors[ 3 ]->get_edge( 6 ), 7 );
                    }
                    // test edge 3 deos not exists
                    if( mEdges[ 11 ] == NULL )
                    {
                       // copy edge 9 from neighbor
                       this->insert_edge( mNeighbors[ 3 ]->get_edge( 9 ), 11 );
                    }
                }

            }

            // test if neighbor 4 exists and is on same level
            if( mNeighbors[ 4 ] != NULL )
            {
                if( mNeighbors[ 4 ]->get_level() == mLevel )
                {
                    // test edge 4 deos not exists
                    if( mEdges[ 0 ] == NULL )
                    {
                       // copy edge 8 from neighbor
                       this->insert_edge( mNeighbors[ 4 ]->get_edge( 8 ), 0 );
                    }
                    // test edge 4 deos not exists
                    if( mEdges[ 1 ] == NULL )
                    {
                       // copy edge 9 from neighbor
                       this->insert_edge( mNeighbors[ 4 ]->get_edge( 9 ), 1 );
                    }
                    // test edge 4 deos not exists
                    if( mEdges[ 2 ] == NULL )
                    {
                       // copy edge 10 from neighbor
                       this->insert_edge( mNeighbors[ 4 ]->get_edge( 10 ), 2 );
                    }
                    // test edge 4 deos not exists
                    if( mEdges[ 3 ] == NULL )
                    {
                       // copy edge 11 from neighbor
                       this->insert_edge( mNeighbors[ 4 ]->get_edge( 11 ), 3 );
                    }
                }

            }

            // test if neighbor 5 exists and is on same level
            if( mNeighbors[ 5 ] != NULL )
            {
                if( mNeighbors[ 5 ]->get_level() == mLevel )
                {
                    // test edge 5 deos not exists
                    if( mEdges[ 8 ] == NULL )
                    {
                       // copy edge 0 from neighbor
                       this->insert_edge( mNeighbors[ 5 ]->get_edge( 0 ), 8 );
                    }
                    // test edge 5 deos not exists
                    if( mEdges[ 9 ] == NULL )
                    {
                       // copy edge 1 from neighbor
                       this->insert_edge( mNeighbors[ 5 ]->get_edge( 1 ), 9 );
                    }
                    // test edge 5 deos not exists
                    if( mEdges[ 10 ] == NULL )
                    {
                       // copy edge 2 from neighbor
                       this->insert_edge( mNeighbors[ 5 ]->get_edge( 2 ), 10 );
                    }
                    // test edge 5 deos not exists
                    if( mEdges[ 11 ] == NULL )
                    {
                       // copy edge 3 from neighbor
                       this->insert_edge( mNeighbors[ 5 ]->get_edge( 3 ), 11 );
                    }
                }

            }

            // step 2: copy existing edges from neighbors 6 to 17

            // test if edge 0 does not exist
            if( mEdges[ 0 ] == NULL )
            {
                // test if neighbor 25 exists and is on same level
                if( mNeighbors[ 6 ] != NULL )
                {
                    if( mNeighbors[ 6 ]->get_level() == mLevel )
                    {
                       // copy edge 10 from neighbor
                       this->insert_edge( mNeighbors[ 6 ]->get_edge( 10 ), 0 );
                    }
                }
            }

            // test if edge 1 does not exist
            if( mEdges[ 1 ] == NULL )
            {
                // test if neighbor 26 exists and is on same level
                if( mNeighbors[ 7 ] != NULL )
                {
                    if( mNeighbors[ 7 ]->get_level() == mLevel )
                    {
                       // copy edge 11 from neighbor
                       this->insert_edge( mNeighbors[ 7 ]->get_edge( 11 ), 1 );
                    }
                }
            }

            // test if edge 2 does not exist
            if( mEdges[ 2 ] == NULL )
            {
                // test if neighbor 27 exists and is on same level
                if( mNeighbors[ 8 ] != NULL )
                {
                    if( mNeighbors[ 8 ]->get_level() == mLevel )
                    {
                       // copy edge 8 from neighbor
                       this->insert_edge( mNeighbors[ 8 ]->get_edge( 8 ), 2 );
                    }
                }
            }

            // test if edge 3 does not exist
            if( mEdges[ 3 ] == NULL )
            {
                // test if neighbor 28 exists and is on same level
                if( mNeighbors[ 9 ] != NULL )
                {
                    if( mNeighbors[ 9 ]->get_level() == mLevel )
                    {
                       // copy edge 9 from neighbor
                       this->insert_edge( mNeighbors[ 9 ]->get_edge( 9 ), 3 );
                    }
                }
            }

            // test if edge 4 does not exist
            if( mEdges[ 4 ] == NULL )
            {
                // test if neighbor 29 exists and is on same level
                if( mNeighbors[ 10 ] != NULL )
                {
                    if( mNeighbors[ 10 ]->get_level() == mLevel )
                    {
                       // copy edge 6 from neighbor
                       this->insert_edge( mNeighbors[ 10 ]->get_edge( 6 ), 4 );
                    }
                }
            }

            // test if edge 5 does not exist
            if( mEdges[ 5 ] == NULL )
            {
                // test if neighbor 30 exists and is on same level
                if( mNeighbors[ 11 ] != NULL )
                {
                    if( mNeighbors[ 11 ]->get_level() == mLevel )
                    {
                       // copy edge 7 from neighbor
                       this->insert_edge( mNeighbors[ 11 ]->get_edge( 7 ), 5 );
                    }
                }
            }

            // test if edge 6 does not exist
            if( mEdges[ 6 ] == NULL )
            {
                // test if neighbor 31 exists and is on same level
                if( mNeighbors[ 12 ] != NULL )
                {
                    if( mNeighbors[ 12 ]->get_level() == mLevel )
                    {
                       // copy edge 4 from neighbor
                       this->insert_edge( mNeighbors[ 12 ]->get_edge( 4 ), 6 );
                    }
                }
            }

            // test if edge 7 does not exist
            if( mEdges[ 7 ] == NULL )
            {
                // test if neighbor 32 exists and is on same level
                if( mNeighbors[ 13 ] != NULL )
                {
                    if( mNeighbors[ 13 ]->get_level() == mLevel )
                    {
                       // copy edge 5 from neighbor
                       this->insert_edge( mNeighbors[ 13 ]->get_edge( 5 ), 7 );
                    }
                }
            }

            // test if edge 8 does not exist
            if( mEdges[ 8 ] == NULL )
            {
                // test if neighbor 33 exists and is on same level
                if( mNeighbors[ 14 ] != NULL )
                {
                    if( mNeighbors[ 14 ]->get_level() == mLevel )
                    {
                       // copy edge 2 from neighbor
                       this->insert_edge( mNeighbors[ 14 ]->get_edge( 2 ), 8 );
                    }
                }
            }

            // test if edge 9 does not exist
            if( mEdges[ 9 ] == NULL )
            {
                // test if neighbor 34 exists and is on same level
                if( mNeighbors[ 15 ] != NULL )
                {
                    if( mNeighbors[ 15 ]->get_level() == mLevel )
                    {
                       // copy edge 3 from neighbor
                       this->insert_edge( mNeighbors[ 15 ]->get_edge( 3 ), 9 );
                    }
                }
            }

            // test if edge 10 does not exist
            if( mEdges[ 10 ] == NULL )
            {
                // test if neighbor 35 exists and is on same level
                if( mNeighbors[ 16 ] != NULL )
                {
                    if( mNeighbors[ 16 ]->get_level() == mLevel )
                    {
                       // copy edge 0 from neighbor
                       this->insert_edge( mNeighbors[ 16 ]->get_edge( 0 ), 10 );
                    }
                }
            }

            // test if edge 11 does not exist
            if( mEdges[ 11 ] == NULL )
            {
                // test if neighbor 36 exists and is on same level
                if( mNeighbors[ 17 ] != NULL )
                {
                    if( mNeighbors[ 17 ]->get_level() == mLevel )
                    {
                       // copy edge 1 from neighbor
                       this->insert_edge( mNeighbors[ 17 ]->get_edge( 1 ), 11 );
                    }
                }
            }

            // step 3: create edges that do not exist

            // test if edge 0 does not exist
            if( mEdges[ 0 ] == NULL )
            {
                // create edge
                mEdges[ 0 ] = new Background_Edge( this, 0 );

                // set owning flag
                mEdgeOwnFlags.set( 0 );

                // test if neighbor 0 exists
                if( mNeighbors[ 0 ] != NULL )
                {
                    if( mNeighbors[ 0 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 0 ]->insert_edge( mEdges[ 0 ], 2 );
                    }
                }

                // test if neighbor 4 exists
                if( mNeighbors[ 4 ] != NULL )
                {
                    if( mNeighbors[ 4 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 4 ]->insert_edge( mEdges[ 0 ], 8 );
                    }
                }

                // test if neighbor 6 exists
                if( mNeighbors[ 6 ] != NULL )
                {
                    if( mNeighbors[ 6 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 6 ]->insert_edge( mEdges[ 0 ], 10 );
                    }
                }
            }
            // test if edge 1 does not exist
            if( mEdges[ 1 ] == NULL )
            {
                // create edge
                mEdges[ 1 ] = new Background_Edge( this, 1 );

                // set owning flag
                mEdgeOwnFlags.set( 1 );

                // test if neighbor 1 exists
                if( mNeighbors[ 1 ] != NULL )
                {
                    if( mNeighbors[ 1 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 1 ]->insert_edge( mEdges[ 1 ], 3 );
                    }
                }

                // test if neighbor 4 exists
                if( mNeighbors[ 4 ] != NULL )
                {
                    if( mNeighbors[ 4 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 4 ]->insert_edge( mEdges[ 1 ], 9 );
                    }
                }

                // test if neighbor 7 exists
                if( mNeighbors[ 7 ] != NULL )
                {
                    if( mNeighbors[ 7 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 7 ]->insert_edge( mEdges[ 1 ], 11 );
                    }
                }
            }
            // test if edge 2 does not exist
            if( mEdges[ 2 ] == NULL )
            {
                // create edge
                mEdges[ 2 ] = new Background_Edge( this, 2 );

                // set owning flag
                mEdgeOwnFlags.set( 2 );

                // test if neighbor 2 exists
                if( mNeighbors[ 2 ] != NULL )
                {
                    if( mNeighbors[ 2 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 2 ]->insert_edge( mEdges[ 2 ], 0 );
                    }
                }

                // test if neighbor 4 exists
                if( mNeighbors[ 4 ] != NULL )
                {
                    if( mNeighbors[ 4 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 4 ]->insert_edge( mEdges[ 2 ], 10 );
                    }
                }

                // test if neighbor 8 exists
                if( mNeighbors[ 8 ] != NULL )
                {
                    if( mNeighbors[ 8 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 8 ]->insert_edge( mEdges[ 2 ], 8 );
                    }
                }
            }
            // test if edge 3 does not exist
            if( mEdges[ 3 ] == NULL )
            {
                // create edge
                mEdges[ 3 ] = new Background_Edge( this, 3 );

                // set owning flag
                mEdgeOwnFlags.set( 3 );

                // test if neighbor 3 exists
                if( mNeighbors[ 3 ] != NULL )
                {
                    if( mNeighbors[ 3 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 3 ]->insert_edge( mEdges[ 3 ], 1 );
                    }
                }

                // test if neighbor 4 exists
                if( mNeighbors[ 4 ] != NULL )
                {
                    if( mNeighbors[ 4 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 4 ]->insert_edge( mEdges[ 3 ], 11 );
                    }
                }

                // test if neighbor 9 exists
                if( mNeighbors[ 9 ] != NULL )
                {
                    if( mNeighbors[ 9 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 9 ]->insert_edge( mEdges[ 3 ], 9 );
                    }
                }
            }
            // test if edge 4 does not exist
            if( mEdges[ 4 ] == NULL )
            {
                // create edge
                mEdges[ 4 ] = new Background_Edge( this, 4 );

                // set owning flag
                mEdgeOwnFlags.set( 4 );

                // test if neighbor 0 exists
                if( mNeighbors[ 0 ] != NULL )
                {
                    if( mNeighbors[ 0 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 0 ]->insert_edge( mEdges[ 4 ], 7 );
                    }
                }

                // test if neighbor 3 exists
                if( mNeighbors[ 3 ] != NULL )
                {
                    if( mNeighbors[ 3 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 3 ]->insert_edge( mEdges[ 4 ], 5 );
                    }
                }

                // test if neighbor 10 exists
                if( mNeighbors[ 10 ] != NULL )
                {
                    if( mNeighbors[ 10 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 10 ]->insert_edge( mEdges[ 4 ], 6 );
                    }
                }
            }
            // test if edge 5 does not exist
            if( mEdges[ 5 ] == NULL )
            {
                // create edge
                mEdges[ 5 ] = new Background_Edge( this, 5 );

                // set owning flag
                mEdgeOwnFlags.set( 5 );

                // test if neighbor 0 exists
                if( mNeighbors[ 0 ] != NULL )
                {
                    if( mNeighbors[ 0 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 0 ]->insert_edge( mEdges[ 5 ], 6 );
                    }
                }

                // test if neighbor 1 exists
                if( mNeighbors[ 1 ] != NULL )
                {
                    if( mNeighbors[ 1 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 1 ]->insert_edge( mEdges[ 5 ], 4 );
                    }
                }

                // test if neighbor 11 exists
                if( mNeighbors[ 11 ] != NULL )
                {
                    if( mNeighbors[ 11 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 11 ]->insert_edge( mEdges[ 5 ], 7 );
                    }
                }
            }
            // test if edge 6 does not exist
            if( mEdges[ 6 ] == NULL )
            {
                // create edge
                mEdges[ 6 ] = new Background_Edge( this, 6 );

                // set owning flag
                mEdgeOwnFlags.set( 6 );

                // test if neighbor 1 exists
                if( mNeighbors[ 1 ] != NULL )
                {
                    if( mNeighbors[ 1 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 1 ]->insert_edge( mEdges[ 6 ], 7 );
                    }
                }

                // test if neighbor 2 exists
                if( mNeighbors[ 2 ] != NULL )
                {
                    if( mNeighbors[ 2 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 2 ]->insert_edge( mEdges[ 6 ], 5 );
                    }
                }

                // test if neighbor 12 exists
                if( mNeighbors[ 12 ] != NULL )
                {
                    if( mNeighbors[ 12 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 12 ]->insert_edge( mEdges[ 6 ], 4 );
                    }
                }
            }
            // test if edge 7 does not exist
            if( mEdges[ 7 ] == NULL )
            {
                // create edge
                mEdges[ 7 ] = new Background_Edge( this, 7 );

                // set owning flag
                mEdgeOwnFlags.set( 7 );

                // test if neighbor 2 exists
                if( mNeighbors[ 2 ] != NULL )
                {
                    if( mNeighbors[ 2 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 2 ]->insert_edge( mEdges[ 7 ], 4 );
                    }
                }

                // test if neighbor 3 exists
                if( mNeighbors[ 3 ] != NULL )
                {
                    if( mNeighbors[ 3 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 3 ]->insert_edge( mEdges[ 7 ], 6 );
                    }
                }

                // test if neighbor 13 exists
                if( mNeighbors[ 13 ] != NULL )
                {
                    if( mNeighbors[ 13 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 13 ]->insert_edge( mEdges[ 7 ], 5 );
                    }
                }
            }
            // test if edge 8 does not exist
            if( mEdges[ 8 ] == NULL )
            {
                // create edge
                mEdges[ 8 ] = new Background_Edge( this, 8 );

                // set owning flag
                mEdgeOwnFlags.set( 8 );

                // test if neighbor 0 exists
                if( mNeighbors[ 0 ] != NULL )
                {
                    if( mNeighbors[ 0 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 0 ]->insert_edge( mEdges[ 8 ], 10 );
                    }
                }

                // test if neighbor 5 exists
                if( mNeighbors[ 5 ] != NULL )
                {
                    if( mNeighbors[ 5 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 5 ]->insert_edge( mEdges[ 8 ], 0 );
                    }
                }

                // test if neighbor 14 exists
                if( mNeighbors[ 14 ] != NULL )
                {
                    if( mNeighbors[ 14 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 14 ]->insert_edge( mEdges[ 8 ], 2 );
                    }
                }
            }
            // test if edge 9 does not exist
            if( mEdges[ 9 ] == NULL )
            {
                // create edge
                mEdges[ 9 ] = new Background_Edge( this, 9 );

                // set owning flag
                mEdgeOwnFlags.set( 9 );

                // test if neighbor 1 exists
                if( mNeighbors[ 1 ] != NULL )
                {
                    if( mNeighbors[ 1 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 1 ]->insert_edge( mEdges[ 9 ], 11 );
                    }
                }

                // test if neighbor 5 exists
                if( mNeighbors[ 5 ] != NULL )
                {
                    if( mNeighbors[ 5 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 5 ]->insert_edge( mEdges[ 9 ], 1 );
                    }
                }

                // test if neighbor 15 exists
                if( mNeighbors[ 15 ] != NULL )
                {
                    if( mNeighbors[ 15 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 15 ]->insert_edge( mEdges[ 9 ], 3 );
                    }
                }
            }
            // test if edge 10 does not exist
            if( mEdges[ 10 ] == NULL )
            {
                // create edge
                mEdges[ 10 ] = new Background_Edge( this, 10 );

                // set owning flag
                mEdgeOwnFlags.set( 10 );

                // test if neighbor 2 exists
                if( mNeighbors[ 2 ] != NULL )
                {
                    if( mNeighbors[ 2 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 2 ]->insert_edge( mEdges[ 10 ], 8 );
                    }
                }

                // test if neighbor 5 exists
                if( mNeighbors[ 5 ] != NULL )
                {
                    if( mNeighbors[ 5 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 5 ]->insert_edge( mEdges[ 10 ], 2 );
                    }
                }

                // test if neighbor 16 exists
                if( mNeighbors[ 16 ] != NULL )
                {
                    if( mNeighbors[ 16 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 16 ]->insert_edge( mEdges[ 10 ], 0 );
                    }
                }
            }
            // test if edge 11 does not exist
            if( mEdges[ 11 ] == NULL )
            {
                // create edge
                mEdges[ 11 ] = new Background_Edge( this, 11 );

                // set owning flag
                mEdgeOwnFlags.set( 11 );

                // test if neighbor 3 exists
                if( mNeighbors[ 3 ] != NULL )
                {
                    if( mNeighbors[ 3 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 3 ]->insert_edge( mEdges[ 11 ], 9 );
                    }
                }

                // test if neighbor 5 exists
                if( mNeighbors[ 5 ] != NULL )
                {
                    if( mNeighbors[ 5 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 5 ]->insert_edge( mEdges[ 11 ], 3 );
                    }
                }

                // test if neighbor 17 exists
                if( mNeighbors[ 17 ] != NULL )
                {
                    if( mNeighbors[ 17 ]->get_level() == mLevel )
                    {
                        // copy edge to neighbor
                        mNeighbors[ 17 ]->insert_edge( mEdges[ 11 ], 1 );
                    }
                }
            }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGES_3D_HPP_ */
