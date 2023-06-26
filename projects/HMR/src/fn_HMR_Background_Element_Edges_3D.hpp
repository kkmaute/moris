/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Background_Element_Edges_3D.hpp
 *
 */

#ifndef SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGES_3D_HPP_
#define SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGES_3D_HPP_

#include "typedefs.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "cl_HMR_Background_Element.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------

    template <>
    inline
    void Background_Element< 3 >::create_edges()
    {
        // step 1: copy existing edges from neighbors 0 to tE5
        moris_index tE0  = 0;
        moris_index tE1  = 1;
        moris_index tE2  = 2;
        moris_index tE3  = 3;
        moris_index tE4  = 8;
        moris_index tE5  = 9;
        moris_index tE6  = 10;
        moris_index tE7  = 11;
        moris_index tE8  = 4;
        moris_index tE9  = 5;
        moris_index tE10 = 6;
        moris_index tE11 = 7;

        // test if neighbor 0 exists and is on same level
        if( mNeighbors[ 0 ] != nullptr )
        {
            if( mNeighbors[ 0 ]->get_level() == mLevel )
            {
                // test edge 0 deos not exists
                if( mEdges[ tE0 ] == nullptr )
                {
                   // copy edge 2 from neighbor
                   this->insert_edge( mNeighbors[ 0 ]->get_edge( tE2 ), tE0 );
                }
                // test edge 0 deos not exists
                if( mEdges[ tE4 ] == nullptr )
                {
                   // copy edge tE7 from neighbor
                   this->insert_edge( mNeighbors[ 0 ]->get_edge( tE7 ), tE4 );
                }
                // test edge 0 deos not exists
                if( mEdges[ tE5 ] == nullptr )
                {
                   // copy edge tE6 from neighbor
                   this->insert_edge( mNeighbors[ 0 ]->get_edge( tE6 ), tE5 );
                }
                // test edge 0 deos not exists
                if( mEdges[ tE8 ] == nullptr )
                {
                   // copy edge tE10 from neighbor
                   this->insert_edge( mNeighbors[ 0 ]->get_edge( tE10 ), tE8 );
                }
            }

        }

        // test if neighbor 1 exists and is on same level
        if( mNeighbors[ 1 ] != nullptr )
        {
            if( mNeighbors[ 1 ]->get_level() == mLevel )
            {
                // test edge 1 deos not exists
                if( mEdges[ tE1 ] == nullptr )
                {
                   // copy edge 3 from neighbor
                   this->insert_edge( mNeighbors[ 1 ]->get_edge( tE3 ), tE1 );
                }
                // test edge 1 deos not exists
                if( mEdges[ tE5 ] == nullptr )
                {
                   // copy edge tE4 from neighbor
                   this->insert_edge( mNeighbors[ 1 ]->get_edge( tE4 ), tE5 );
                }
                // test edge 1 deos not exists
                if( mEdges[ tE6 ] == nullptr )
                {
                   // copy edge tE7 from neighbor
                   this->insert_edge( mNeighbors[ 1 ]->get_edge( tE7 ), tE6 );
                }
                // test edge 1 deos not exists
                if( mEdges[ tE9 ] == nullptr )
                {
                   // copy edge tE11 from neighbor
                   this->insert_edge( mNeighbors[ 1 ]->get_edge( tE11 ), tE9 );
                }
            }

        }

        // test if neighbor 2 exists and is on same level
        if( mNeighbors[ 2 ] != nullptr )
        {
            if( mNeighbors[ 2 ]->get_level() == mLevel )
            {
                // test edge 2 deos not exists
                if( mEdges[ tE2 ] == nullptr )
                {
                   // copy edge 0 from neighbor
                   this->insert_edge( mNeighbors[ 2 ]->get_edge( tE0 ), tE2 );
                }
                // test edge 2 deos not exists
                if( mEdges[ tE6 ] == nullptr )
                {
                   // copy edge tE5 from neighbor
                   this->insert_edge( mNeighbors[ 2 ]->get_edge( tE5 ), tE6 );
                }
                // test edge 2 deos not exists
                if( mEdges[ tE7 ] == nullptr )
                {
                   // copy edge tE4 from neighbor
                   this->insert_edge( mNeighbors[ 2 ]->get_edge( tE4 ), tE7 );
                }
                // test edge 2 deos not exists
                if( mEdges[ tE10 ] == nullptr )
                {
                   // copy edge tE8 from neighbor
                   this->insert_edge( mNeighbors[ 2 ]->get_edge( tE8 ), tE10 );
                }
            }

        }

        // test if neighbor 3 exists and is on same level
        if( mNeighbors[ 3 ] != nullptr )
        {
            if( mNeighbors[ 3 ]->get_level() == mLevel )
            {
                // test edge 3 deos not exists
                if( mEdges[ tE3 ] == nullptr )
                {
                   // copy edge 1 from neighbor
                   this->insert_edge( mNeighbors[ 3 ]->get_edge( tE1 ), tE3 );
                }
                // test edge 3 deos not exists
                if( mEdges[ tE4 ] == nullptr )
                {
                   // copy edge tE5 from neighbor
                   this->insert_edge( mNeighbors[ 3 ]->get_edge( tE5 ), tE4 );
                }
                // test edge 3 deos not exists
                if( mEdges[ tE7 ] == nullptr )
                {
                   // copy edge tE6 from neighbor
                   this->insert_edge( mNeighbors[ 3 ]->get_edge( tE6 ), tE7 );
                }
                // test edge 3 deos not exists
                if( mEdges[ tE11 ] == nullptr )
                {
                   // copy edge tE9 from neighbor
                   this->insert_edge( mNeighbors[ 3 ]->get_edge( tE9 ), tE11 );
                }
            }

        }

        // test if neighbor tE4 exists and is on same level
        if( mNeighbors[ 4 ] != nullptr )
        {
            if( mNeighbors[ 4 ]->get_level() == mLevel )
            {
                // test edge tE4 deos not exists
                if( mEdges[ tE0 ] == nullptr )
                {
                   // copy edge tE8 from neighbor
                   this->insert_edge( mNeighbors[ 4 ]->get_edge( tE8 ), tE0 );
                }
                // test edge tE4 deos not exists
                if( mEdges[ tE1 ] == nullptr )
                {
                   // copy edge tE9 from neighbor
                   this->insert_edge( mNeighbors[ 4 ]->get_edge( tE9 ), tE1 );
                }
                // test edge tE4 deos not exists
                if( mEdges[ tE2 ] == nullptr )
                {
                   // copy edge tE10 from neighbor
                   this->insert_edge( mNeighbors[ 4 ]->get_edge( tE10 ), tE2 );
                }
                // test edge tE4 deos not exists
                if( mEdges[ tE3 ] == nullptr )
                {
                   // copy edge tE11 from neighbor
                   this->insert_edge( mNeighbors[ 4 ]->get_edge( tE11 ), tE3 );
                }
            }

        }

        // test if neighbor tE5 exists and is on same level
        if( mNeighbors[ 5 ] != nullptr )
        {
            if( mNeighbors[ 5 ]->get_level() == mLevel )
            {
                // test edge tE5 deos not exists
                if( mEdges[ tE8 ] == nullptr )
                {
                   // copy edge 0 from neighbor
                   this->insert_edge( mNeighbors[ 5 ]->get_edge( tE0 ), tE8 );
                }
                // test edge tE5 deos not exists
                if( mEdges[ tE9 ] == nullptr )
                {
                   // copy edge 1 from neighbor
                   this->insert_edge( mNeighbors[ 5 ]->get_edge( tE1 ), tE9 );
                }
                // test edge tE5 deos not exists
                if( mEdges[ tE10 ] == nullptr )
                {
                   // copy edge 2 from neighbor
                   this->insert_edge( mNeighbors[ 5 ]->get_edge( tE2 ), tE10 );
                }
                // test edge tE5 deos not exists
                if( mEdges[ tE11 ] == nullptr )
                {
                   // copy edge 3 from neighbor
                   this->insert_edge( mNeighbors[ 5 ]->get_edge( tE3 ), tE11 );
                }
            }

        }

        // step 2: copy existing edges from neighbors tE6 to 17

        // test if edge 0 does not exist
        if( mEdges[ tE0 ] == nullptr )
        {
            // test if neighbor 25 exists and is on same level
            if( mNeighbors[ 6 ] != nullptr )
            {
                if( mNeighbors[ 6 ]->get_level() == mLevel )
                {
                   // copy edge tE10 from neighbor
                   this->insert_edge( mNeighbors[ 6 ]->get_edge( tE10 ), tE0 );
                }
            }
        }

        // test if edge 1 does not exist
        if( mEdges[ tE1 ] == nullptr )
        {
            // test if neighbor 26 exists and is on same level
            if( mNeighbors[ 7 ] != nullptr )
            {
                if( mNeighbors[ 7 ]->get_level() == mLevel )
                {
                   // copy edge tE11 from neighbor
                   this->insert_edge( mNeighbors[ 7 ]->get_edge( tE11 ), tE1 );
                }
            }
        }

        // test if edge 2 does not exist
        if( mEdges[ tE2 ] == nullptr )
        {
            // test if neighbor 27 exists and is on same level
            if( mNeighbors[ 8 ] != nullptr )
            {
                if( mNeighbors[ 8 ]->get_level() == mLevel )
                {
                   // copy edge tE8 from neighbor
                   this->insert_edge( mNeighbors[ 8 ]->get_edge( tE8 ), tE2 );
                }
            }
        }

        // test if edge 3 does not exist
        if( mEdges[ tE3 ] == nullptr )
        {
            // test if neighbor 28 exists and is on same level
            if( mNeighbors[ 9 ] != nullptr )
            {
                if( mNeighbors[ 9 ]->get_level() == mLevel )
                {
                   // copy edge tE9 from neighbor
                   this->insert_edge( mNeighbors[ 9 ]->get_edge( tE9 ), tE3 );
                }
            }
        }

        // test if edge tE4 does not exist
        if( mEdges[ tE4 ] == nullptr )
        {
            // test if neighbor 29 exists and is on same level
            if( mNeighbors[ 10 ] != nullptr )
            {
                if( mNeighbors[ 10 ]->get_level() == mLevel )
                {
                   // copy edge tE6 from neighbor
                   this->insert_edge( mNeighbors[ 10 ]->get_edge( tE6 ), tE4 );
                }
            }
        }

        // test if edge tE5 does not exist
        if( mEdges[ tE5 ] == nullptr )
        {
            // test if neighbor 30 exists and is on same level
            if( mNeighbors[ 11 ] != nullptr )
            {
                if( mNeighbors[ 11 ]->get_level() == mLevel )
                {
                   // copy edge tE7 from neighbor
                   this->insert_edge( mNeighbors[ 11 ]->get_edge( tE7 ), tE5 );
                }
            }
        }

        // test if edge tE6 does not exist
        if( mEdges[ tE6 ] == nullptr )
        {
            // test if neighbor 31 exists and is on same level
            if( mNeighbors[ 12 ] != nullptr )
            {
                if( mNeighbors[ 12 ]->get_level() == mLevel )
                {
                   // copy edge tE4 from neighbor
                   this->insert_edge( mNeighbors[ 12 ]->get_edge( tE4 ), tE6 );
                }
            }
        }

        // test if edge tE7 does not exist
        if( mEdges[ tE7 ] == nullptr )
        {
            // test if neighbor 32 exists and is on same level
            if( mNeighbors[ 13 ] != nullptr )
            {
                if( mNeighbors[ 13 ]->get_level() == mLevel )
                {
                   // copy edge tE5 from neighbor
                   this->insert_edge( mNeighbors[ 13 ]->get_edge( tE5 ), tE7 );
                }
            }
        }

        // test if edge tE8 does not exist
        if( mEdges[ tE8 ] == nullptr )
        {
            // test if neighbor 33 exists and is on same level
            if( mNeighbors[ 14 ] != nullptr )
            {
                if( mNeighbors[ 14 ]->get_level() == mLevel )
                {
                   // copy edge 2 from neighbor
                   this->insert_edge( mNeighbors[ 14 ]->get_edge( tE2 ), tE8 );
                }
            }
        }

        // test if edge tE9 does not exist
        if( mEdges[ tE9 ] == nullptr )
        {
            // test if neighbor 34 exists and is on same level
            if( mNeighbors[ 15 ] != nullptr )
            {
                if( mNeighbors[ 15 ]->get_level() == mLevel )
                {
                   // copy edge 3 from neighbor
                   this->insert_edge( mNeighbors[ 15 ]->get_edge( tE3 ), tE9 );
                }
            }
        }

        // test if edge tE10 does not exist
        if( mEdges[ tE10 ] == nullptr )
        {
            // test if neighbor 35 exists and is on same level
            if( mNeighbors[ 16 ] != nullptr )
            {
                if( mNeighbors[ 16 ]->get_level() == mLevel )
                {
                   // copy edge 0 from neighbor
                   this->insert_edge( mNeighbors[ 16 ]->get_edge( tE0 ), tE10 );
                }
            }
        }

        // test if edge tE11 does not exist
        if( mEdges[ tE11 ] == nullptr )
        {
            // test if neighbor 36 exists and is on same level
            if( mNeighbors[ 17 ] != nullptr )
            {
                if( mNeighbors[ 17 ]->get_level() == mLevel )
                {
                   // copy edge 1 from neighbor
                   this->insert_edge( mNeighbors[ 17 ]->get_edge( tE1 ), tE11 );
                }
            }
        }

        // step 3: create edges that do not exist

        // test if edge 0 does not exist
        if( mEdges[ tE0 ] == nullptr )
        {
            // create edge
            mEdges[ tE0 ] = new Background_Edge( this, tE0 );

            // set owning flag
            mEdgeOwnFlags.set( tE0 );

            // test if neighbor 0 exists
            if( mNeighbors[ 0 ] != nullptr )
            {
                if( mNeighbors[ 0 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 0 ]->insert_edge( mEdges[ tE0 ], tE2 );
                }
            }

            // test if neighbor tE4 exists
            if( mNeighbors[ 4 ] != nullptr )
            {
                if( mNeighbors[ 4 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 4 ]->insert_edge( mEdges[ tE0 ], tE8 );
                }
            }

            // test if neighbor tE6 exists
            if( mNeighbors[ 6 ] != nullptr )
            {
                if( mNeighbors[ 6 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 6 ]->insert_edge( mEdges[ tE0 ], tE10 );
                }
            }
        }
        // test if edge 1 does not exist
        if( mEdges[ tE1 ] == nullptr )
        {
            // create edge
            mEdges[ tE1 ] = new Background_Edge( this, tE1 );

            // set owning flag
            mEdgeOwnFlags.set( tE1 );

            // test if neighbor 1 exists
            if( mNeighbors[ 1 ] != nullptr )
            {
                if( mNeighbors[ 1 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 1 ]->insert_edge( mEdges[ tE1 ], tE3 );
                }
            }

            // test if neighbor tE4 exists
            if( mNeighbors[ 4 ] != nullptr )
            {
                if( mNeighbors[ 4 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 4 ]->insert_edge( mEdges[ tE1 ], tE9 );
                }
            }

            // test if neighbor tE7 exists
            if( mNeighbors[ 7 ] != nullptr )
            {
                if( mNeighbors[ 7 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 7 ]->insert_edge( mEdges[ tE1 ], tE11 );
                }
            }
        }
        // test if edge 2 does not exist
        if( mEdges[ tE2 ] == nullptr )
        {
            // create edge
            mEdges[ tE2 ] = new Background_Edge( this, tE2 );

            // set owning flag
            mEdgeOwnFlags.set( tE2 );

            // test if neighbor tE2 exists
            if( mNeighbors[ 2 ] != nullptr )
            {
                if( mNeighbors[ 2 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 2 ]->insert_edge( mEdges[ tE2 ], tE0 );
                }
            }

            // test if neighbor tE4 exists
            if( mNeighbors[ 4 ] != nullptr )
            {
                if( mNeighbors[ 4 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 4 ]->insert_edge( mEdges[ tE2 ], tE10 );
                }
            }

            // test if neighbor tE8 exists
            if( mNeighbors[ 8 ] != nullptr )
            {
                if( mNeighbors[ 8 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 8 ]->insert_edge( mEdges[ tE2 ], tE8 );
                }
            }
        }
        // test if edge 3 does not exist
        if( mEdges[ tE3 ] == nullptr )
        {
            // create edge
            mEdges[ tE3 ] = new Background_Edge( this, tE3 );

            // set owning flag
            mEdgeOwnFlags.set( tE3 );

            // test if neighbor 3 exists
            if( mNeighbors[ 3 ] != nullptr )
            {
                if( mNeighbors[ 3 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 3 ]->insert_edge( mEdges[ tE3 ], tE1 );
                }
            }

            // test if neighbor tE4 exists
            if( mNeighbors[ 4 ] != nullptr )
            {
                if( mNeighbors[ 4 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 4 ]->insert_edge( mEdges[ tE3 ], tE11 );
                }
            }

            // test if neighbor tE9 exists
            if( mNeighbors[ 9 ] != nullptr )
            {
                if( mNeighbors[ 9 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 9 ]->insert_edge( mEdges[ tE3 ], tE9 );
                }
            }
        }
        // test if edge tE4 does not exist
        if( mEdges[ tE4 ] == nullptr )
        {
            // create edge
            mEdges[ tE4 ] = new Background_Edge( this, tE4 );

            // set owning flag
            mEdgeOwnFlags.set( tE4 );

            // test if neighbor 0 exists
            if( mNeighbors[ 0 ] != nullptr )
            {
                if( mNeighbors[ 0 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 0 ]->insert_edge( mEdges[ tE4 ], tE7 );
                }
            }

            // test if neighbor 3 exists
            if( mNeighbors[ 3 ] != nullptr )
            {
                if( mNeighbors[ 3 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 3 ]->insert_edge( mEdges[ tE4 ], tE5 );
                }
            }

            // test if neighbor tE10 exists
            if( mNeighbors[ 10 ] != nullptr )
            {
                if( mNeighbors[ 10 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 10 ]->insert_edge( mEdges[ tE4 ], tE6 );
                }
            }
        }
        // test if edge tE5 does not exist
        if( mEdges[ tE5 ] == nullptr )
        {
            // create edge
            mEdges[ tE5 ] = new Background_Edge( this, tE5 );

            // set owning flag
            mEdgeOwnFlags.set( tE5 );

            // test if neighbor 0 exists
            if( mNeighbors[ 0 ] != nullptr )
            {
                if( mNeighbors[ 0 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 0 ]->insert_edge( mEdges[ tE5 ], tE6 );
                }
            }

            // test if neighbor 1 exists
            if( mNeighbors[ 1 ] != nullptr )
            {
                if( mNeighbors[ 1 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 1 ]->insert_edge( mEdges[ tE5 ], tE4 );
                }
            }

            // test if neighbor tE11 exists
            if( mNeighbors[ 11 ] != nullptr )
            {
                if( mNeighbors[ 11 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 11 ]->insert_edge( mEdges[ tE5 ], tE7 );
                }
            }
        }
        // test if edge tE6 does not exist
        if( mEdges[ tE6 ] == nullptr )
        {
            // create edge
            mEdges[ tE6 ] = new Background_Edge( this, tE6 );

            // set owning flag
            mEdgeOwnFlags.set( tE6 );

            // test if neighbor 1 exists
            if( mNeighbors[ 1 ] != nullptr )
            {
                if( mNeighbors[ 1 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 1 ]->insert_edge( mEdges[ tE6 ], tE7 );
                }
            }

            // test if neighbor 2 exists
            if( mNeighbors[ 2 ] != nullptr )
            {
                if( mNeighbors[ 2 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 2 ]->insert_edge( mEdges[ tE6 ], tE5 );
                }
            }

            // test if neighbor 12 exists
            if( mNeighbors[ 12 ] != nullptr )
            {
                if( mNeighbors[ 12 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 12 ]->insert_edge( mEdges[ tE6 ], tE4 );
                }
            }
        }
        // test if edge tE7 does not exist
        if( mEdges[ tE7 ] == nullptr )
        {
            // create edge
            mEdges[ tE7 ] = new Background_Edge( this, tE7 );

            // set owning flag
            mEdgeOwnFlags.set( tE7 );

            // test if neighbor 2 exists
            if( mNeighbors[ 2 ] != nullptr )
            {
                if( mNeighbors[ 2 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 2 ]->insert_edge( mEdges[ tE7 ], tE4 );
                }
            }

            // test if neighbor 3 exists
            if( mNeighbors[ 3 ] != nullptr )
            {
                if( mNeighbors[ 3 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 3 ]->insert_edge( mEdges[ tE7 ], tE6 );
                }
            }

            // test if neighbor 13 exists
            if( mNeighbors[ 13 ] != nullptr )
            {
                if( mNeighbors[ 13 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 13 ]->insert_edge( mEdges[ tE7 ], tE5 );
                }
            }
        }
        // test if edge tE8 does not exist
        if( mEdges[ tE8 ] == nullptr )
        {
            // create edge
            mEdges[ tE8 ] = new Background_Edge( this, tE8 );

            // set owning flag
            mEdgeOwnFlags.set( tE8 );

            // test if neighbor 0 exists
            if( mNeighbors[ 0 ] != nullptr )
            {
                if( mNeighbors[ 0 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 0 ]->insert_edge( mEdges[ tE8 ], tE10 );
                }
            }

            // test if neighbor tE5 exists
            if( mNeighbors[ 5 ] != nullptr )
            {
                if( mNeighbors[ 5 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 5 ]->insert_edge( mEdges[ tE8 ], tE0 );
                }
            }

            // test if neighbor 14 exists
            if( mNeighbors[ 14 ] != nullptr )
            {
                if( mNeighbors[ 14 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 14 ]->insert_edge( mEdges[ tE8 ], tE2 );
                }
            }
        }
        // test if edge tE9 does not exist
        if( mEdges[ tE9 ] == nullptr )
        {
            // create edge
            mEdges[ tE9 ] = new Background_Edge( this, tE9 );

            // set owning flag
            mEdgeOwnFlags.set( tE9 );

            // test if neighbor 1 exists
            if( mNeighbors[ 1 ] != nullptr )
            {
                if( mNeighbors[ 1 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 1 ]->insert_edge( mEdges[ tE9 ], tE11 );
                }
            }

            // test if neighbor tE5 exists
            if( mNeighbors[ 5 ] != nullptr )
            {
                if( mNeighbors[ 5 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 5 ]->insert_edge( mEdges[ tE9 ], tE1 );
                }
            }

            // test if neighbor 15 exists
            if( mNeighbors[ 15 ] != nullptr )
            {
                if( mNeighbors[ 15 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 15 ]->insert_edge( mEdges[ tE9 ], tE3 );
                }
            }
        }
        // test if edge tE10 does not exist
        if( mEdges[ tE10 ] == nullptr )
        {
            // create edge
            mEdges[ tE10 ] = new Background_Edge( this, tE10 );

            // set owning flag
            mEdgeOwnFlags.set( tE10 );

            // test if neighbor 2 exists
            if( mNeighbors[ 2 ] != nullptr )
            {
                if( mNeighbors[ 2 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 2 ]->insert_edge( mEdges[ tE10 ], tE8 );
                }
            }

            // test if neighbor tE5 exists
            if( mNeighbors[ 5 ] != nullptr )
            {
                if( mNeighbors[ 5 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 5 ]->insert_edge( mEdges[ tE10 ], tE2 );
                }
            }

            // test if neighbor 16 exists
            if( mNeighbors[ 16 ] != nullptr )
            {
                if( mNeighbors[ 16 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 16 ]->insert_edge( mEdges[ tE10 ], tE0 );
                }
            }
        }
        // test if edge tE11 does not exist
        if( mEdges[ tE11 ] == nullptr )
        {
            // create edge
            mEdges[ tE11 ] = new Background_Edge( this, tE11 );

            // set owning flag
            mEdgeOwnFlags.set( tE11 );

            // test if neighbor 3 exists
            if( mNeighbors[ 3 ] != nullptr )
            {
                if( mNeighbors[ 3 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 3 ]->insert_edge( mEdges[ tE11 ], tE9 );
                }
            }

            // test if neighbor tE5 exists
            if( mNeighbors[ 5 ] != nullptr )
            {
                if( mNeighbors[ 5 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 5 ]->insert_edge( mEdges[ tE11 ], tE3 );
                }
            }

            // test if neighbor 17 exists
            if( mNeighbors[ 17 ] != nullptr )
            {
                if( mNeighbors[ 17 ]->get_level() == mLevel )
                {
                    // copy edge to neighbor
                    mNeighbors[ 17 ]->insert_edge( mEdges[ tE11 ], tE1 );
                }
            }
        }
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_EDGES_3D_HPP_ */

