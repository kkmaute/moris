/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Background_Element_Neighbors_3D.hpp
 *
 */

#pragma once

#include "cl_HMR_Background_Element.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::hmr
{
    template <>
    inline
    void Background_Element< 3 >::get_neighbors_from_same_level(
            uint aOrder,
            Vector< Background_Element_Base * > & aNeighbors )
    {
        // make sure order is not too big
        MORIS_ERROR( 0 < aOrder && aOrder <= 3,
                "Invalid refinement buffer size specified; valid refinement buffer sizes are 1, 2, and 3.");

        // array that contains max size
        uint tArraySize[ 4 ] = { 0, 26, 124, 342 };

        // initialize temporary neighbor array
        Vector< Background_Element_Base * >
        tNeighbors( tArraySize[ aOrder ], nullptr );

        // fill first frame
        for ( uint k=0; k<26; ++k)
        {
            // test if neighbor exists
            if ( mNeighbors[ k ] != nullptr )
            {
                // test if neighbor is on same level
                if ( mNeighbors[ k ]->get_level() == mLevel )
                {
                    // copy neighbor into array
                    tNeighbors( k ) = mNeighbors[ k ];
                }
            }
        }

        // temporary variable containing a neighbor
        Background_Element_Base * tNeighbor = nullptr;

        if ( aOrder >= 2 )
        {
            // test if neighbor 0 exists
            if ( tNeighbors( 0 ) != nullptr )
            {
                // get neighbor 0 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 52 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 52 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 54 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 54 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 84 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 84 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 86 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 86 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 1 exists
            if ( tNeighbors( 1 ) != nullptr )
            {
                // get neighbor 1 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 57 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 57 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 61 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 61 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 89 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 89 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 93 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 93 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 2 exists
            if ( tNeighbors( 2 ) != nullptr )
            {
                // get neighbor 2 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 65 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 65 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 63 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 63 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 97 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 97 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 95 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 95 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 3 exists
            if ( tNeighbors( 3 ) != nullptr )
            {
                // get neighbor 3 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 56 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 56 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 60 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 60 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 88 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 88 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 92 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 92 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 4 exists
            if ( tNeighbors( 4 ) != nullptr )
            {
                // get neighbor 4 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 32 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 32 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 34 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 34 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 44 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 44 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 42 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 42 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 5 exists
            if ( tNeighbors( 5 ) != nullptr )
            {
                // get neighbor 5 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 105 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 105 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 107 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 107 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 117 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 117 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 115 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 115 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 6 exists
            if ( tNeighbors( 6 ) != nullptr )
            {
                // get neighbor 0 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 28 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 28 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 34 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 34 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 32 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 32 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 52 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 52 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 54 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 54 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 27 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 27 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 29 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 29 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 7 exists
            if ( tNeighbors( 7 ) != nullptr )
            {
                // get neighbor 1 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 34 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 34 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 40 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 40 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 44 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 44 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 57 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 57 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 61 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 61 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 35 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 35 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 45 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 45 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 8 exists
            if ( tNeighbors( 8 ) != nullptr )
            {
                // get neighbor 2 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 44 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 44 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 48 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 48 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 42 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 42 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 65 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 65 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 63 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 63 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 49 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 49 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 47 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 47 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 8
                tNeighbor =  tNeighbors( 8 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 9 exists
            if ( tNeighbors( 9 ) != nullptr )
            {
                // get neighbor 3 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 32 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 32 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 42 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 42 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 36 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 36 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 56 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 56 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 60 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 60 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 31 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 31 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 41 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 41 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 9
                tNeighbor =  tNeighbors( 9 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 10 exists
            if ( tNeighbors( 10 ) != nullptr )
            {
                // get neighbor 0 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 52 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 52 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 56 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 56 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 67 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 67 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 84 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 84 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 88 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 88 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 51 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 51 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 83 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 83 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 10
                tNeighbor =  tNeighbors( 10 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 11 exists
            if ( tNeighbors( 11 ) != nullptr )
            {
                // get neighbor 0 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

                // get neighbor 1 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 54 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 54 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 57 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 57 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 71 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 71 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 86 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 86 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 89 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 89 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 55 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 55 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 87 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 87 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 11
                tNeighbor =  tNeighbors( 11 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 12 exists
            if ( tNeighbors( 12 ) != nullptr )
            {
                // get neighbor 1 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

                // get neighbor 2 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 61 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 61 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 65 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 65 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 82 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 82 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 93 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 93 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 97 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 97 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 66 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 66 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 98 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 98 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 12
                tNeighbor =  tNeighbors( 12 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 13 exists
            if ( tNeighbors( 13 ) != nullptr )
            {
                // get neighbor 2 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 63 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 63 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 60 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 60 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 78 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 78 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 95 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 95 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 92 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 92 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 62 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 62 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 13
                tNeighbor =  tNeighbors( 13 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 94 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 94 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 14 exists
            if ( tNeighbors( 14 ) != nullptr )
            {
                // get neighbor 0 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 84 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 84 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 86 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 86 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 101 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 101 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 107 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 107 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 105 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 105 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 100 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 100 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 102 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 102 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 14
                tNeighbor =  tNeighbors( 14 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 15 exists
            if ( tNeighbors( 15 ) != nullptr )
            {
                // get neighbor 1 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 89 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 89 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 93 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 93 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 107 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 107 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 113 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 113 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 117 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 117 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 108 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 108 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 118 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 118 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 15
                tNeighbor =  tNeighbors( 15 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 16 exists
            if ( tNeighbors( 16 ) != nullptr )
            {
                // get neighbor 2 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 97 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 97 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 95 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 95 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 117 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 117 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 121 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 121 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 115 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 115 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 122 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 122 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 16
                tNeighbor =  tNeighbors( 16 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 120 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 120 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 17 exists
            if ( tNeighbors( 17 ) != nullptr )
            {
                // get neighbor 3 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 88 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 88 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 92 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 92 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 105 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 105 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 115 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 115 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 109 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 109 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 104 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 104 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 17
                tNeighbor =  tNeighbors( 17 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 114 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 114 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 18 exists
            if ( tNeighbors( 18 ) != nullptr )
            {
                // get neighbor 0 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 52 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 52 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 56 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 56 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 32 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 32 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 27 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 27 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 31 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 31 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 51 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 51 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 26 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 26 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 28 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 28 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 36 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 36 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 67 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 67 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 18
                tNeighbor =  tNeighbors( 18 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 19 exists
            if ( tNeighbors( 19 ) != nullptr )
            {
                // get neighbor 0 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 54 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 54 ) = tNeighbor;
                    }
                }

                // get neighbor 1 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 57 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 57 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 34 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 34 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 29 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 29 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 35 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 35 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 33 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 33 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 53 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 53 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 55 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 55 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 28 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 28 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 30 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 30 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 40 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 40 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 71 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 71 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 19
                tNeighbor =  tNeighbors( 19 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 20 exists
            if ( tNeighbors( 20 ) != nullptr )
            {
                // get neighbor 1 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 61 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 61 ) = tNeighbor;
                    }
                }

                // get neighbor 2 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 65 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 65 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 44 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 44 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 39 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 39 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 45 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 45 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 49 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 49 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 59 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 59 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 66 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 66 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 40 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 40 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 50 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 50 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 48 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 48 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 82 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 82 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 20
                tNeighbor =  tNeighbors( 20 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 21 exists
            if ( tNeighbors( 21 ) != nullptr )
            {
                // get neighbor 2 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 63 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 63 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 60 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 60 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 42 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 42 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 37 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 37 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 43 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 43 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 47 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 47 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 41 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 41 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 58 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 58 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 64 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 64 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 62 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 62 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 36 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 36 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 38 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 38 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 48 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 48 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 46 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 46 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 21
                tNeighbor =  tNeighbors( 21 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 78 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 78 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 22 exists
            if ( tNeighbors( 22 ) != nullptr )
            {
                // get neighbor 0 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 84 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 84 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 88 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 88 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 105 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 105 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 68 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 68 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 72 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 72 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 83 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 83 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 100 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 100 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 104 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 104 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 67 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 67 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 99 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 99 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 101 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 101 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 22
                tNeighbor =  tNeighbors( 22 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 109 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 109 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 23 exists
            if ( tNeighbors( 23 ) != nullptr )
            {
                // get neighbor 0 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 86 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 86 ) = tNeighbor;
                    }
                }

                // get neighbor 1 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 89 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 89 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 107 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 107 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 70 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 70 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 73 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 73 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 85 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 85 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 87 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 87 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 102 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 102 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 108 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 108 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 106 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 106 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 69 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 69 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 71 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 71 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 101 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 101 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 103 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 103 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 113 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 113 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 23
                tNeighbor =  tNeighbors( 23 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 24 exists
            if ( tNeighbors( 24 ) != nullptr )
            {
                // get neighbor 1 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 93 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 93 ) = tNeighbor;
                    }
                }

                // get neighbor 2 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 97 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 97 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 117 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 117 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 77 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 77 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 81 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 81 ) = tNeighbor;
                    }
                }

                // get neighbor 11 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 11 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 91 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 91 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 98 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 98 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 112 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 112 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 118 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 118 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 122 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 122 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

                // get neighbor 19 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 19 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 75 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 75 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 82 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 82 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 113 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 113 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 123 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 123 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 24
                tNeighbor =  tNeighbors( 24 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 121 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 121 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 25 exists
            if ( tNeighbors( 25 ) != nullptr )
            {
                // get neighbor 2 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 95 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 95 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 92 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 92 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 115 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 115 ) = tNeighbor;
                    }
                }

                // get neighbor 8 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 8 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 79 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 79 ) = tNeighbor;
                    }
                }

                // get neighbor 9 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 9 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 76 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 76 ) = tNeighbor;
                    }
                }

                // get neighbor 10 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 10 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 90 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 90 ) = tNeighbor;
                    }
                }

                // get neighbor 12 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 12 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 96 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 96 ) = tNeighbor;
                    }
                }

                // get neighbor 13 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 13 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 94 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 94 ) = tNeighbor;
                    }
                }

                // get neighbor 14 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 14 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 110 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 110 ) = tNeighbor;
                    }
                }

                // get neighbor 15 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 15 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 116 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 116 ) = tNeighbor;
                    }
                }

                // get neighbor 16 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 16 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 120 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 120 ) = tNeighbor;
                    }
                }

                // get neighbor 17 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 17 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 114 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 114 ) = tNeighbor;
                    }
                }

                // get neighbor 18 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 18 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 74 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 74 ) = tNeighbor;
                    }
                }

                // get neighbor 20 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 20 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 80 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 80 ) = tNeighbor;
                    }
                }

                // get neighbor 21 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 21 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 78 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 78 ) = tNeighbor;
                    }
                }

                // get neighbor 22 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 22 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 109 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 109 ) = tNeighbor;
                    }
                }

                // get neighbor 23 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 23 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 111 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 111 ) = tNeighbor;
                    }
                }

                // get neighbor 24 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 24 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 121 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 121 ) = tNeighbor;
                    }
                }

                // get neighbor 25 of neighbor 25
                tNeighbor =  tNeighbors( 25 )->get_neighbor( 25 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 119 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 119 ) = tNeighbor;
                    }
                }

            }
            if ( aOrder >= 3 )
            {
                // test if neighbor 26 exists
                if ( tNeighbors( 26 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 174 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 174 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 180 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 180 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 132 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 132 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 125 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 125 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 131 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 131 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 173 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 173 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 124 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 124 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 126 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 126 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 138 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 138 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 197 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 197 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 26
                    tNeighbor =  tNeighbors( 26 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 27 exists
                if ( tNeighbors( 27 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 126 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 126 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 132 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 132 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 174 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 174 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 125 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 125 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 127 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 127 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 27
                    tNeighbor =  tNeighbors( 27 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 28 exists
                if ( tNeighbors( 28 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 127 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 127 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 126 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 126 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 128 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 128 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 28
                    tNeighbor =  tNeighbors( 28 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 29 exists
                if ( tNeighbors( 29 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 128 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 128 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 136 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 136 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 178 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 178 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 127 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 127 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 129 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 129 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 29
                    tNeighbor =  tNeighbors( 29 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 30 exists
                if ( tNeighbors( 30 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 178 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 178 ) = tNeighbor;
                        }
                    }

                    // get neighbor 1 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 181 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 181 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 136 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 136 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 129 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 129 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 137 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 137 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 179 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 179 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 128 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 128 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 130 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 130 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 144 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 144 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 203 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 203 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 30
                    tNeighbor =  tNeighbors( 30 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 31 exists
                if ( tNeighbors( 31 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 132 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 132 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 138 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 138 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 180 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 180 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 131 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 131 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 145 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 145 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 31
                    tNeighbor =  tNeighbors( 31 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 32 exists
                if ( tNeighbors( 32 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 132 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 132 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 32
                    tNeighbor =  tNeighbors( 32 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 33 exists
                if ( tNeighbors( 33 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 133 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 133 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 33
                    tNeighbor =  tNeighbors( 33 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 34 exists
                if ( tNeighbors( 34 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 134 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 134 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 136 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 136 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 34
                    tNeighbor =  tNeighbors( 34 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 35 exists
                if ( tNeighbors( 35 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 136 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 136 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 144 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 144 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 181 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 181 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 135 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 135 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 137 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 137 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 151 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 151 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 35
                    tNeighbor =  tNeighbors( 35 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 36 exists
                if ( tNeighbors( 36 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 145 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 145 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 138 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 138 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 152 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 152 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 36
                    tNeighbor =  tNeighbors( 36 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 37 exists
                if ( tNeighbors( 37 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 139 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 139 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 37
                    tNeighbor =  tNeighbors( 37 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 38 exists
                if ( tNeighbors( 38 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 140 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 140 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 38
                    tNeighbor =  tNeighbors( 38 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 39 exists
                if ( tNeighbors( 39 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 141 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 141 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 39
                    tNeighbor =  tNeighbors( 39 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 40 exists
                if ( tNeighbors( 40 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 143 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 143 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 151 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 151 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 142 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 142 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 144 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 144 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 158 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 158 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 40
                    tNeighbor =  tNeighbors( 40 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 41 exists
                if ( tNeighbors( 41 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 160 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 160 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 152 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 152 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 188 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 188 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 145 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 145 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 159 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 159 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 41
                    tNeighbor =  tNeighbors( 41 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 42 exists
                if ( tNeighbors( 42 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 146 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 146 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 42
                    tNeighbor =  tNeighbors( 42 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 160 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 160 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 43 exists
                if ( tNeighbors( 43 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 147 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 147 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 43
                    tNeighbor =  tNeighbors( 43 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 44 exists
                if ( tNeighbors( 44 ) != nullptr )
                {
                    // get neighbor 4 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 148 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 148 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 164 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 164 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 44
                    tNeighbor =  tNeighbors( 44 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 45 exists
                if ( tNeighbors( 45 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 150 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 150 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 158 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 158 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 164 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 164 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 189 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 189 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 149 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 149 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 151 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 151 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 165 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 165 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 45
                    tNeighbor =  tNeighbors( 45 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 46 exists
                if ( tNeighbors( 46 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 191 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 191 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 188 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 188 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 160 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 160 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 167 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 167 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 159 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 159 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 190 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 190 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 152 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 152 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 168 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 168 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 166 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 166 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 46
                    tNeighbor =  tNeighbors( 46 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 214 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 214 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 47 exists
                if ( tNeighbors( 47 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 168 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 168 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 160 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 160 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 191 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 191 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 153 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 153 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 169 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 169 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 167 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 167 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 47
                    tNeighbor =  tNeighbors( 47 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 48 exists
                if ( tNeighbors( 48 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 169 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 169 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 161 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 161 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 154 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 154 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 170 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 170 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 168 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 168 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 48
                    tNeighbor =  tNeighbors( 48 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 49 exists
                if ( tNeighbors( 49 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 164 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 164 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 170 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 170 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 162 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 162 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 195 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 195 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 155 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 155 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 171 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 171 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 169 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 169 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 49
                    tNeighbor =  tNeighbors( 49 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 50 exists
                if ( tNeighbors( 50 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 189 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 189 ) = tNeighbor;
                        }
                    }

                    // get neighbor 2 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 195 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 195 ) = tNeighbor;
                        }
                    }

                    // get neighbor 4 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 164 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 164 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 157 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 157 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 165 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 165 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 171 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 171 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 163 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 163 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 196 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 196 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 156 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 156 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 158 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 158 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 172 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 172 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 170 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 170 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 220 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 220 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 50
                    tNeighbor =  tNeighbors( 50 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 51 exists
                if ( tNeighbors( 51 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 174 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 174 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 180 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 180 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 197 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 197 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 173 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 173 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 221 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 221 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 51
                    tNeighbor =  tNeighbors( 51 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 52 exists
                if ( tNeighbors( 52 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 174 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 174 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 52
                    tNeighbor =  tNeighbors( 52 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 53 exists
                if ( tNeighbors( 53 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 175 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 175 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 53
                    tNeighbor =  tNeighbors( 53 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 54 exists
                if ( tNeighbors( 54 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 176 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 176 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 178 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 178 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 54
                    tNeighbor =  tNeighbors( 54 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 55 exists
                if ( tNeighbors( 55 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                    // get neighbor 1 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 178 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 178 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 181 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 181 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 203 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 203 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 177 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 177 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 179 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 179 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 227 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 227 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 55
                    tNeighbor =  tNeighbors( 55 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 56 exists
                if ( tNeighbors( 56 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 180 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 180 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 56
                    tNeighbor =  tNeighbors( 56 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 57 exists
                if ( tNeighbors( 57 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 181 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 181 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 57
                    tNeighbor =  tNeighbors( 57 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 58 exists
                if ( tNeighbors( 58 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 182 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 182 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 58
                    tNeighbor =  tNeighbors( 58 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 59 exists
                if ( tNeighbors( 59 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 183 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 183 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 59
                    tNeighbor =  tNeighbors( 59 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 60 exists
                if ( tNeighbors( 60 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 184 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 184 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 188 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 188 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 60
                    tNeighbor =  tNeighbors( 60 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 61 exists
                if ( tNeighbors( 61 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 185 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 185 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 189 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 189 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 61
                    tNeighbor =  tNeighbors( 61 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 62 exists
                if ( tNeighbors( 62 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 191 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 191 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 188 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 188 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 214 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 214 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 186 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 186 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 190 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 190 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 62
                    tNeighbor =  tNeighbors( 62 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 238 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 238 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 63 exists
                if ( tNeighbors( 63 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 191 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 191 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 63
                    tNeighbor =  tNeighbors( 63 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 64 exists
                if ( tNeighbors( 64 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 192 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 192 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 64
                    tNeighbor =  tNeighbors( 64 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 65 exists
                if ( tNeighbors( 65 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 195 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 195 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 193 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 193 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 65
                    tNeighbor =  tNeighbors( 65 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 66 exists
                if ( tNeighbors( 66 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                    // get neighbor 2 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 189 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 189 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 195 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 195 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 220 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 220 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 187 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 187 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 196 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 196 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 194 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 194 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 244 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 244 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 66
                    tNeighbor =  tNeighbors( 66 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 67 exists
                if ( tNeighbors( 67 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 221 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 221 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 197 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 197 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 245 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 245 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 67
                    tNeighbor =  tNeighbors( 67 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 68 exists
                if ( tNeighbors( 68 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 198 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 198 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 68
                    tNeighbor =  tNeighbors( 68 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 69 exists
                if ( tNeighbors( 69 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 199 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 199 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 69
                    tNeighbor =  tNeighbors( 69 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 70 exists
                if ( tNeighbors( 70 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 200 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 200 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 70
                    tNeighbor =  tNeighbors( 70 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 71 exists
                if ( tNeighbors( 71 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                    // get neighbor 1 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 202 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 202 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 227 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 227 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 201 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 201 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 203 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 203 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 251 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 251 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 71
                    tNeighbor =  tNeighbors( 71 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 72 exists
                if ( tNeighbors( 72 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 204 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 204 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 72
                    tNeighbor =  tNeighbors( 72 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 73 exists
                if ( tNeighbors( 73 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 205 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 205 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 73
                    tNeighbor =  tNeighbors( 73 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 74 exists
                if ( tNeighbors( 74 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 206 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 206 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 74
                    tNeighbor =  tNeighbors( 74 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 75 exists
                if ( tNeighbors( 75 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 207 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 207 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 75
                    tNeighbor =  tNeighbors( 75 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 76 exists
                if ( tNeighbors( 76 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 208 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 208 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 76
                    tNeighbor =  tNeighbors( 76 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 77 exists
                if ( tNeighbors( 77 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 209 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 209 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 77
                    tNeighbor =  tNeighbors( 77 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 78 exists
                if ( tNeighbors( 78 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 212 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 212 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 238 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 238 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 210 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 210 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 214 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 214 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 78
                    tNeighbor =  tNeighbors( 78 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 262 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 262 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 79 exists
                if ( tNeighbors( 79 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 215 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 215 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 79
                    tNeighbor =  tNeighbors( 79 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 80 exists
                if ( tNeighbors( 80 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 216 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 216 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 80
                    tNeighbor =  tNeighbors( 80 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 81 exists
                if ( tNeighbors( 81 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 217 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 217 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 81
                    tNeighbor =  tNeighbors( 81 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 82 exists
                if ( tNeighbors( 82 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                    // get neighbor 2 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 213 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 213 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 219 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 219 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 244 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 244 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 211 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 211 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 220 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 220 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 218 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 218 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 268 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 268 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 82
                    tNeighbor =  tNeighbors( 82 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 83 exists
                if ( tNeighbors( 83 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 245 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 245 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 270 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 270 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 276 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 276 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 221 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 221 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 269 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 269 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 83
                    tNeighbor =  tNeighbors( 83 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 84 exists
                if ( tNeighbors( 84 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 222 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 222 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 270 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 270 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 84
                    tNeighbor =  tNeighbors( 84 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 85 exists
                if ( tNeighbors( 85 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 223 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 223 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 85
                    tNeighbor =  tNeighbors( 85 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 86 exists
                if ( tNeighbors( 86 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 224 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 224 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 86
                    tNeighbor =  tNeighbors( 86 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 274 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 274 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 87 exists
                if ( tNeighbors( 87 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                    // get neighbor 1 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 226 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 226 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 251 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 251 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 274 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 274 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 277 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 277 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 225 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 225 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 227 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 227 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 275 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 275 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 87
                    tNeighbor =  tNeighbors( 87 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 88 exists
                if ( tNeighbors( 88 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 228 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 228 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 276 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 276 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 88
                    tNeighbor =  tNeighbors( 88 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 89 exists
                if ( tNeighbors( 89 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 229 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 229 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 277 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 277 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 89
                    tNeighbor =  tNeighbors( 89 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 90 exists
                if ( tNeighbors( 90 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 230 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 230 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 90
                    tNeighbor =  tNeighbors( 90 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 91 exists
                if ( tNeighbors( 91 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 231 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 231 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 91
                    tNeighbor =  tNeighbors( 91 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 92 exists
                if ( tNeighbors( 92 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 232 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 232 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 92
                    tNeighbor =  tNeighbors( 92 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 284 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 284 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 93 exists
                if ( tNeighbors( 93 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 233 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 233 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 93
                    tNeighbor =  tNeighbors( 93 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 285 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 285 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 94 exists
                if ( tNeighbors( 94 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 236 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 236 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 262 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 262 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 287 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 287 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 284 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 284 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 234 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 234 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 238 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 238 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 94
                    tNeighbor =  tNeighbors( 94 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 286 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 286 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 95 exists
                if ( tNeighbors( 95 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 239 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 239 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 95
                    tNeighbor =  tNeighbors( 95 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 287 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 287 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 96 exists
                if ( tNeighbors( 96 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 240 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 240 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 96
                    tNeighbor =  tNeighbors( 96 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 97 exists
                if ( tNeighbors( 97 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 241 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 241 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 291 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 291 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 97
                    tNeighbor =  tNeighbors( 97 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 98 exists
                if ( tNeighbors( 98 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                    // get neighbor 2 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 237 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 237 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 243 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 243 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 268 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 268 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 285 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 285 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 291 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 291 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 235 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 235 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 244 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 244 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 242 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 242 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 292 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 292 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 98
                    tNeighbor =  tNeighbors( 98 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 99 exists
                if ( tNeighbors( 99 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 270 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 270 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 276 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 276 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 301 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 301 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 269 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 269 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 294 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 294 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 300 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 300 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 245 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 245 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 293 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 293 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 295 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 295 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 99
                    tNeighbor =  tNeighbors( 99 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 307 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 307 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 100 exists
                if ( tNeighbors( 100 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 270 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 270 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 295 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 295 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 301 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 301 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 246 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 246 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 294 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 294 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 296 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 296 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 100
                    tNeighbor =  tNeighbors( 100 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 101 exists
                if ( tNeighbors( 101 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 271 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 271 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 296 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 296 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 247 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 247 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 295 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 295 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 297 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 297 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 101
                    tNeighbor =  tNeighbors( 101 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 102 exists
                if ( tNeighbors( 102 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 272 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 272 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 274 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 274 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 297 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 297 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 305 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 305 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 248 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 248 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 296 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 296 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 298 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 298 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 102
                    tNeighbor =  tNeighbors( 102 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 103 exists
                if ( tNeighbors( 103 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 274 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 274 ) = tNeighbor;
                        }
                    }

                    // get neighbor 1 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 277 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 277 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 305 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 305 ) = tNeighbor;
                        }
                    }

                    // get neighbor 6 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 6 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 250 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 250 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 273 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 273 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 275 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 275 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 298 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 298 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 306 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 306 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 249 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 249 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 251 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 251 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 297 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 297 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 299 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 299 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 313 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 313 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 103
                    tNeighbor =  tNeighbors( 103 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 104 exists
                if ( tNeighbors( 104 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 276 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 276 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 301 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 301 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 307 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 307 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 252 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 252 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 300 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 300 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 104
                    tNeighbor =  tNeighbors( 104 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 314 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 314 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 105 exists
                if ( tNeighbors( 105 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 301 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 301 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 105
                    tNeighbor =  tNeighbors( 105 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 106 exists
                if ( tNeighbors( 106 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 302 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 302 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 106
                    tNeighbor =  tNeighbors( 106 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 107 exists
                if ( tNeighbors( 107 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 303 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 303 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 305 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 305 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 107
                    tNeighbor =  tNeighbors( 107 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 108 exists
                if ( tNeighbors( 108 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 277 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 277 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 305 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 305 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 313 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 313 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 253 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 253 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 304 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 304 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 306 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 306 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 320 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 320 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 108
                    tNeighbor =  tNeighbors( 108 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 109 exists
                if ( tNeighbors( 109 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 278 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 278 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 314 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 314 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 254 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 254 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 307 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 307 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 109
                    tNeighbor =  tNeighbors( 109 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 321 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 321 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 110 exists
                if ( tNeighbors( 110 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 308 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 308 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 110
                    tNeighbor =  tNeighbors( 110 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 111 exists
                if ( tNeighbors( 111 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 309 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 309 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 111
                    tNeighbor =  tNeighbors( 111 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 112 exists
                if ( tNeighbors( 112 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 310 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 310 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 112
                    tNeighbor =  tNeighbors( 112 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 113 exists
                if ( tNeighbors( 113 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 279 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 279 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 312 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 312 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 320 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 320 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 255 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 255 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 311 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 311 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 313 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 313 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 327 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 327 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 113
                    tNeighbor =  tNeighbors( 113 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 114 exists
                if ( tNeighbors( 114 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 280 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 280 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 284 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 284 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 329 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 329 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 321 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 321 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 256 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 256 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 314 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 314 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 114
                    tNeighbor =  tNeighbors( 114 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 328 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 328 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 115 exists
                if ( tNeighbors( 115 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 315 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 315 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 115
                    tNeighbor =  tNeighbors( 115 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 329 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 329 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 116 exists
                if ( tNeighbors( 116 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 316 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 316 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 116
                    tNeighbor =  tNeighbors( 116 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 117 exists
                if ( tNeighbors( 117 ) != nullptr )
                {
                    // get neighbor 5 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 317 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 317 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 333 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 333 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 117
                    tNeighbor =  tNeighbors( 117 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 118 exists
                if ( tNeighbors( 118 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 281 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 281 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 285 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 285 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 319 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 319 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 327 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 327 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 333 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 333 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 257 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 257 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 318 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 318 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 320 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 320 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 334 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 334 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 118
                    tNeighbor =  tNeighbors( 118 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 119 exists
                if ( tNeighbors( 119 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 287 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 287 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 3 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 284 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 284 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 329 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 329 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                    // get neighbor 9 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 9 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 260 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 260 ) = tNeighbor;
                        }
                    }

                    // get neighbor 10 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 10 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 282 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 282 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 286 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 286 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 336 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 336 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 328 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 328 ) = tNeighbor;
                        }
                    }

                    // get neighbor 18 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 18 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 258 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 258 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 262 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 262 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 321 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 321 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 337 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 337 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 119
                    tNeighbor =  tNeighbors( 119 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 335 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 335 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 120 exists
                if ( tNeighbors( 120 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 287 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 287 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 337 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 337 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 329 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 329 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 263 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 263 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 322 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 322 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 338 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 338 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 120
                    tNeighbor =  tNeighbors( 120 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 336 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 336 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 121 exists
                if ( tNeighbors( 121 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 288 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 288 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 338 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 338 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 330 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 330 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 264 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 264 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 323 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 323 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 339 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 339 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 121
                    tNeighbor =  tNeighbors( 121 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 337 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 337 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 122 exists
                if ( tNeighbors( 122 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 291 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 291 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 289 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 289 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 333 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 333 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 339 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 339 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 331 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 331 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 265 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 265 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 324 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 324 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 340 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 340 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 122
                    tNeighbor =  tNeighbors( 122 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 338 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 338 ) = tNeighbor;
                        }
                    }

                }
                // test if neighbor 123 exists
                if ( tNeighbors( 123 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 1 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 285 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 285 ) = tNeighbor;
                        }
                    }

                    // get neighbor 2 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 2 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 291 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 291 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 5 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 333 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 333 ) = tNeighbor;
                        }
                    }

                    // get neighbor 7 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 7 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 261 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 261 ) = tNeighbor;
                        }
                    }

                    // get neighbor 8 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 8 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 267 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 267 ) = tNeighbor;
                        }
                    }

                    // get neighbor 11 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 11 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 283 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 283 ) = tNeighbor;
                        }
                    }

                    // get neighbor 12 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 12 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 292 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 292 ) = tNeighbor;
                        }
                    }

                    // get neighbor 13 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 13 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 290 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 290 ) = tNeighbor;
                        }
                    }

                    // get neighbor 14 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 14 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 326 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 326 ) = tNeighbor;
                        }
                    }

                    // get neighbor 15 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 15 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 334 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 334 ) = tNeighbor;
                        }
                    }

                    // get neighbor 16 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 16 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 340 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 340 ) = tNeighbor;
                        }
                    }

                    // get neighbor 17 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 17 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 332 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 332 ) = tNeighbor;
                        }
                    }

                    // get neighbor 19 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 19 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 259 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 259 ) = tNeighbor;
                        }
                    }

                    // get neighbor 20 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 20 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 268 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 268 ) = tNeighbor;
                        }
                    }

                    // get neighbor 21 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 21 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 266 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 266 ) = tNeighbor;
                        }
                    }

                    // get neighbor 22 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 22 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 325 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 325 ) = tNeighbor;
                        }
                    }

                    // get neighbor 23 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 23 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 327 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 327 ) = tNeighbor;
                        }
                    }

                    // get neighbor 24 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 24 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 341 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 341 ) = tNeighbor;
                        }
                    }

                    // get neighbor 25 of neighbor 123
                    tNeighbor =  tNeighbors( 123 )->get_neighbor( 25 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 339 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 339 ) = tNeighbor;
                        }
                    }

                }
            } // end order 3
        } // end order 2

        // initialize element counter
        uint tCount = 0;

        // count number of existing elements
        for( auto tNeighbor : tNeighbors )
        {
            if ( tNeighbor != nullptr )
            {
                ++tCount;
            }
        }

        // allocate output Cell
        aNeighbors.resize( tCount, nullptr );

        // reset counter
        tCount = 0;

        // copy existing elements
        for( auto tNeighbor : tNeighbors )
        {
            if ( tNeighbor != nullptr )
            {
                aNeighbors( tCount++ ) = tNeighbor;
            }
        }
    }

    // ----------------------------------------------------------------------------
} /* namespace moris */
