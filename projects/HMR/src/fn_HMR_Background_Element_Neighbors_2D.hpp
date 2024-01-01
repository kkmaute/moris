/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_Background_Element_Neighbors_2D.hpp
 *
 */

#ifndef SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_
#define SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_

#include "cl_HMR_Background_Element.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::hmr
{
    // ----------------------------------------------------------------------------

    template <>
    inline
    void Background_Element< 2 >::get_neighbors_from_same_level(
            uint aOrder,
            Vector< Background_Element_Base * > & aNeighbors )
    {
        // make sure order is not too big
        MORIS_ERROR( 0 < aOrder && aOrder <= 3,
                "Invalid refinement buffer size specified; valid refinement buffer sizes are 1, 2, and 3.");

        // array that contains max size
        uint tArraySize[ 4 ] = { 0, 8, 24, 48 };

        // initialize temporary neighbor array
        Vector< Background_Element_Base * >
        tNeighbors( tArraySize[ aOrder ], nullptr );

        // fill first frame
        for ( uint k = 0; k < 8; ++k)
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
                if ( tNeighbor != nullptr && tNeighbors( 10 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 10 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 9 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 9 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 0
                tNeighbor =  tNeighbors( 0 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 11 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 11 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 1 exists
            if ( tNeighbors( 1 ) != nullptr )
            {
                // get neighbor 1 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 16 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 16 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 14 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 14 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 1
                tNeighbor =  tNeighbors( 1 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 18 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 18 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 2 exists
            if ( tNeighbors( 2 ) != nullptr )
            {
                // get neighbor 2 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 21 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 21 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 22 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 22 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 2
                tNeighbor =  tNeighbors( 2 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 20 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 20 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 3 exists
            if ( tNeighbors( 3 ) != nullptr )
            {
                // get neighbor 3 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 15 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 15 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 13 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 13 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 3
                tNeighbor =  tNeighbors( 3 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 17 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 17 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 4 exists
            if ( tNeighbors( 4 ) != nullptr )
            {
                // get neighbor 0 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 9 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 9 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 13 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 13 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 8 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 8 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 10 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 10 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 4
                tNeighbor =  tNeighbors( 4 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 15 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 15 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 5 exists
            if ( tNeighbors( 5 ) != nullptr )
            {
                // get neighbor 0 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 0 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 11 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 11 ) = tNeighbor;
                    }
                }

                // get neighbor 1 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 14 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 14 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 10 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 10 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 12 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 12 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 5
                tNeighbor =  tNeighbors( 5 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 16 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 16 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 6 exists
            if ( tNeighbors( 6 ) != nullptr )
            {
                // get neighbor 1 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 1 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 18 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 18 ) = tNeighbor;
                    }
                }

                // get neighbor 2 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 22 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 22 ) = tNeighbor;
                    }
                }

                // get neighbor 5 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 5 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 16 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 16 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 23 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 23 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 6
                tNeighbor =  tNeighbors( 6 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 21 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 21 ) = tNeighbor;
                    }
                }

            }
            // test if neighbor 7 exists
            if ( tNeighbors( 7 ) != nullptr )
            {
                // get neighbor 2 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 2 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 20 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 20 ) = tNeighbor;
                    }
                }

                // get neighbor 3 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 3 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 17 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 17 ) = tNeighbor;
                    }
                }

                // get neighbor 4 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 4 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 15 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 15 ) = tNeighbor;
                    }
                }

                // get neighbor 6 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 6 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 21 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 21 ) = tNeighbor;
                    }
                }

                // get neighbor 7 of neighbor 7
                tNeighbor =  tNeighbors( 7 )->get_neighbor( 7 );

                // test if neighbor exists and was not copied yet
                if ( tNeighbor != nullptr && tNeighbors( 19 ) == nullptr )
                {
                    // test if neighbor is on same level
                    if ( tNeighbor->get_level() == mLevel )
                    {
                        // copy pointer in big array
                        tNeighbors( 19 ) = tNeighbor;
                    }
                }

            }
            if ( aOrder >= 3 )
            {
                // test if neighbor 8 exists
                if ( tNeighbors( 8 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 8
                    tNeighbor =  tNeighbors( 8 )->get_neighbor( 0 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 25 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 25 ) = tNeighbor;
                        }
                    }

                    // get neighbor 3 of neighbor 8
                    tNeighbor =  tNeighbors( 8 )->get_neighbor( 3 );

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

                    // get neighbor 4 of neighbor 8
                    tNeighbor =  tNeighbors( 8 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 24 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 24 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 8
                    tNeighbor =  tNeighbors( 8 )->get_neighbor( 5 );

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

                    // get neighbor 7 of neighbor 8
                    tNeighbor =  tNeighbors( 8 )->get_neighbor( 7 );

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

                }
                // test if neighbor 9 exists
                if ( tNeighbors( 9 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 9
                    tNeighbor =  tNeighbors( 9 )->get_neighbor( 0 );

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

                    // get neighbor 4 of neighbor 9
                    tNeighbor =  tNeighbors( 9 )->get_neighbor( 4 );

                    // test if neighbor exists and was not copied yet
                    if ( tNeighbor != nullptr && tNeighbors( 25 ) == nullptr )
                    {
                        // test if neighbor is on same level
                        if ( tNeighbor->get_level() == mLevel )
                        {
                            // copy pointer in big array
                            tNeighbors( 25 ) = tNeighbor;
                        }
                    }

                    // get neighbor 5 of neighbor 9
                    tNeighbor =  tNeighbors( 9 )->get_neighbor( 5 );

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

                }
                // test if neighbor 10 exists
                if ( tNeighbors( 10 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 10
                    tNeighbor =  tNeighbors( 10 )->get_neighbor( 0 );

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

                    // get neighbor 4 of neighbor 10
                    tNeighbor =  tNeighbors( 10 )->get_neighbor( 4 );

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

                    // get neighbor 5 of neighbor 10
                    tNeighbor =  tNeighbors( 10 )->get_neighbor( 5 );

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

                }
                // test if neighbor 11 exists
                if ( tNeighbors( 11 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 11
                    tNeighbor =  tNeighbors( 11 )->get_neighbor( 0 );

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

                    // get neighbor 4 of neighbor 11
                    tNeighbor =  tNeighbors( 11 )->get_neighbor( 4 );

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

                    // get neighbor 5 of neighbor 11
                    tNeighbor =  tNeighbors( 11 )->get_neighbor( 5 );

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

                }
                // test if neighbor 12 exists
                if ( tNeighbors( 12 ) != nullptr )
                {
                    // get neighbor 0 of neighbor 12
                    tNeighbor =  tNeighbors( 12 )->get_neighbor( 0 );

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

                    // get neighbor 1 of neighbor 12
                    tNeighbor =  tNeighbors( 12 )->get_neighbor( 1 );

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

                    // get neighbor 4 of neighbor 12
                    tNeighbor =  tNeighbors( 12 )->get_neighbor( 4 );

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

                    // get neighbor 5 of neighbor 12
                    tNeighbor =  tNeighbors( 12 )->get_neighbor( 5 );

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

                    // get neighbor 6 of neighbor 12
                    tNeighbor =  tNeighbors( 12 )->get_neighbor( 6 );

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

                }
                // test if neighbor 13 exists
                if ( tNeighbors( 13 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 13
                    tNeighbor =  tNeighbors( 13 )->get_neighbor( 3 );

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

                    // get neighbor 4 of neighbor 13
                    tNeighbor =  tNeighbors( 13 )->get_neighbor( 4 );

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

                    // get neighbor 7 of neighbor 13
                    tNeighbor =  tNeighbors( 13 )->get_neighbor( 7 );

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

                }
                // test if neighbor 14 exists
                if ( tNeighbors( 14 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 14
                    tNeighbor =  tNeighbors( 14 )->get_neighbor( 1 );

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

                    // get neighbor 5 of neighbor 14
                    tNeighbor =  tNeighbors( 14 )->get_neighbor( 5 );

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

                    // get neighbor 6 of neighbor 14
                    tNeighbor =  tNeighbors( 14 )->get_neighbor( 6 );

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

                }
                // test if neighbor 15 exists
                if ( tNeighbors( 15 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 15
                    tNeighbor =  tNeighbors( 15 )->get_neighbor( 3 );

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

                    // get neighbor 4 of neighbor 15
                    tNeighbor =  tNeighbors( 15 )->get_neighbor( 4 );

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

                    // get neighbor 7 of neighbor 15
                    tNeighbor =  tNeighbors( 15 )->get_neighbor( 7 );

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

                }
                // test if neighbor 16 exists
                if ( tNeighbors( 16 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 16
                    tNeighbor =  tNeighbors( 16 )->get_neighbor( 1 );

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

                    // get neighbor 5 of neighbor 16
                    tNeighbor =  tNeighbors( 16 )->get_neighbor( 5 );

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

                    // get neighbor 6 of neighbor 16
                    tNeighbor =  tNeighbors( 16 )->get_neighbor( 6 );

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

                }
                // test if neighbor 17 exists
                if ( tNeighbors( 17 ) != nullptr )
                {
                    // get neighbor 3 of neighbor 17
                    tNeighbor =  tNeighbors( 17 )->get_neighbor( 3 );

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

                    // get neighbor 4 of neighbor 17
                    tNeighbor =  tNeighbors( 17 )->get_neighbor( 4 );

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

                    // get neighbor 7 of neighbor 17
                    tNeighbor =  tNeighbors( 17 )->get_neighbor( 7 );

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

                }
                // test if neighbor 18 exists
                if ( tNeighbors( 18 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 18
                    tNeighbor =  tNeighbors( 18 )->get_neighbor( 1 );

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

                    // get neighbor 5 of neighbor 18
                    tNeighbor =  tNeighbors( 18 )->get_neighbor( 5 );

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

                    // get neighbor 6 of neighbor 18
                    tNeighbor =  tNeighbors( 18 )->get_neighbor( 6 );

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

                }
                // test if neighbor 19 exists
                if ( tNeighbors( 19 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 19
                    tNeighbor =  tNeighbors( 19 )->get_neighbor( 2 );

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

                    // get neighbor 3 of neighbor 19
                    tNeighbor =  tNeighbors( 19 )->get_neighbor( 3 );

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

                    // get neighbor 4 of neighbor 19
                    tNeighbor =  tNeighbors( 19 )->get_neighbor( 4 );

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

                    // get neighbor 6 of neighbor 19
                    tNeighbor =  tNeighbors( 19 )->get_neighbor( 6 );

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

                    // get neighbor 7 of neighbor 19
                    tNeighbor =  tNeighbors( 19 )->get_neighbor( 7 );

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

                }
                // test if neighbor 20 exists
                if ( tNeighbors( 20 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 20
                    tNeighbor =  tNeighbors( 20 )->get_neighbor( 2 );

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

                    // get neighbor 6 of neighbor 20
                    tNeighbor =  tNeighbors( 20 )->get_neighbor( 6 );

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

                    // get neighbor 7 of neighbor 20
                    tNeighbor =  tNeighbors( 20 )->get_neighbor( 7 );

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
                // test if neighbor 21 exists
                if ( tNeighbors( 21 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 21
                    tNeighbor =  tNeighbors( 21 )->get_neighbor( 2 );

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

                    // get neighbor 6 of neighbor 21
                    tNeighbor =  tNeighbors( 21 )->get_neighbor( 6 );

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

                }
                // test if neighbor 22 exists
                if ( tNeighbors( 22 ) != nullptr )
                {
                    // get neighbor 2 of neighbor 22
                    tNeighbor =  tNeighbors( 22 )->get_neighbor( 2 );

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

                    // get neighbor 6 of neighbor 22
                    tNeighbor =  tNeighbors( 22 )->get_neighbor( 6 );

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

                    // get neighbor 7 of neighbor 22
                    tNeighbor =  tNeighbors( 22 )->get_neighbor( 7 );

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

                }
                // test if neighbor 23 exists
                if ( tNeighbors( 23 ) != nullptr )
                {
                    // get neighbor 1 of neighbor 23
                    tNeighbor =  tNeighbors( 23 )->get_neighbor( 1 );

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

                    // get neighbor 2 of neighbor 23
                    tNeighbor =  tNeighbors( 23 )->get_neighbor( 2 );

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

                    // get neighbor 5 of neighbor 23
                    tNeighbor =  tNeighbors( 23 )->get_neighbor( 5 );

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

                    // get neighbor 6 of neighbor 23
                    tNeighbor =  tNeighbors( 23 )->get_neighbor( 6 );

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

                    // get neighbor 7 of neighbor 23
                    tNeighbor =  tNeighbors( 23 )->get_neighbor( 7 );

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

#endif /* SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_ */

