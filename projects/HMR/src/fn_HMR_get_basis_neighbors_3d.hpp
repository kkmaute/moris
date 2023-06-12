/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_get_basis_neighbors_3d.hpp
 *
 */

#ifndef SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_
#define SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_

#include "cl_HMR_Basis.hpp"
#include "typedefs.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        void get_basis_neighbors_3d(       Basis  * aBasis,
                                     uint aOrder,
                                           Basis ** aNeighbors )
        {
             // make sure order is not too big
             MORIS_ASSERT( 0 < aOrder && aOrder <= 3, "Neighbor order too big.");

             // array that contains max size
             uint tArraySize[ 4 ] = { 0, 26, 124, 342 };

             if ( aOrder >= 2 )
             {
                 // fill first frame
                 for ( uint k=0; k<26; ++k)
                 {
                     aNeighbors[ k ] = aBasis->get_neighbor( k );
                 }

                 // initialize higher order neighbors with null
                 for( uint k=26; k<tArraySize[ aOrder ]; ++k )
                 {
                     aNeighbors[ k ] = nullptr;

                 }

                 // placeholder for basis neighbor
                 Basis* tNeighbor;

                 // test if neighbor 0 exists
                 if ( aNeighbors[ 0 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 52 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 54 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 84 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 86 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 1 exists
                 if ( aNeighbors[ 1 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 57 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 61 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 89 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 93 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 2 exists
                 if ( aNeighbors[ 2 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 65 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 63 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 97 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 95 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 3 exists
                 if ( aNeighbors[ 3 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 56 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 60 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 88 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 92 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 4 exists
                 if ( aNeighbors[ 4 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 5 exists
                 if ( aNeighbors[ 5 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 105 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 107 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 117 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 115 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 6 exists
                 if ( aNeighbors[ 6 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 52 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 54 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 27 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 27 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 29 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 29 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 7 exists
                 if ( aNeighbors[ 7 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 40 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 57 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 61 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 35 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 35 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 45 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 45 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 8 exists
                 if ( aNeighbors[ 8 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 48 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 65 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 63 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 49 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 49 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 47 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 47 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 9 exists
                 if ( aNeighbors[ 9 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 56 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 60 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 31 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 31 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 41 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 41 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 10 exists
                 if ( aNeighbors[ 10 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 52 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 56 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 67 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 84 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 88 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 51 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 51 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 83 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 83 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 11 exists
                 if ( aNeighbors[ 11 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 54 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 57 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 71 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 86 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 89 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 55 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 55 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 87 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 87 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 12 exists
                 if ( aNeighbors[ 12 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 61 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 65 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 82 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 93 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 97 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 66 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 66 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 98 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 98 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 13 exists
                 if ( aNeighbors[ 13 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 63 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 60 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 78 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 95 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 92 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 62 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 62 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 94 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 94 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 14 exists
                 if ( aNeighbors[ 14 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 84 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 86 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 101 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 107 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 105 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 100 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 100 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 102 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 102 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 15 exists
                 if ( aNeighbors[ 15 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 89 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 93 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 107 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 113 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 117 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 108 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 108 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 118 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 118 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 16 exists
                 if ( aNeighbors[ 16 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 97 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 95 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 117 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 121 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 115 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 122 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 122 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 120 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 120 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 17 exists
                 if ( aNeighbors[ 17 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 88 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 92 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 105 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 115 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 109 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 104 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 104 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 114 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 114 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 18 exists
                 if ( aNeighbors[ 18 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 52 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 56 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 27 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 27 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 31 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 31 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 51 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 51 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 26 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 26 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 67 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 19 exists
                 if ( aNeighbors[ 19 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 54 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 57 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 29 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 29 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 35 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 35 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 53 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 55 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 55 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 30 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 30 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 40 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 71 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 20 exists
                 if ( aNeighbors[ 20 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 61 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 65 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 45 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 45 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 49 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 49 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 59 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 66 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 66 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 40 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 50 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 50 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 48 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 82 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 21 exists
                 if ( aNeighbors[ 21 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 63 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 60 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 47 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 47 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 41 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 41 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 58 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 64 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 62 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 62 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 48 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 46 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 46 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 78 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 22 exists
                 if ( aNeighbors[ 22 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 84 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 88 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 105 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 68 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 72 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 83 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 83 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 100 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 100 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 104 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 104 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 67 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 99 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 99 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 101 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 109 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 23 exists
                 if ( aNeighbors[ 23 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 86 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 89 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 107 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 70 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 73 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 85 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 87 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 87 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 102 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 102 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 108 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 108 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 106 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 69 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 71 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 101 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 103 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 103 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 113 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 24 exists
                 if ( aNeighbors[ 24 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 93 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 97 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 117 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 77 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 81 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 91 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 98 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 98 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 112 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 118 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 118 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 122 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 122 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 75 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 82 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 113 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 123 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 123 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 121 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 25 exists
                 if ( aNeighbors[ 25 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 95 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 92 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 115 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 79 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 76 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 90 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 96 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 94 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 94 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 110 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 116 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 120 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 120 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 114 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 114 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 74 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 80 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 78 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 109 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 111 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 121 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 119 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 119 ] = tNeighbor;
                     }
                 }
             if ( aOrder >= 3 )
             {

                 // test if neighbor 26 exists
                 if ( aNeighbors[ 26 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 174 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 180 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 132 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 125 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 125 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 131 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 131 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 173 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 173 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 124 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 124 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 126 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 138 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 197 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 27 exists
                 if ( aNeighbors[ 27 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 126 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 132 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 174 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 125 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 125 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 127 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 28 exists
                 if ( aNeighbors[ 28 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 127 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 126 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 128 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 29 exists
                 if ( aNeighbors[ 29 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 128 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 136 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 178 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 127 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 129 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 129 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 30 exists
                 if ( aNeighbors[ 30 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 178 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 181 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 136 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 129 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 129 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 137 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 137 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 179 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 179 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 128 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 130 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 130 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 144 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 203 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 31 exists
                 if ( aNeighbors[ 31 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 132 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 138 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 180 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 131 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 131 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 145 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 32 exists
                 if ( aNeighbors[ 32 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 132 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 33 exists
                 if ( aNeighbors[ 33 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 133 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 34 exists
                 if ( aNeighbors[ 34 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 134 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 136 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 35 exists
                 if ( aNeighbors[ 35 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 136 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 144 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 181 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 135 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 137 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 137 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 151 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 36 exists
                 if ( aNeighbors[ 36 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 145 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 138 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 152 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 37 exists
                 if ( aNeighbors[ 37 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 139 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 38 exists
                 if ( aNeighbors[ 38 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 140 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 39 exists
                 if ( aNeighbors[ 39 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 141 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 40 exists
                 if ( aNeighbors[ 40 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 143 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 151 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 142 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 144 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 158 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 41 exists
                 if ( aNeighbors[ 41 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 160 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 152 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 188 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 145 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 159 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 159 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 42 exists
                 if ( aNeighbors[ 42 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 146 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 160 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 43 exists
                 if ( aNeighbors[ 43 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 147 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 44 exists
                 if ( aNeighbors[ 44 ] != nullptr )
                 {
                     // get neighbor 4 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 148 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 164 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 45 exists
                 if ( aNeighbors[ 45 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 150 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 158 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 164 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 189 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 149 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 151 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 165 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 165 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 46 exists
                 if ( aNeighbors[ 46 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 191 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 188 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 160 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 167 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 167 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 159 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 159 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 190 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 190 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 152 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 168 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 166 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 166 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 214 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 47 exists
                 if ( aNeighbors[ 47 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 168 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 160 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 191 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 153 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 169 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 167 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 167 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 48 exists
                 if ( aNeighbors[ 48 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 169 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 161 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 154 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 170 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 168 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 49 exists
                 if ( aNeighbors[ 49 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 164 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 170 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 162 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 195 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 155 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 171 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 171 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 169 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 50 exists
                 if ( aNeighbors[ 50 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 189 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 195 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 164 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 157 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 165 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 165 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 171 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 171 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 163 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 196 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 196 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 156 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 158 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 172 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 172 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 170 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 220 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 51 exists
                 if ( aNeighbors[ 51 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 174 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 180 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 197 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 173 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 173 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 221 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 52 exists
                 if ( aNeighbors[ 52 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 174 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 53 exists
                 if ( aNeighbors[ 53 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 175 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 54 exists
                 if ( aNeighbors[ 54 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 176 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 178 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 55 exists
                 if ( aNeighbors[ 55 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 178 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 181 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 203 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 177 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 179 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 179 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 227 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 56 exists
                 if ( aNeighbors[ 56 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 180 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 57 exists
                 if ( aNeighbors[ 57 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 181 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 58 exists
                 if ( aNeighbors[ 58 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 182 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 59 exists
                 if ( aNeighbors[ 59 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 183 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 60 exists
                 if ( aNeighbors[ 60 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 184 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 188 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 61 exists
                 if ( aNeighbors[ 61 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 185 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 189 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 62 exists
                 if ( aNeighbors[ 62 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 191 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 188 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 214 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 186 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 190 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 190 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 238 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 63 exists
                 if ( aNeighbors[ 63 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 191 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 64 exists
                 if ( aNeighbors[ 64 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 192 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 65 exists
                 if ( aNeighbors[ 65 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 195 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 193 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 66 exists
                 if ( aNeighbors[ 66 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 189 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 195 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 220 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 187 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 196 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 196 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 194 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 244 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 67 exists
                 if ( aNeighbors[ 67 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 221 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 197 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 245 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 68 exists
                 if ( aNeighbors[ 68 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 198 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 69 exists
                 if ( aNeighbors[ 69 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 199 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 70 exists
                 if ( aNeighbors[ 70 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 200 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 71 exists
                 if ( aNeighbors[ 71 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 202 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 227 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 201 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 203 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 251 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 72 exists
                 if ( aNeighbors[ 72 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 204 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 73 exists
                 if ( aNeighbors[ 73 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 205 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 74 exists
                 if ( aNeighbors[ 74 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 206 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 75 exists
                 if ( aNeighbors[ 75 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 207 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 76 exists
                 if ( aNeighbors[ 76 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 208 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 77 exists
                 if ( aNeighbors[ 77 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 209 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 78 exists
                 if ( aNeighbors[ 78 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 212 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 238 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 210 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 214 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 262 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 79 exists
                 if ( aNeighbors[ 79 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 215 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 80 exists
                 if ( aNeighbors[ 80 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 216 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 81 exists
                 if ( aNeighbors[ 81 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 217 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 82 exists
                 if ( aNeighbors[ 82 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 213 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 219 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 244 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 211 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 220 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 218 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 268 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 83 exists
                 if ( aNeighbors[ 83 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 245 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 270 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 276 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 221 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 269 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 269 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 84 exists
                 if ( aNeighbors[ 84 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 222 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 270 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 85 exists
                 if ( aNeighbors[ 85 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 223 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 86 exists
                 if ( aNeighbors[ 86 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 224 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 274 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 87 exists
                 if ( aNeighbors[ 87 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 226 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 251 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 274 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 277 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 225 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 227 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 275 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 275 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 88 exists
                 if ( aNeighbors[ 88 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 228 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 276 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 89 exists
                 if ( aNeighbors[ 89 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 229 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 277 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 90 exists
                 if ( aNeighbors[ 90 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 230 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 91 exists
                 if ( aNeighbors[ 91 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 231 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 92 exists
                 if ( aNeighbors[ 92 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 232 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 284 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 93 exists
                 if ( aNeighbors[ 93 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 233 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 285 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 94 exists
                 if ( aNeighbors[ 94 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 236 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 262 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 287 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 284 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 234 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 238 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 286 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 286 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 95 exists
                 if ( aNeighbors[ 95 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 239 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 287 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 96 exists
                 if ( aNeighbors[ 96 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 240 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 97 exists
                 if ( aNeighbors[ 97 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 241 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 291 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 98 exists
                 if ( aNeighbors[ 98 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 237 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 243 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 268 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 285 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 291 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 235 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 244 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 242 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 292 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 292 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 99 exists
                 if ( aNeighbors[ 99 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 270 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 276 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 301 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 269 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 269 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 294 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 294 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 300 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 300 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 245 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 293 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 293 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 295 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 307 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 100 exists
                 if ( aNeighbors[ 100 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 270 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 295 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 301 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 246 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 294 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 294 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 296 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 101 exists
                 if ( aNeighbors[ 101 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 271 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 296 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 247 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 295 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 297 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 102 exists
                 if ( aNeighbors[ 102 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 272 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 274 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 297 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 305 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 248 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 296 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 298 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 298 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 103 exists
                 if ( aNeighbors[ 103 ] != nullptr )
                 {
                     // get neighbor 0 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 274 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 277 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 305 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 250 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 273 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 275 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 275 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 298 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 298 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 306 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 306 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 249 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 251 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 297 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 299 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 299 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 313 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 104 exists
                 if ( aNeighbors[ 104 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 276 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 301 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 307 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 252 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 300 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 300 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 314 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 105 exists
                 if ( aNeighbors[ 105 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 301 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 106 exists
                 if ( aNeighbors[ 106 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 302 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 107 exists
                 if ( aNeighbors[ 107 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 303 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 305 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 108 exists
                 if ( aNeighbors[ 108 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 277 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 305 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 313 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 253 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 304 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 306 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 306 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 320 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 109 exists
                 if ( aNeighbors[ 109 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 278 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 314 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 254 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 307 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 321 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 110 exists
                 if ( aNeighbors[ 110 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 308 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 111 exists
                 if ( aNeighbors[ 111 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 309 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 112 exists
                 if ( aNeighbors[ 112 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 310 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 113 exists
                 if ( aNeighbors[ 113 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 279 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 312 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 320 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 255 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 311 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 313 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 327 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 114 exists
                 if ( aNeighbors[ 114 ] != nullptr )
                 {
                     // get neighbor 3 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 280 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 284 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 329 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 321 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 256 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 314 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 328 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 328 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 115 exists
                 if ( aNeighbors[ 115 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 315 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 329 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 116 exists
                 if ( aNeighbors[ 116 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 316 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 117 exists
                 if ( aNeighbors[ 117 ] != nullptr )
                 {
                     // get neighbor 5 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 317 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 333 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 118 exists
                 if ( aNeighbors[ 118 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 281 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 285 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 319 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 327 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 333 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 257 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 318 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 320 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 334 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 334 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 119 exists
                 if ( aNeighbors[ 119 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 287 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 284 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 329 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 260 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 282 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 286 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 286 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 336 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 336 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 328 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 328 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 258 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 262 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 321 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 337 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 335 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 335 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 120 exists
                 if ( aNeighbors[ 120 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 287 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 337 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 329 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 263 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 322 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 338 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 336 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 336 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 121 exists
                 if ( aNeighbors[ 121 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 288 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 338 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 330 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 264 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 323 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 339 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 339 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 337 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 122 exists
                 if ( aNeighbors[ 122 ] != nullptr )
                 {
                     // get neighbor 2 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 291 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 289 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 333 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 339 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 339 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 331 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 265 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 324 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 340 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 340 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 338 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 123 exists
                 if ( aNeighbors[ 123 ] != nullptr )
                 {
                     // get neighbor 1 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 285 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 291 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 333 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 261 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 267 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 283 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 292 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 292 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 290 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 326 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 334 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 334 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 340 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 340 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 332 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 259 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 268 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 266 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 325 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 327 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 341 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 341 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != nullptr && aNeighbors[ 339 ] == nullptr )
                     {
                         // copy pointer into big array
                         aNeighbors[ 339 ] = tNeighbor;
                     }
                 }
             } // end order 3
             } // end order 2
         }
// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
#endif /* SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_ */

