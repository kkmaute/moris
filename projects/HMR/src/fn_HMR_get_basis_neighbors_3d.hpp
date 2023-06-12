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
                 if ( aNeighbors[ 0 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 52 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 54 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 84 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 86 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 1 exists
                 if ( aNeighbors[ 1 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 57 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 61 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 89 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 93 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 2 exists
                 if ( aNeighbors[ 2 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 65 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 63 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 97 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 95 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 3 exists
                 if ( aNeighbors[ 3 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 56 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 60 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 88 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 92 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 4 exists
                 if ( aNeighbors[ 4 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 32 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 34 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 44 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 42 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 5 exists
                 if ( aNeighbors[ 5 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 105 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 107 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 117 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 115 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 6 exists
                 if ( aNeighbors[ 6 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 28 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 34 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 32 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 52 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 54 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 27 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 27 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 29 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 29 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 7 exists
                 if ( aNeighbors[ 7 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 34 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 40 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 44 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 57 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 61 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 35 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 35 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 45 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 45 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 8 exists
                 if ( aNeighbors[ 8 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 44 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 48 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 42 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 65 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 63 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 49 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 49 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 47 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 47 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 8
                     tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 9 exists
                 if ( aNeighbors[ 9 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 32 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 42 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 36 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 56 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 60 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 31 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 31 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 41 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 41 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 9
                     tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 10 exists
                 if ( aNeighbors[ 10 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 52 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 56 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 67 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 84 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 88 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 51 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 51 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 83 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 83 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 10
                     tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 11 exists
                 if ( aNeighbors[ 11 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 54 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 57 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 71 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 86 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 89 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 55 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 55 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 87 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 87 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 11
                     tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 12 exists
                 if ( aNeighbors[ 12 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 61 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 65 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 82 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 93 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 97 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 66 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 66 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 98 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 98 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 12
                     tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 13 exists
                 if ( aNeighbors[ 13 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 63 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 60 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 78 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 95 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 92 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 62 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 62 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 13
                     tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 94 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 94 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 14 exists
                 if ( aNeighbors[ 14 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 84 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 86 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 101 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 107 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 105 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 100 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 100 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 102 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 102 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 14
                     tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 15 exists
                 if ( aNeighbors[ 15 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 89 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 93 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 107 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 113 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 117 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 108 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 108 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 118 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 118 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 15
                     tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 16 exists
                 if ( aNeighbors[ 16 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 97 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 95 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 117 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 121 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 115 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 122 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 122 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 16
                     tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 120 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 120 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 17 exists
                 if ( aNeighbors[ 17 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 88 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 92 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 105 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 115 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 109 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 104 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 104 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 17
                     tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 114 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 114 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 18 exists
                 if ( aNeighbors[ 18 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 52 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 52 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 56 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 56 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 32 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 32 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 27 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 27 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 31 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 31 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 51 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 51 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 26 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 26 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 28 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 36 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 67 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 18
                     tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 19 exists
                 if ( aNeighbors[ 19 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 54 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 54 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 57 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 57 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 34 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 34 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 29 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 29 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 35 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 35 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 33 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 33 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 53 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 53 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 55 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 55 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 28 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 28 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 30 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 30 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 40 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 71 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 19
                     tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 20 exists
                 if ( aNeighbors[ 20 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 61 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 61 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 65 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 65 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 44 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 44 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 39 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 39 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 45 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 45 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 49 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 49 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 59 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 59 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 66 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 66 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 40 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 40 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 50 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 50 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 48 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 82 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 20
                     tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 21 exists
                 if ( aNeighbors[ 21 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 63 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 63 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 60 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 60 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 42 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 42 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 37 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 37 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 43 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 43 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 47 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 47 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 41 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 41 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 58 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 58 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 64 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 64 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 62 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 62 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 36 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 36 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 38 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 38 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 48 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 48 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 46 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 46 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 21
                     tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 78 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 22 exists
                 if ( aNeighbors[ 22 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 84 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 84 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 88 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 88 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 105 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 105 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 68 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 68 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 72 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 72 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 83 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 83 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 100 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 100 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 104 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 104 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 67 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 67 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 99 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 99 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 101 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 22
                     tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 109 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 23 exists
                 if ( aNeighbors[ 23 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 86 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 86 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 89 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 89 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 107 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 107 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 70 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 70 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 73 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 73 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 85 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 85 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 87 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 87 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 102 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 102 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 108 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 108 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 106 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 106 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 69 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 69 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 71 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 71 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 101 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 101 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 103 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 103 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 113 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 23
                     tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 24 exists
                 if ( aNeighbors[ 24 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 93 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 93 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 97 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 97 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 117 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 117 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 77 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 77 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 81 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 81 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 91 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 91 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 98 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 98 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 112 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 112 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 118 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 118 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 122 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 122 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 75 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 75 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 82 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 82 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 113 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 113 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 123 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 123 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 24
                     tNeighbor =  aNeighbors[ 24 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 121 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 25 exists
                 if ( aNeighbors[ 25 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 95 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 95 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 92 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 92 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 115 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 115 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 79 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 79 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 76 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 76 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 90 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 90 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 96 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 96 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 94 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 94 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 110 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 110 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 116 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 116 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 120 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 120 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 114 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 114 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 74 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 74 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 80 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 80 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 78 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 78 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 109 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 109 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 111 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 111 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 121 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 121 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 25
                     tNeighbor =  aNeighbors[ 25 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 119 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 119 ] = tNeighbor;
                     }
                 }
             if ( aOrder >= 3 )
             {

                 // test if neighbor 26 exists
                 if ( aNeighbors[ 26 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 174 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 180 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 132 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 125 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 125 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 131 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 131 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 173 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 173 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 124 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 124 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 126 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 138 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 197 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 26
                     tNeighbor =  aNeighbors[ 26 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 27 exists
                 if ( aNeighbors[ 27 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 126 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 132 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 174 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 125 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 125 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 127 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 27
                     tNeighbor =  aNeighbors[ 27 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 28 exists
                 if ( aNeighbors[ 28 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 127 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 126 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 126 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 128 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 28
                     tNeighbor =  aNeighbors[ 28 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 29 exists
                 if ( aNeighbors[ 29 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 128 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 136 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 178 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 127 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 127 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 129 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 129 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 29
                     tNeighbor =  aNeighbors[ 29 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 30 exists
                 if ( aNeighbors[ 30 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 178 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 181 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 136 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 129 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 129 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 137 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 137 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 179 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 179 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 128 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 128 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 130 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 130 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 144 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 203 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 30
                     tNeighbor =  aNeighbors[ 30 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 31 exists
                 if ( aNeighbors[ 31 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 132 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 138 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 180 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 131 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 131 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 145 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 31
                     tNeighbor =  aNeighbors[ 31 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 32 exists
                 if ( aNeighbors[ 32 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 132 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 132 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 32
                     tNeighbor =  aNeighbors[ 32 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 33 exists
                 if ( aNeighbors[ 33 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 133 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 133 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 33
                     tNeighbor =  aNeighbors[ 33 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 34 exists
                 if ( aNeighbors[ 34 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 134 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 134 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 136 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 34
                     tNeighbor =  aNeighbors[ 34 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 35 exists
                 if ( aNeighbors[ 35 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 136 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 136 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 144 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 181 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 135 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 135 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 137 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 137 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 151 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 35
                     tNeighbor =  aNeighbors[ 35 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 36 exists
                 if ( aNeighbors[ 36 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 145 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 138 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 138 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 152 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 36
                     tNeighbor =  aNeighbors[ 36 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 37 exists
                 if ( aNeighbors[ 37 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 139 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 139 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 37
                     tNeighbor =  aNeighbors[ 37 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 38 exists
                 if ( aNeighbors[ 38 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 140 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 140 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 38
                     tNeighbor =  aNeighbors[ 38 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 39 exists
                 if ( aNeighbors[ 39 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 141 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 141 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 39
                     tNeighbor =  aNeighbors[ 39 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 40 exists
                 if ( aNeighbors[ 40 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 143 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 143 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 151 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 142 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 142 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 144 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 144 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 158 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 40
                     tNeighbor =  aNeighbors[ 40 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 41 exists
                 if ( aNeighbors[ 41 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 160 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 152 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 188 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 145 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 145 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 159 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 159 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 41
                     tNeighbor =  aNeighbors[ 41 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 42 exists
                 if ( aNeighbors[ 42 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 146 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 146 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 42
                     tNeighbor =  aNeighbors[ 42 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 160 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 43 exists
                 if ( aNeighbors[ 43 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 147 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 147 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 43
                     tNeighbor =  aNeighbors[ 43 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 44 exists
                 if ( aNeighbors[ 44 ] != NULL )
                 {
                     // get neighbor 4 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 148 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 148 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 164 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 44
                     tNeighbor =  aNeighbors[ 44 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 45 exists
                 if ( aNeighbors[ 45 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 150 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 150 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 158 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 164 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 189 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 149 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 149 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 151 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 151 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 165 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 165 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 45
                     tNeighbor =  aNeighbors[ 45 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 46 exists
                 if ( aNeighbors[ 46 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 191 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 188 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 160 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 167 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 167 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 159 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 159 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 190 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 190 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 152 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 152 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 168 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 166 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 166 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 46
                     tNeighbor =  aNeighbors[ 46 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 214 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 47 exists
                 if ( aNeighbors[ 47 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 168 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 160 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 160 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 191 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 153 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 153 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 169 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 167 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 167 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 47
                     tNeighbor =  aNeighbors[ 47 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 48 exists
                 if ( aNeighbors[ 48 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 169 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 161 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 161 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 154 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 154 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 170 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 168 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 168 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 48
                     tNeighbor =  aNeighbors[ 48 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 49 exists
                 if ( aNeighbors[ 49 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 164 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 170 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 162 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 162 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 195 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 155 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 155 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 171 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 171 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 169 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 169 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 49
                     tNeighbor =  aNeighbors[ 49 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 50 exists
                 if ( aNeighbors[ 50 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 189 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 195 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 164 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 164 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 157 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 157 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 165 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 165 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 171 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 171 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 163 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 163 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 196 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 196 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 156 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 156 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 158 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 158 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 172 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 172 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 170 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 170 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 220 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 50
                     tNeighbor =  aNeighbors[ 50 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 51 exists
                 if ( aNeighbors[ 51 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 174 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 180 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 197 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 173 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 173 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 221 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 51
                     tNeighbor =  aNeighbors[ 51 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 52 exists
                 if ( aNeighbors[ 52 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 174 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 174 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 52
                     tNeighbor =  aNeighbors[ 52 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 53 exists
                 if ( aNeighbors[ 53 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 175 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 175 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 53
                     tNeighbor =  aNeighbors[ 53 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 54 exists
                 if ( aNeighbors[ 54 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 176 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 176 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 178 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 54
                     tNeighbor =  aNeighbors[ 54 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 55 exists
                 if ( aNeighbors[ 55 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 178 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 178 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 181 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 203 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 177 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 177 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 179 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 179 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 227 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 55
                     tNeighbor =  aNeighbors[ 55 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 56 exists
                 if ( aNeighbors[ 56 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 180 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 180 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 56
                     tNeighbor =  aNeighbors[ 56 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 57 exists
                 if ( aNeighbors[ 57 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 181 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 181 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 57
                     tNeighbor =  aNeighbors[ 57 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 58 exists
                 if ( aNeighbors[ 58 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 182 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 182 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 58
                     tNeighbor =  aNeighbors[ 58 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 59 exists
                 if ( aNeighbors[ 59 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 183 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 183 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 59
                     tNeighbor =  aNeighbors[ 59 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 60 exists
                 if ( aNeighbors[ 60 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 184 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 184 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 188 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 60
                     tNeighbor =  aNeighbors[ 60 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 61 exists
                 if ( aNeighbors[ 61 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 185 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 185 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 189 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 61
                     tNeighbor =  aNeighbors[ 61 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 62 exists
                 if ( aNeighbors[ 62 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 191 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 188 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 188 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 214 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 186 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 186 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 190 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 190 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 62
                     tNeighbor =  aNeighbors[ 62 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 238 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 63 exists
                 if ( aNeighbors[ 63 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 191 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 191 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 63
                     tNeighbor =  aNeighbors[ 63 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 64 exists
                 if ( aNeighbors[ 64 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 192 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 192 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 64
                     tNeighbor =  aNeighbors[ 64 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 65 exists
                 if ( aNeighbors[ 65 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 195 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 193 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 193 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 65
                     tNeighbor =  aNeighbors[ 65 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 66 exists
                 if ( aNeighbors[ 66 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 189 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 189 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 195 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 195 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 220 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 187 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 187 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 196 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 196 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 194 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 194 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 244 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 66
                     tNeighbor =  aNeighbors[ 66 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 67 exists
                 if ( aNeighbors[ 67 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 221 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 197 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 197 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 245 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 67
                     tNeighbor =  aNeighbors[ 67 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 68 exists
                 if ( aNeighbors[ 68 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 198 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 198 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 68
                     tNeighbor =  aNeighbors[ 68 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 69 exists
                 if ( aNeighbors[ 69 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 199 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 199 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 69
                     tNeighbor =  aNeighbors[ 69 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 70 exists
                 if ( aNeighbors[ 70 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 200 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 200 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 70
                     tNeighbor =  aNeighbors[ 70 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 71 exists
                 if ( aNeighbors[ 71 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 202 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 202 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 227 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 201 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 201 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 203 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 203 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 251 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 71
                     tNeighbor =  aNeighbors[ 71 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 72 exists
                 if ( aNeighbors[ 72 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 204 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 204 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 72
                     tNeighbor =  aNeighbors[ 72 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 73 exists
                 if ( aNeighbors[ 73 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 205 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 205 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 73
                     tNeighbor =  aNeighbors[ 73 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 74 exists
                 if ( aNeighbors[ 74 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 206 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 206 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 74
                     tNeighbor =  aNeighbors[ 74 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 75 exists
                 if ( aNeighbors[ 75 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 207 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 207 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 75
                     tNeighbor =  aNeighbors[ 75 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 76 exists
                 if ( aNeighbors[ 76 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 208 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 208 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 76
                     tNeighbor =  aNeighbors[ 76 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 77 exists
                 if ( aNeighbors[ 77 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 209 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 209 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 77
                     tNeighbor =  aNeighbors[ 77 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 78 exists
                 if ( aNeighbors[ 78 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 212 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 212 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 238 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 210 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 210 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 214 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 214 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 78
                     tNeighbor =  aNeighbors[ 78 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 262 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 79 exists
                 if ( aNeighbors[ 79 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 215 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 215 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 79
                     tNeighbor =  aNeighbors[ 79 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 80 exists
                 if ( aNeighbors[ 80 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 216 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 216 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 80
                     tNeighbor =  aNeighbors[ 80 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 81 exists
                 if ( aNeighbors[ 81 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 217 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 217 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 81
                     tNeighbor =  aNeighbors[ 81 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 82 exists
                 if ( aNeighbors[ 82 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 213 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 213 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 219 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 219 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 244 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 211 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 211 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 220 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 220 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 218 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 218 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 268 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 82
                     tNeighbor =  aNeighbors[ 82 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 83 exists
                 if ( aNeighbors[ 83 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 245 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 270 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 276 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 221 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 221 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 269 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 269 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 83
                     tNeighbor =  aNeighbors[ 83 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 84 exists
                 if ( aNeighbors[ 84 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 222 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 222 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 270 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 84
                     tNeighbor =  aNeighbors[ 84 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 85 exists
                 if ( aNeighbors[ 85 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 223 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 223 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 85
                     tNeighbor =  aNeighbors[ 85 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 86 exists
                 if ( aNeighbors[ 86 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 224 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 224 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 86
                     tNeighbor =  aNeighbors[ 86 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 274 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 87 exists
                 if ( aNeighbors[ 87 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 226 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 226 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 251 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 274 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 277 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 225 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 225 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 227 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 227 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 275 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 275 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 87
                     tNeighbor =  aNeighbors[ 87 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 88 exists
                 if ( aNeighbors[ 88 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 228 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 228 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 276 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 88
                     tNeighbor =  aNeighbors[ 88 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 89 exists
                 if ( aNeighbors[ 89 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 229 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 229 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 277 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 89
                     tNeighbor =  aNeighbors[ 89 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 90 exists
                 if ( aNeighbors[ 90 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 230 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 230 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 90
                     tNeighbor =  aNeighbors[ 90 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 91 exists
                 if ( aNeighbors[ 91 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 231 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 231 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 91
                     tNeighbor =  aNeighbors[ 91 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 92 exists
                 if ( aNeighbors[ 92 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 232 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 232 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 92
                     tNeighbor =  aNeighbors[ 92 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 284 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 93 exists
                 if ( aNeighbors[ 93 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 233 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 233 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 93
                     tNeighbor =  aNeighbors[ 93 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 285 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 94 exists
                 if ( aNeighbors[ 94 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 236 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 236 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 262 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 287 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 284 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 234 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 234 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 238 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 238 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 94
                     tNeighbor =  aNeighbors[ 94 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 286 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 286 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 95 exists
                 if ( aNeighbors[ 95 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 239 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 239 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 95
                     tNeighbor =  aNeighbors[ 95 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 287 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 96 exists
                 if ( aNeighbors[ 96 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 240 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 240 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 96
                     tNeighbor =  aNeighbors[ 96 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 97 exists
                 if ( aNeighbors[ 97 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 241 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 241 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 291 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 97
                     tNeighbor =  aNeighbors[ 97 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 98 exists
                 if ( aNeighbors[ 98 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 237 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 237 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 243 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 243 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 268 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 285 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 291 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 235 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 235 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 244 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 244 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 242 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 242 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 292 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 292 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 98
                     tNeighbor =  aNeighbors[ 98 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 99 exists
                 if ( aNeighbors[ 99 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 270 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 276 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 301 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 269 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 269 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 294 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 294 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 300 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 300 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 245 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 245 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 293 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 293 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 295 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 99
                     tNeighbor =  aNeighbors[ 99 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 307 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 100 exists
                 if ( aNeighbors[ 100 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 270 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 270 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 295 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 301 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 246 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 246 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 294 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 294 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 296 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 100
                     tNeighbor =  aNeighbors[ 100 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 101 exists
                 if ( aNeighbors[ 101 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 271 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 271 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 296 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 247 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 247 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 295 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 295 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 297 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 101
                     tNeighbor =  aNeighbors[ 101 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 102 exists
                 if ( aNeighbors[ 102 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 272 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 272 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 274 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 297 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 305 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 248 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 248 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 296 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 296 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 298 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 298 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 102
                     tNeighbor =  aNeighbors[ 102 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 103 exists
                 if ( aNeighbors[ 103 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 274 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 274 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 277 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 305 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 250 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 250 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 273 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 273 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 275 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 275 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 298 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 298 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 306 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 306 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 249 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 249 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 251 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 251 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 297 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 297 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 299 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 299 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 313 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 103
                     tNeighbor =  aNeighbors[ 103 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 104 exists
                 if ( aNeighbors[ 104 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 276 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 276 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 301 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 307 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 252 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 252 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 300 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 300 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 104
                     tNeighbor =  aNeighbors[ 104 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 314 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 105 exists
                 if ( aNeighbors[ 105 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 301 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 301 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 105
                     tNeighbor =  aNeighbors[ 105 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 106 exists
                 if ( aNeighbors[ 106 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 302 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 302 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 106
                     tNeighbor =  aNeighbors[ 106 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 107 exists
                 if ( aNeighbors[ 107 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 303 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 303 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 305 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 107
                     tNeighbor =  aNeighbors[ 107 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 108 exists
                 if ( aNeighbors[ 108 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 277 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 277 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 305 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 305 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 313 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 253 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 253 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 304 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 304 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 306 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 306 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 320 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 108
                     tNeighbor =  aNeighbors[ 108 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 109 exists
                 if ( aNeighbors[ 109 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 278 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 278 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 314 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 254 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 254 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 307 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 307 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 109
                     tNeighbor =  aNeighbors[ 109 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 321 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 110 exists
                 if ( aNeighbors[ 110 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 308 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 308 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 110
                     tNeighbor =  aNeighbors[ 110 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 111 exists
                 if ( aNeighbors[ 111 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 309 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 309 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 111
                     tNeighbor =  aNeighbors[ 111 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 112 exists
                 if ( aNeighbors[ 112 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 310 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 310 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 112
                     tNeighbor =  aNeighbors[ 112 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 113 exists
                 if ( aNeighbors[ 113 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 279 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 279 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 312 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 312 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 320 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 255 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 255 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 311 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 311 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 313 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 313 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 327 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 113
                     tNeighbor =  aNeighbors[ 113 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 114 exists
                 if ( aNeighbors[ 114 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 280 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 280 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 284 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 329 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 321 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 256 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 256 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 314 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 314 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 114
                     tNeighbor =  aNeighbors[ 114 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 328 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 328 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 115 exists
                 if ( aNeighbors[ 115 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 315 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 315 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 115
                     tNeighbor =  aNeighbors[ 115 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 329 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 116 exists
                 if ( aNeighbors[ 116 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 316 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 316 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 116
                     tNeighbor =  aNeighbors[ 116 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 117 exists
                 if ( aNeighbors[ 117 ] != NULL )
                 {
                     // get neighbor 5 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 317 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 317 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 333 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 117
                     tNeighbor =  aNeighbors[ 117 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 118 exists
                 if ( aNeighbors[ 118 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 281 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 281 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 285 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 319 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 319 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 327 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 333 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 257 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 257 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 318 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 318 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 320 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 320 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 334 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 334 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 118
                     tNeighbor =  aNeighbors[ 118 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 119 exists
                 if ( aNeighbors[ 119 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 287 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 284 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 284 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 329 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 9 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 9 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 260 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 260 ] = tNeighbor;
                     }
                     // get neighbor 10 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 10 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 282 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 282 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 286 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 286 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 336 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 336 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 328 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 328 ] = tNeighbor;
                     }
                     // get neighbor 18 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 18 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 258 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 258 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 262 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 262 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 321 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 321 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 337 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 119
                     tNeighbor =  aNeighbors[ 119 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 335 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 335 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 120 exists
                 if ( aNeighbors[ 120 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 287 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 287 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 337 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 329 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 329 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 263 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 263 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 322 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 322 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 338 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 120
                     tNeighbor =  aNeighbors[ 120 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 336 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 336 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 121 exists
                 if ( aNeighbors[ 121 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 288 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 288 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 338 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 330 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 330 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 264 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 264 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 323 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 323 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 339 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 339 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 121
                     tNeighbor =  aNeighbors[ 121 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 337 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 337 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 122 exists
                 if ( aNeighbors[ 122 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 291 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 289 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 289 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 333 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 339 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 339 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 331 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 331 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 265 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 265 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 324 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 324 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 340 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 340 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 122
                     tNeighbor =  aNeighbors[ 122 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 338 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 338 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 123 exists
                 if ( aNeighbors[ 123 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 285 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 285 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 291 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 291 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 333 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 333 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 261 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 261 ] = tNeighbor;
                     }
                     // get neighbor 8 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 8 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 267 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 267 ] = tNeighbor;
                     }
                     // get neighbor 11 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 11 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 283 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 283 ] = tNeighbor;
                     }
                     // get neighbor 12 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 12 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 292 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 292 ] = tNeighbor;
                     }
                     // get neighbor 13 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 13 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 290 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 290 ] = tNeighbor;
                     }
                     // get neighbor 14 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 14 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 326 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 326 ] = tNeighbor;
                     }
                     // get neighbor 15 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 15 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 334 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 334 ] = tNeighbor;
                     }
                     // get neighbor 16 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 16 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 340 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 340 ] = tNeighbor;
                     }
                     // get neighbor 17 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 17 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 332 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 332 ] = tNeighbor;
                     }
                     // get neighbor 19 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 19 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 259 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 259 ] = tNeighbor;
                     }
                     // get neighbor 20 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 20 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 268 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 268 ] = tNeighbor;
                     }
                     // get neighbor 21 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 21 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 266 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 266 ] = tNeighbor;
                     }
                     // get neighbor 22 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 22 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 325 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 325 ] = tNeighbor;
                     }
                     // get neighbor 23 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 23 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 327 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 327 ] = tNeighbor;
                     }
                     // get neighbor 24 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 24 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 341 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 341 ] = tNeighbor;
                     }
                     // get neighbor 25 of neighbor 123
                     tNeighbor =  aNeighbors[ 123 ]->get_neighbor( 25 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 339 ] == NULL )
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

