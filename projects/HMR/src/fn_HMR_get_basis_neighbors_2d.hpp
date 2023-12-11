/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_HMR_get_basis_neighbors_2d.hpp
 *
 */

#ifndef SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_
#define SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_

#include "cl_HMR_Basis.hpp"
#include "typedefs.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------

    void get_basis_neighbors_2d(       Basis_Function  * aBasis,
                                 uint aOrder,
                                       Basis_Function ** aNeighbors )
    {
         // make sure order is not too big
         MORIS_ASSERT( 0 < aOrder && aOrder <= 3, "Neighbor order too big.");

         // array that contains max size
         uint tArraySize[ 4 ] = { 0, 8, 24, 48 };

         if ( aOrder >= 2 )
         {
             // fill first frame
             for ( uint k=0; k<8; ++k)
             {
                 aNeighbors[ k ] = aBasis->get_neighbor( k );
             }

             // initialize higher order neighbors with null
             for( uint k=8; k<tArraySize[ aOrder ]; ++k )
             {
                 aNeighbors[ k ] = nullptr;
             }

             // placeholder for basis neighbor
             Basis_Function* tNeighbor;

             // test if neighbor 0 exists
             if ( aNeighbors[ 0 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 0
                 tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 10 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 10 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 0
                 tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 9 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 9 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 0
                 tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 11 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 11 ] = tNeighbor;
                 }
             }

             // test if neighbor 1 exists
             if ( aNeighbors[ 1 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 1
                 tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 16 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 16 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 1
                 tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 14 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 14 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 1
                 tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 18 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 18 ] = tNeighbor;
                 }
             }

             // test if neighbor 2 exists
             if ( aNeighbors[ 2 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 2
                 tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 21 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 21 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 2
                 tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 22 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 22 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 2
                 tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 20 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 20 ] = tNeighbor;
                 }
             }

             // test if neighbor 3 exists
             if ( aNeighbors[ 3 ] != nullptr )
             {
                 // get neighbor 3 of neighbor 3
                 tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 15 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 15 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 3
                 tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 13 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 13 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 3
                 tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 17 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 17 ] = tNeighbor;
                 }
             }

             // test if neighbor 4 exists
             if ( aNeighbors[ 4 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 4
                 tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 9 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 9 ] = tNeighbor;
                 }
                 // get neighbor 3 of neighbor 4
                 tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 13 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 13 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 4
                 tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 8 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 8 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 4
                 tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 10 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 10 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 4
                 tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 15 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 15 ] = tNeighbor;
                 }
             }

             // test if neighbor 5 exists
             if ( aNeighbors[ 5 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 5
                 tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 11 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 11 ] = tNeighbor;
                 }
                 // get neighbor 1 of neighbor 5
                 tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 14 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 14 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 5
                 tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 10 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 10 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 5
                 tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 12 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 12 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 5
                 tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 16 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 16 ] = tNeighbor;
                 }
             }

             // test if neighbor 6 exists
             if ( aNeighbors[ 6 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 6
                 tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 18 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 18 ] = tNeighbor;
                 }
                 // get neighbor 2 of neighbor 6
                 tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 22 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 22 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 6
                 tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 16 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 16 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 6
                 tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 23 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 23 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 6
                 tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 21 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 21 ] = tNeighbor;
                 }
             }

             // test if neighbor 7 exists
             if ( aNeighbors[ 7 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 7
                 tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 20 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 20 ] = tNeighbor;
                 }
                 // get neighbor 3 of neighbor 7
                 tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 17 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 17 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 7
                 tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 15 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 15 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 7
                 tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 21 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 21 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 7
                 tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 19 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 19 ] = tNeighbor;
                 }
             }

         if ( aOrder >= 3 )
         {
             // test if neighbor 8 exists
             if ( aNeighbors[ 8 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 8
                 tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 25 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 25 ] = tNeighbor;
                 }
                 // get neighbor 3 of neighbor 8
                 tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 31 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 31 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 8
                 tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 24 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 24 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 8
                 tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 26 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 26 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 8
                 tNeighbor =  aNeighbors[ 8 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 33 ] = tNeighbor;
                 }
             }

             // test if neighbor 9 exists
             if ( aNeighbors[ 9 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 9
                 tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 26 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 26 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 9
                 tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 25 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 25 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 9
                 tNeighbor =  aNeighbors[ 9 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 27 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 27 ] = tNeighbor;
                 }
             }

             // test if neighbor 10 exists
             if ( aNeighbors[ 10 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 10
                 tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 27 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 27 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 10
                 tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 26 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 26 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 10
                 tNeighbor =  aNeighbors[ 10 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 28 ] = tNeighbor;
                 }
             }

             // test if neighbor 11 exists
             if ( aNeighbors[ 11 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 11
                 tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 28 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 11
                 tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 27 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 27 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 11
                 tNeighbor =  aNeighbors[ 11 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 29 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 29 ] = tNeighbor;
                 }
             }

             // test if neighbor 12 exists
             if ( aNeighbors[ 12 ] != nullptr )
             {
                 // get neighbor 0 of neighbor 12
                 tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 0 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 29 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 29 ] = tNeighbor;
                 }
                 // get neighbor 1 of neighbor 12
                 tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 32 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 12
                 tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 28 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 28 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 12
                 tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 30 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 30 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 12
                 tNeighbor =  aNeighbors[ 12 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 34 ] = tNeighbor;
                 }
             }

             // test if neighbor 13 exists
             if ( aNeighbors[ 13 ] != nullptr )
             {
                 // get neighbor 3 of neighbor 13
                 tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 33 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 13
                 tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 31 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 31 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 13
                 tNeighbor =  aNeighbors[ 13 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 35 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 35 ] = tNeighbor;
                 }
             }

             // test if neighbor 14 exists
             if ( aNeighbors[ 14 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 14
                 tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 34 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 14
                 tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 32 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 32 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 14
                 tNeighbor =  aNeighbors[ 14 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 36 ] = tNeighbor;
                 }
             }

             // test if neighbor 15 exists
             if ( aNeighbors[ 15 ] != nullptr )
             {
                 // get neighbor 3 of neighbor 15
                 tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 35 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 35 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 15
                 tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 33 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 33 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 15
                 tNeighbor =  aNeighbors[ 15 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 37 ] = tNeighbor;
                 }
             }

             // test if neighbor 16 exists
             if ( aNeighbors[ 16 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 16
                 tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 36 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 16
                 tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 34 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 34 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 16
                 tNeighbor =  aNeighbors[ 16 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 38 ] = tNeighbor;
                 }
             }

             // test if neighbor 17 exists
             if ( aNeighbors[ 17 ] != nullptr )
             {
                 // get neighbor 3 of neighbor 17
                 tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 37 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 17
                 tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 35 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 35 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 17
                 tNeighbor =  aNeighbors[ 17 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 39 ] = tNeighbor;
                 }
             }

             // test if neighbor 18 exists
             if ( aNeighbors[ 18 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 18
                 tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 38 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 18
                 tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 36 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 36 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 18
                 tNeighbor =  aNeighbors[ 18 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 40 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 40 ] = tNeighbor;
                 }
             }

             // test if neighbor 19 exists
             if ( aNeighbors[ 19 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 19
                 tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 42 ] = tNeighbor;
                 }
                 // get neighbor 3 of neighbor 19
                 tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 3 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 39 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 39 ] = tNeighbor;
                 }
                 // get neighbor 4 of neighbor 19
                 tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 4 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 37 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 37 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 19
                 tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 43 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 19
                 tNeighbor =  aNeighbors[ 19 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 41 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 41 ] = tNeighbor;
                 }
             }

             // test if neighbor 20 exists
             if ( aNeighbors[ 20 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 20
                 tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 43 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 20
                 tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 44 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 20
                 tNeighbor =  aNeighbors[ 20 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 42 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 42 ] = tNeighbor;
                 }
             }

             // test if neighbor 21 exists
             if ( aNeighbors[ 21 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 21
                 tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 44 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 21
                 tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 45 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 45 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 21
                 tNeighbor =  aNeighbors[ 21 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 43 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 43 ] = tNeighbor;
                 }
             }

             // test if neighbor 22 exists
             if ( aNeighbors[ 22 ] != nullptr )
             {
                 // get neighbor 2 of neighbor 22
                 tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 45 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 45 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 22
                 tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 46 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 46 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 22
                 tNeighbor =  aNeighbors[ 22 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 44 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 44 ] = tNeighbor;
                 }
             }

             // test if neighbor 23 exists
             if ( aNeighbors[ 23 ] != nullptr )
             {
                 // get neighbor 1 of neighbor 23
                 tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 1 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 40 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 40 ] = tNeighbor;
                 }
                 // get neighbor 2 of neighbor 23
                 tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 2 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 46 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 46 ] = tNeighbor;
                 }
                 // get neighbor 5 of neighbor 23
                 tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 5 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 38 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 38 ] = tNeighbor;
                 }
                 // get neighbor 6 of neighbor 23
                 tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 6 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 47 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 47 ] = tNeighbor;
                 }
                 // get neighbor 7 of neighbor 23
                 tNeighbor =  aNeighbors[ 23 ]->get_neighbor( 7 );

                 // test if neighbor exists and was not copied yet
                 if ( tNeighbor != nullptr && aNeighbors[ 45 ] == nullptr )
                 {
                     // copy pointer into big array
                     aNeighbors[ 45 ] = tNeighbor;
                 }
             }
         } // end order 3
         } // end order 2
     }
// ----------------------------------------------------------------------------
} /* namespace moris */
#endif /* SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_ */

