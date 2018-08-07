/*
 * fn_HMR_get_basis_neighbors_3d.hpp
 *
 *  Created on: July 26, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_
#define SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_

#include "typedefs.hpp" //COR/src
#include "cl_HMR_Basis.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        void
        get_basis_neighbors_3d(
            Basis          * aBasis,
            const uint     & aOrder,
            Basis          ** aNeighbors )
        {
             // make sure order is not too big
             MORIS_ASSERT( 0 < aOrder && aOrder <= 2, "Neighbor order too big.");

             // array that contains max size
             uint tArraySize[ 3 ] = { 0, 26, 124 };

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
             } // end order 2
         }
// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
#endif /* SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_3D_HPP_ */
