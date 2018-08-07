/*
 * fn_HMR_get_basis_neighbors_2d.hpp
 *
 *  Created on: July 26, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_
#define SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_

#include "typedefs.hpp" //COR/src
#include "cl_HMR_Basis.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        void
        get_basis_neighbors_2d(
            Basis          * aBasis,
            const uint     & aOrder,
            Basis          ** aNeighbors )
        {
             // make sure order is not too big
             MORIS_ASSERT( 0 < aOrder && aOrder <= 2, "Neighbor order too big.");

             // array that contains max size
             uint tArraySize[ 3 ] = { 0, 8, 24 };

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
                 Basis* tNeighbor;

                 // test if neighbor 0 exists
                 if ( aNeighbors[ 0 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 10 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 10 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 9 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 9 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 0
                     tNeighbor =  aNeighbors[ 0 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 11 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 11 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 1 exists
                 if ( aNeighbors[ 1 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 16 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 16 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 14 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 14 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 1
                     tNeighbor =  aNeighbors[ 1 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 18 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 18 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 2 exists
                 if ( aNeighbors[ 2 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 21 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 21 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 22 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 22 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 2
                     tNeighbor =  aNeighbors[ 2 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 20 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 20 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 3 exists
                 if ( aNeighbors[ 3 ] != NULL )
                 {
                     // get neighbor 3 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 15 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 15 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 13 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 13 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 3
                     tNeighbor =  aNeighbors[ 3 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 17 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 17 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 4 exists
                 if ( aNeighbors[ 4 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 9 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 9 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 13 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 13 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 8 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 8 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 10 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 10 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 4
                     tNeighbor =  aNeighbors[ 4 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 15 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 15 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 5 exists
                 if ( aNeighbors[ 5 ] != NULL )
                 {
                     // get neighbor 0 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 11 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 11 ] = tNeighbor;
                     }
                     // get neighbor 1 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 14 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 14 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 10 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 10 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 12 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 12 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 5
                     tNeighbor =  aNeighbors[ 5 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 16 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 16 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 6 exists
                 if ( aNeighbors[ 6 ] != NULL )
                 {
                     // get neighbor 1 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 18 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 18 ] = tNeighbor;
                     }
                     // get neighbor 2 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 22 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 22 ] = tNeighbor;
                     }
                     // get neighbor 5 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 16 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 16 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 23 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 23 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 6
                     tNeighbor =  aNeighbors[ 6 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 21 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 21 ] = tNeighbor;
                     }
                 }

                 // test if neighbor 7 exists
                 if ( aNeighbors[ 7 ] != NULL )
                 {
                     // get neighbor 2 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 20 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 20 ] = tNeighbor;
                     }
                     // get neighbor 3 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 17 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 17 ] = tNeighbor;
                     }
                     // get neighbor 4 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 15 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 15 ] = tNeighbor;
                     }
                     // get neighbor 6 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 6 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 21 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 21 ] = tNeighbor;
                     }
                     // get neighbor 7 of neighbor 7
                     tNeighbor =  aNeighbors[ 7 ]->get_neighbor( 7 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && aNeighbors[ 19 ] == NULL )
                     {
                         // copy pointer into big array
                         aNeighbors[ 19 ] = tNeighbor;
                     }
                 }
             } // end order 2
         }
// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
#endif /* SRC_HMR_FN_HMR_GET_BASIS_NEIGHBORS_2D_HPP_ */
