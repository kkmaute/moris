/*
 * fn_HMR_Background_Element_Neighbors_3D.hpp
 *
 *  Created on: July 26, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_3D_HPP_
#define SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_3D_HPP_

#include "typedefs.hpp" //COR/src
#include "cl_Cell.hpp" //CON/src
#include "cl_HMR_Background_Element_Base.hpp" //HMR/src
#include "cl_HMR_Background_Element.hpp" //HMR/src

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        template <>
        void
        Background_Element< 3, 8, 26 >::get_neighbors_from_same_level(
            const uint                        & aOrder,
            Cell< Background_Element_Base * > & aNeighbors )
        {
             // make sure order is not too big
             MORIS_ASSERT( 0 < aOrder && aOrder <= 2, "Neighbor order too big.");

             // array that contains max size
             uint tArraySize[ 3 ] = { 0, 26, 124 };

             // initialize temporary neighbor array
             Cell< Background_Element_Base * >
             tNeighbors( tArraySize[ aOrder ], nullptr );

             // fill first frame
             for ( uint k=0; k<26; ++k)
             {
                 // test if neighbor exists
                 if ( mNeighbors[ k ] != NULL )
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
                 if ( tNeighbors( 0 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 0
                     tNeighbor =  tNeighbors( 0 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 52 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 54 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 84 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 86 ) == NULL )
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
                 if ( tNeighbors( 1 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 1
                     tNeighbor =  tNeighbors( 1 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 57 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 61 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 89 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 93 ) == NULL )
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
                 if ( tNeighbors( 2 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 2
                     tNeighbor =  tNeighbors( 2 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 65 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 63 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 97 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 95 ) == NULL )
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
                 if ( tNeighbors( 3 ) != NULL )
                 {
                     // get neighbor 3 of neighbor 3
                     tNeighbor =  tNeighbors( 3 )->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 56 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 60 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 88 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 92 ) == NULL )
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
                 if ( tNeighbors( 4 ) != NULL )
                 {
                     // get neighbor 4 of neighbor 4
                     tNeighbor =  tNeighbors( 4 )->get_neighbor( 4 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 32 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 34 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 44 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 42 ) == NULL )
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
                 if ( tNeighbors( 5 ) != NULL )
                 {
                     // get neighbor 5 of neighbor 5
                     tNeighbor =  tNeighbors( 5 )->get_neighbor( 5 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 105 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 107 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 117 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 115 ) == NULL )
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
                 if ( tNeighbors( 6 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 6
                     tNeighbor =  tNeighbors( 6 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 28 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 34 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 32 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 52 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 54 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 27 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 29 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                 if ( tNeighbors( 7 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 7
                     tNeighbor =  tNeighbors( 7 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 34 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 40 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 44 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 57 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 61 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 35 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 45 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                 if ( tNeighbors( 8 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 8
                     tNeighbor =  tNeighbors( 8 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 44 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 48 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 42 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 65 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 63 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 49 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 47 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                 if ( tNeighbors( 9 ) != NULL )
                 {
                     // get neighbor 3 of neighbor 9
                     tNeighbor =  tNeighbors( 9 )->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 32 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 42 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 36 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 56 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 60 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 31 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 41 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                 if ( tNeighbors( 10 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 10
                     tNeighbor =  tNeighbors( 10 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 52 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 56 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 67 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 84 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 88 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 51 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 83 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                 if ( tNeighbors( 11 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 11
                     tNeighbor =  tNeighbors( 11 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 54 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 57 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 71 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 86 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 89 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 55 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 87 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                 if ( tNeighbors( 12 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 12
                     tNeighbor =  tNeighbors( 12 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 61 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 65 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 82 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 93 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 97 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 66 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 98 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                 if ( tNeighbors( 13 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 13
                     tNeighbor =  tNeighbors( 13 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 63 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 60 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 78 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 95 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 92 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 62 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 94 ) == NULL )
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
                 if ( tNeighbors( 14 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 14
                     tNeighbor =  tNeighbors( 14 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 84 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 86 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 101 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 107 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 105 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 100 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 102 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                 if ( tNeighbors( 15 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 15
                     tNeighbor =  tNeighbors( 15 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 89 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 93 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 107 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 113 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 117 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 108 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 118 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                 if ( tNeighbors( 16 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 16
                     tNeighbor =  tNeighbors( 16 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 97 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 95 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 117 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 121 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 115 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 122 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 120 ) == NULL )
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
                 if ( tNeighbors( 17 ) != NULL )
                 {
                     // get neighbor 3 of neighbor 17
                     tNeighbor =  tNeighbors( 17 )->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 88 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 92 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 105 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 115 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 109 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 104 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 114 ) == NULL )
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
                 if ( tNeighbors( 18 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 18
                     tNeighbor =  tNeighbors( 18 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 52 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 56 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 32 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 27 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 31 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 51 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 26 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 28 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 36 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 67 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                 if ( tNeighbors( 19 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 19
                     tNeighbor =  tNeighbors( 19 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 54 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 57 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 34 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 29 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 35 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 33 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 53 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 55 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 28 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 30 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 40 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 71 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                 if ( tNeighbors( 20 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 20
                     tNeighbor =  tNeighbors( 20 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 61 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 65 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 44 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 39 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 45 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 49 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 59 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 66 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 40 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 50 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 48 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 82 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                 if ( tNeighbors( 21 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 21
                     tNeighbor =  tNeighbors( 21 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 63 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 60 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 42 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 37 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 43 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 47 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 41 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 58 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 64 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 62 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 36 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 38 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 48 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 46 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 78 ) == NULL )
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
                 if ( tNeighbors( 22 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 22
                     tNeighbor =  tNeighbors( 22 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 84 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 88 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 105 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 68 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 72 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 83 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 100 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 104 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 67 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 99 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 101 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 109 ) == NULL )
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
                 if ( tNeighbors( 23 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 23
                     tNeighbor =  tNeighbors( 23 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 86 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 89 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 107 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 70 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 73 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 85 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 87 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 102 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 108 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 106 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 69 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 71 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 101 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 103 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 113 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                 if ( tNeighbors( 24 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 24
                     tNeighbor =  tNeighbors( 24 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 93 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 97 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 117 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 77 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 81 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 91 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 98 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 112 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 118 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 122 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 75 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 82 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 113 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 123 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 121 ) == NULL )
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
                 if ( tNeighbors( 25 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 25
                     tNeighbor =  tNeighbors( 25 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 95 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 92 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 115 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 79 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 76 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 90 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 96 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 94 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 110 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 116 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 120 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 114 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 74 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 80 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 78 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 109 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 111 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 121 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 119 ) == NULL )
                     {
                         // test if neighbor is on same level
                         if ( tNeighbor->get_level() == mLevel )
                         {
                             // copy pointer in big array
                             tNeighbors( 119 ) = tNeighbor;
                         }
                     }

                 }
             } // end order 2

             // initialize element counter
             uint tCount = 0;

             // count number of existing elements
             for( auto tNeighbor : tNeighbors )
             {
                 if ( tNeighbor != NULL )
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
                 if ( tNeighbor != NULL ) 
                 {
                     aNeighbors( tCount++ ) = tNeighbor;
                 }
             }
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_3D_HPP_ */