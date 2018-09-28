/*
 * fn_HMR_Background_Element_Neighbors_2D.hpp
 *
 *  Created on: July 26, 2018
 *  using MATLAB
 */
 
#ifndef SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_
#define SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_

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
        Background_Element< 2, 4, 8, 4, 0 >::get_neighbors_from_same_level(
            const uint                        & aOrder,
            Cell< Background_Element_Base * > & aNeighbors )
        {
             // make sure order is not too big
             MORIS_ASSERT( 0 < aOrder && aOrder <= 2, "Neighbor order too big.");

             // array that contains max size
             uint tArraySize[ 3 ] = { 0, 8, 24 };

             // initialize temporary neighbor array
             Cell< Background_Element_Base * >
             tNeighbors( tArraySize[ aOrder ], nullptr );

             // fill first frame
             for ( uint k=0; k<8; ++k)
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
                     if ( tNeighbor != NULL && tNeighbors( 10 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 9 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 11 ) == NULL )
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
                 if ( tNeighbors( 1 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 1
                     tNeighbor =  tNeighbors( 1 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 16 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 14 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 18 ) == NULL )
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
                 if ( tNeighbors( 2 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 2
                     tNeighbor =  tNeighbors( 2 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 21 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 22 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 20 ) == NULL )
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
                 if ( tNeighbors( 3 ) != NULL )
                 {
                     // get neighbor 3 of neighbor 3
                     tNeighbor =  tNeighbors( 3 )->get_neighbor( 3 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 15 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 13 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 17 ) == NULL )
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
                 if ( tNeighbors( 4 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 4
                     tNeighbor =  tNeighbors( 4 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 9 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 13 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 8 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 10 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 15 ) == NULL )
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
                 if ( tNeighbors( 5 ) != NULL )
                 {
                     // get neighbor 0 of neighbor 5
                     tNeighbor =  tNeighbors( 5 )->get_neighbor( 0 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 11 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 14 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 10 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 12 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 16 ) == NULL )
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
                 if ( tNeighbors( 6 ) != NULL )
                 {
                     // get neighbor 1 of neighbor 6
                     tNeighbor =  tNeighbors( 6 )->get_neighbor( 1 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 18 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 22 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 16 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 23 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 21 ) == NULL )
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
                 if ( tNeighbors( 7 ) != NULL )
                 {
                     // get neighbor 2 of neighbor 7
                     tNeighbor =  tNeighbors( 7 )->get_neighbor( 2 );

                     // test if neighbor exists and was not copied yet
                     if ( tNeighbor != NULL && tNeighbors( 20 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 17 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 15 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 21 ) == NULL )
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
                     if ( tNeighbor != NULL && tNeighbors( 19 ) == NULL )
                     {
                         // test if neighbor is on same level
                         if ( tNeighbor->get_level() == mLevel )
                         {
                             // copy pointer in big array
                             tNeighbors( 19 ) = tNeighbor;
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

#endif /* SRC_HMR_FN_HMR_BACKGROUND_ELEMENT_NEIGHBORS_2D_HPP_ */
