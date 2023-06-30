/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element_Quad16.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD16_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD16_HPP_

#include "cl_HMR_BSpline_Element.hpp"
#include "fn_HMR_get_basis_neighbors_2d.hpp"

namespace moris::hmr
{
    /**
    * Returns the geometry type of this element
    *
    * @return mtk::Geometry_Type
    */
    template<>
    mtk::Geometry_Type
    BSpline_Element< 3, 3, 0 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::QUAD;
    }

// ----------------------------------------------------------------------------

    /**
     * returns the ijk position of a given basis
     *
     * @param[in]  aBasisNumber   element local number of basis
     * @param[out] aIJK           proc local ijk position of this basis
     *
     * @return void
     *
     */
    template<>
    void
    BSpline_Element< 3, 3, 0 >::get_ijk_of_basis(
        uint aBasisNumber,
        luint      * aIJK )
    {
        // get element local coordinate
        switch ( aBasisNumber )
        {
            case(  0 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case(  1 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case(  2 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 3 ;
                break;
            }
            case(  3 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 3 ;
                break;
            }
            case(  4 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case(  5 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case(  6 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case(  7 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case(  8 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 3 ;
                break;
            }
            case(  9 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 3 ;
                break;
            }
            case( 10 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case( 11 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case( 12 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case( 13 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case( 14 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case( 15 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
        }

        // get position of element on background mesh
        const luint * tElIJK = mElement->get_ijk();

        // add element offset
        aIJK[ 0 ] += tElIJK[ 0 ];
        aIJK[ 1 ] += tElIJK[ 1 ];
    }

// ----------------------------------------------------------------------------

    /**
    * Links each basis of an element with neighbor basis.
    *
    * @param[inout] aAllElementsOnProc   cell containing all B-Spline
    *                                    elements including the aura
    * @return void
    */
    template<>
    void
    BSpline_Element< 3, 3, 0 >::link_basis_with_neighbors(
          moris::Cell< Element* > & aAllElementsOnProc )
    {
         // initialize frame of basis around basis from this element
         Basis* tBasis[ 20 ] = { nullptr };

         // get pointer to neighbor  0
         Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

         // test if neighbor  0 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   1 ] = tNeighbor->get_basis(   0 );
             tBasis[   2 ] = tNeighbor->get_basis(   4 );
             tBasis[   3 ] = tNeighbor->get_basis(   5 );
             tBasis[   4 ] = tNeighbor->get_basis(   1 );
         }

         // get pointer to neighbor  1
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

         // test if neighbor  1 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   7 ] = tNeighbor->get_basis(   1 );
             tBasis[   9 ] = tNeighbor->get_basis(   6 );
             tBasis[  11 ] = tNeighbor->get_basis(   7 );
             tBasis[  13 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  2
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

         // test if neighbor  2 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  15 ] = tNeighbor->get_basis(   3 );
             tBasis[  16 ] = tNeighbor->get_basis(   9 );
             tBasis[  17 ] = tNeighbor->get_basis(   8 );
             tBasis[  18 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  3
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

         // test if neighbor  3 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   6 ] = tNeighbor->get_basis(   0 );
             tBasis[   8 ] = tNeighbor->get_basis(  11 );
             tBasis[  10 ] = tNeighbor->get_basis(  10 );
             tBasis[  12 ] = tNeighbor->get_basis(   3 );
         }

         // get pointer to neighbor  4
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

         // test if neighbor  4 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   0 ] = tNeighbor->get_basis(   0 );
             tBasis[   1 ] = tNeighbor->get_basis(   4 );
             tBasis[   2 ] = tNeighbor->get_basis(   5 );
             tBasis[   3 ] = tNeighbor->get_basis(   1 );
             tBasis[   6 ] = tNeighbor->get_basis(  11 );
             tBasis[   8 ] = tNeighbor->get_basis(  10 );
             tBasis[  10 ] = tNeighbor->get_basis(   3 );
         }

         // get pointer to neighbor  5
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

         // test if neighbor  5 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   2 ] = tNeighbor->get_basis(   0 );
             tBasis[   3 ] = tNeighbor->get_basis(   4 );
             tBasis[   4 ] = tNeighbor->get_basis(   5 );
             tBasis[   5 ] = tNeighbor->get_basis(   1 );
             tBasis[   7 ] = tNeighbor->get_basis(   6 );
             tBasis[   9 ] = tNeighbor->get_basis(   7 );
             tBasis[  11 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  6
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

         // test if neighbor  6 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   9 ] = tNeighbor->get_basis(   1 );
             tBasis[  11 ] = tNeighbor->get_basis(   6 );
             tBasis[  13 ] = tNeighbor->get_basis(   7 );
             tBasis[  16 ] = tNeighbor->get_basis(   3 );
             tBasis[  17 ] = tNeighbor->get_basis(   9 );
             tBasis[  18 ] = tNeighbor->get_basis(   8 );
             tBasis[  19 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  7
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

         // test if neighbor  7 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   8 ] = tNeighbor->get_basis(   0 );
             tBasis[  10 ] = tNeighbor->get_basis(  11 );
             tBasis[  12 ] = tNeighbor->get_basis(  10 );
             tBasis[  14 ] = tNeighbor->get_basis(   3 );
             tBasis[  15 ] = tNeighbor->get_basis(   9 );
             tBasis[  16 ] = tNeighbor->get_basis(   8 );
             tBasis[  17 ] = tNeighbor->get_basis(   2 );
         }

         // test if basis 0 exists
         if ( mBasis[   0 ] != nullptr )
         {
             // test if basis 0 has been processed
             if ( ! mBasis[   0 ]->is_flagged() )
             {
                 // link neighbors of basis 0
                 mBasis[   0 ]->insert_neighbor(  0, tBasis[   1 ] );
                 mBasis[   0 ]->insert_neighbor(  1, mBasis[   4 ] );
                 mBasis[   0 ]->insert_neighbor(  2, mBasis[  11 ] );
                 mBasis[   0 ]->insert_neighbor(  3, tBasis[   6 ] );
                 mBasis[   0 ]->insert_neighbor(  4, tBasis[   0 ] );
                 mBasis[   0 ]->insert_neighbor(  5, tBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor(  6, mBasis[  12 ] );
                 mBasis[   0 ]->insert_neighbor(  7, tBasis[   8 ] );

                 // flag this basis
                 mBasis[   0 ]->flag();
             }

         }

         // test if basis 1 exists
         if ( mBasis[   1 ] != nullptr )
         {
             // test if basis 1 has been processed
             if ( ! mBasis[   1 ]->is_flagged() )
             {
                 // link neighbors of basis 1
                 mBasis[   1 ]->insert_neighbor(  0, tBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor(  1, tBasis[   7 ] );
                 mBasis[   1 ]->insert_neighbor(  2, mBasis[   6 ] );
                 mBasis[   1 ]->insert_neighbor(  3, mBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor(  4, tBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor(  5, tBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor(  6, tBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor(  7, mBasis[  13 ] );

                 // flag this basis
                 mBasis[   1 ]->flag();
             }

         }

         // test if basis 2 exists
         if ( mBasis[   2 ] != nullptr )
         {
             // test if basis 2 has been processed
             if ( ! mBasis[   2 ]->is_flagged() )
             {
                 // link neighbors of basis 2
                 mBasis[   2 ]->insert_neighbor(  0, mBasis[   7 ] );
                 mBasis[   2 ]->insert_neighbor(  1, tBasis[  13 ] );
                 mBasis[   2 ]->insert_neighbor(  2, tBasis[  18 ] );
                 mBasis[   2 ]->insert_neighbor(  3, mBasis[   8 ] );
                 mBasis[   2 ]->insert_neighbor(  4, mBasis[  14 ] );
                 mBasis[   2 ]->insert_neighbor(  5, tBasis[  11 ] );
                 mBasis[   2 ]->insert_neighbor(  6, tBasis[  19 ] );
                 mBasis[   2 ]->insert_neighbor(  7, tBasis[  17 ] );

                 // flag this basis
                 mBasis[   2 ]->flag();
             }

         }

         // test if basis 3 exists
         if ( mBasis[   3 ] != nullptr )
         {
             // test if basis 3 has been processed
             if ( ! mBasis[   3 ]->is_flagged() )
             {
                 // link neighbors of basis 3
                 mBasis[   3 ]->insert_neighbor(  0, mBasis[  10 ] );
                 mBasis[   3 ]->insert_neighbor(  1, mBasis[   9 ] );
                 mBasis[   3 ]->insert_neighbor(  2, tBasis[  15 ] );
                 mBasis[   3 ]->insert_neighbor(  3, tBasis[  12 ] );
                 mBasis[   3 ]->insert_neighbor(  4, tBasis[  10 ] );
                 mBasis[   3 ]->insert_neighbor(  5, mBasis[  15 ] );
                 mBasis[   3 ]->insert_neighbor(  6, tBasis[  16 ] );
                 mBasis[   3 ]->insert_neighbor(  7, tBasis[  14 ] );

                 // flag this basis
                 mBasis[   3 ]->flag();
             }

         }

         // test if basis 4 exists
         if ( mBasis[   4 ] != nullptr )
         {
             // test if basis 4 has been processed
             if ( ! mBasis[   4 ]->is_flagged() )
             {
                 // link neighbors of basis 4
                 mBasis[   4 ]->insert_neighbor(  0, tBasis[   2 ] );
                 mBasis[   4 ]->insert_neighbor(  1, mBasis[   5 ] );
                 mBasis[   4 ]->insert_neighbor(  2, mBasis[  12 ] );
                 mBasis[   4 ]->insert_neighbor(  3, mBasis[   0 ] );
                 mBasis[   4 ]->insert_neighbor(  4, tBasis[   1 ] );
                 mBasis[   4 ]->insert_neighbor(  5, tBasis[   3 ] );
                 mBasis[   4 ]->insert_neighbor(  6, mBasis[  13 ] );
                 mBasis[   4 ]->insert_neighbor(  7, mBasis[  11 ] );

                 // flag this basis
                 mBasis[   4 ]->flag();
             }

         }

         // test if basis 5 exists
         if ( mBasis[   5 ] != nullptr )
         {
             // test if basis 5 has been processed
             if ( ! mBasis[   5 ]->is_flagged() )
             {
                 // link neighbors of basis 5
                 mBasis[   5 ]->insert_neighbor(  0, tBasis[   3 ] );
                 mBasis[   5 ]->insert_neighbor(  1, mBasis[   1 ] );
                 mBasis[   5 ]->insert_neighbor(  2, mBasis[  13 ] );
                 mBasis[   5 ]->insert_neighbor(  3, mBasis[   4 ] );
                 mBasis[   5 ]->insert_neighbor(  4, tBasis[   2 ] );
                 mBasis[   5 ]->insert_neighbor(  5, tBasis[   4 ] );
                 mBasis[   5 ]->insert_neighbor(  6, mBasis[   6 ] );
                 mBasis[   5 ]->insert_neighbor(  7, mBasis[  12 ] );

                 // flag this basis
                 mBasis[   5 ]->flag();
             }

         }

         // test if basis 6 exists
         if ( mBasis[   6 ] != nullptr )
         {
             // test if basis 6 has been processed
             if ( ! mBasis[   6 ]->is_flagged() )
             {
                 // link neighbors of basis 6
                 mBasis[   6 ]->insert_neighbor(  0, mBasis[   1 ] );
                 mBasis[   6 ]->insert_neighbor(  1, tBasis[   9 ] );
                 mBasis[   6 ]->insert_neighbor(  2, mBasis[   7 ] );
                 mBasis[   6 ]->insert_neighbor(  3, mBasis[  13 ] );
                 mBasis[   6 ]->insert_neighbor(  4, mBasis[   5 ] );
                 mBasis[   6 ]->insert_neighbor(  5, tBasis[   7 ] );
                 mBasis[   6 ]->insert_neighbor(  6, tBasis[  11 ] );
                 mBasis[   6 ]->insert_neighbor(  7, mBasis[  14 ] );

                 // flag this basis
                 mBasis[   6 ]->flag();
             }

         }

         // test if basis 7 exists
         if ( mBasis[   7 ] != nullptr )
         {
             // test if basis 7 has been processed
             if ( ! mBasis[   7 ]->is_flagged() )
             {
                 // link neighbors of basis 7
                 mBasis[   7 ]->insert_neighbor(  0, mBasis[   6 ] );
                 mBasis[   7 ]->insert_neighbor(  1, tBasis[  11 ] );
                 mBasis[   7 ]->insert_neighbor(  2, mBasis[   2 ] );
                 mBasis[   7 ]->insert_neighbor(  3, mBasis[  14 ] );
                 mBasis[   7 ]->insert_neighbor(  4, mBasis[  13 ] );
                 mBasis[   7 ]->insert_neighbor(  5, tBasis[   9 ] );
                 mBasis[   7 ]->insert_neighbor(  6, tBasis[  13 ] );
                 mBasis[   7 ]->insert_neighbor(  7, mBasis[   8 ] );

                 // flag this basis
                 mBasis[   7 ]->flag();
             }

         }

         // test if basis 8 exists
         if ( mBasis[   8 ] != nullptr )
         {
             // test if basis 8 has been processed
             if ( ! mBasis[   8 ]->is_flagged() )
             {
                 // link neighbors of basis 8
                 mBasis[   8 ]->insert_neighbor(  0, mBasis[  14 ] );
                 mBasis[   8 ]->insert_neighbor(  1, mBasis[   2 ] );
                 mBasis[   8 ]->insert_neighbor(  2, tBasis[  17 ] );
                 mBasis[   8 ]->insert_neighbor(  3, mBasis[   9 ] );
                 mBasis[   8 ]->insert_neighbor(  4, mBasis[  15 ] );
                 mBasis[   8 ]->insert_neighbor(  5, mBasis[   7 ] );
                 mBasis[   8 ]->insert_neighbor(  6, tBasis[  18 ] );
                 mBasis[   8 ]->insert_neighbor(  7, tBasis[  16 ] );

                 // flag this basis
                 mBasis[   8 ]->flag();
             }

         }

         // test if basis 9 exists
         if ( mBasis[   9 ] != nullptr )
         {
             // test if basis 9 has been processed
             if ( ! mBasis[   9 ]->is_flagged() )
             {
                 // link neighbors of basis 9
                 mBasis[   9 ]->insert_neighbor(  0, mBasis[  15 ] );
                 mBasis[   9 ]->insert_neighbor(  1, mBasis[   8 ] );
                 mBasis[   9 ]->insert_neighbor(  2, tBasis[  16 ] );
                 mBasis[   9 ]->insert_neighbor(  3, mBasis[   3 ] );
                 mBasis[   9 ]->insert_neighbor(  4, mBasis[  10 ] );
                 mBasis[   9 ]->insert_neighbor(  5, mBasis[  14 ] );
                 mBasis[   9 ]->insert_neighbor(  6, tBasis[  17 ] );
                 mBasis[   9 ]->insert_neighbor(  7, tBasis[  15 ] );

                 // flag this basis
                 mBasis[   9 ]->flag();
             }

         }

         // test if basis 10 exists
         if ( mBasis[  10 ] != nullptr )
         {
             // test if basis 10 has been processed
             if ( ! mBasis[  10 ]->is_flagged() )
             {
                 // link neighbors of basis 10
                 mBasis[  10 ]->insert_neighbor(  0, mBasis[  11 ] );
                 mBasis[  10 ]->insert_neighbor(  1, mBasis[  15 ] );
                 mBasis[  10 ]->insert_neighbor(  2, mBasis[   3 ] );
                 mBasis[  10 ]->insert_neighbor(  3, tBasis[  10 ] );
                 mBasis[  10 ]->insert_neighbor(  4, tBasis[   8 ] );
                 mBasis[  10 ]->insert_neighbor(  5, mBasis[  12 ] );
                 mBasis[  10 ]->insert_neighbor(  6, mBasis[   9 ] );
                 mBasis[  10 ]->insert_neighbor(  7, tBasis[  12 ] );

                 // flag this basis
                 mBasis[  10 ]->flag();
             }

         }

         // test if basis 11 exists
         if ( mBasis[  11 ] != nullptr )
         {
             // test if basis 11 has been processed
             if ( ! mBasis[  11 ]->is_flagged() )
             {
                 // link neighbors of basis 11
                 mBasis[  11 ]->insert_neighbor(  0, mBasis[   0 ] );
                 mBasis[  11 ]->insert_neighbor(  1, mBasis[  12 ] );
                 mBasis[  11 ]->insert_neighbor(  2, mBasis[  10 ] );
                 mBasis[  11 ]->insert_neighbor(  3, tBasis[   8 ] );
                 mBasis[  11 ]->insert_neighbor(  4, tBasis[   6 ] );
                 mBasis[  11 ]->insert_neighbor(  5, mBasis[   4 ] );
                 mBasis[  11 ]->insert_neighbor(  6, mBasis[  15 ] );
                 mBasis[  11 ]->insert_neighbor(  7, tBasis[  10 ] );

                 // flag this basis
                 mBasis[  11 ]->flag();
             }

         }

         // test if basis 12 exists
         if ( mBasis[  12 ] != nullptr )
         {
             // test if basis 12 has been processed
             if ( ! mBasis[  12 ]->is_flagged() )
             {
                 // link neighbors of basis 12
                 mBasis[  12 ]->insert_neighbor(  0, mBasis[   4 ] );
                 mBasis[  12 ]->insert_neighbor(  1, mBasis[  13 ] );
                 mBasis[  12 ]->insert_neighbor(  2, mBasis[  15 ] );
                 mBasis[  12 ]->insert_neighbor(  3, mBasis[  11 ] );
                 mBasis[  12 ]->insert_neighbor(  4, mBasis[   0 ] );
                 mBasis[  12 ]->insert_neighbor(  5, mBasis[   5 ] );
                 mBasis[  12 ]->insert_neighbor(  6, mBasis[  14 ] );
                 mBasis[  12 ]->insert_neighbor(  7, mBasis[  10 ] );

                 // flag this basis
                 mBasis[  12 ]->flag();
             }

         }

         // test if basis 13 exists
         if ( mBasis[  13 ] != nullptr )
         {
             // test if basis 13 has been processed
             if ( ! mBasis[  13 ]->is_flagged() )
             {
                 // link neighbors of basis 13
                 mBasis[  13 ]->insert_neighbor(  0, mBasis[   5 ] );
                 mBasis[  13 ]->insert_neighbor(  1, mBasis[   6 ] );
                 mBasis[  13 ]->insert_neighbor(  2, mBasis[  14 ] );
                 mBasis[  13 ]->insert_neighbor(  3, mBasis[  12 ] );
                 mBasis[  13 ]->insert_neighbor(  4, mBasis[   4 ] );
                 mBasis[  13 ]->insert_neighbor(  5, mBasis[   1 ] );
                 mBasis[  13 ]->insert_neighbor(  6, mBasis[   7 ] );
                 mBasis[  13 ]->insert_neighbor(  7, mBasis[  15 ] );

                 // flag this basis
                 mBasis[  13 ]->flag();
             }

         }

         // test if basis 14 exists
         if ( mBasis[  14 ] != nullptr )
         {
             // test if basis 14 has been processed
             if ( ! mBasis[  14 ]->is_flagged() )
             {
                 // link neighbors of basis 14
                 mBasis[  14 ]->insert_neighbor(  0, mBasis[  13 ] );
                 mBasis[  14 ]->insert_neighbor(  1, mBasis[   7 ] );
                 mBasis[  14 ]->insert_neighbor(  2, mBasis[   8 ] );
                 mBasis[  14 ]->insert_neighbor(  3, mBasis[  15 ] );
                 mBasis[  14 ]->insert_neighbor(  4, mBasis[  12 ] );
                 mBasis[  14 ]->insert_neighbor(  5, mBasis[   6 ] );
                 mBasis[  14 ]->insert_neighbor(  6, mBasis[   2 ] );
                 mBasis[  14 ]->insert_neighbor(  7, mBasis[   9 ] );

                 // flag this basis
                 mBasis[  14 ]->flag();
             }

         }

         // test if basis 15 exists
         if ( mBasis[  15 ] != nullptr )
         {
             // test if basis 15 has been processed
             if ( ! mBasis[  15 ]->is_flagged() )
             {
                 // link neighbors of basis 15
                 mBasis[  15 ]->insert_neighbor(  0, mBasis[  12 ] );
                 mBasis[  15 ]->insert_neighbor(  1, mBasis[  14 ] );
                 mBasis[  15 ]->insert_neighbor(  2, mBasis[   9 ] );
                 mBasis[  15 ]->insert_neighbor(  3, mBasis[  10 ] );
                 mBasis[  15 ]->insert_neighbor(  4, mBasis[  11 ] );
                 mBasis[  15 ]->insert_neighbor(  5, mBasis[  13 ] );
                 mBasis[  15 ]->insert_neighbor(  6, mBasis[   8 ] );
                 mBasis[  15 ]->insert_neighbor(  7, mBasis[   3 ] );

                 // flag this basis
                 mBasis[  15 ]->flag();
             }

         }

    }
// ----------------------------------------------------------------------------

    /**
     * Refines a basis of this element
     *
     * @param aBasisNumber Index of the basis to refine
     * @return Number of created bases
     */
    template<>
    luint BSpline_Element< 3, 3, 0 >::refine_basis( uint aBasisNumber )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // get pointer to basis
        Basis* tBasis = mBasis[ aBasisNumber ];

        // test if basis exists
        if ( tBasis != nullptr )
        {
            // test if basis has been refined already
            if ( ! tBasis->has_children() )
            {
                // initialize neighbor container
                Basis* tNeighbors[ 24 ] = { nullptr };

                // populate neighbor container
                get_basis_neighbors_2d( tBasis, 2, tNeighbors );

                // initialize temporary child container
                Basis* tChildren[ 25 ] = { nullptr };

                // test if neighbor 0 exists
                if( tNeighbors[ 0 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 0 ];

                    // test if neighbor 0 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  10 );
                       tChildren[   1 ] = tNeighbor->get_child(  11 );
                       tChildren[   2 ] = tNeighbor->get_child(  12 );
                       tChildren[   3 ] = tNeighbor->get_child(  13 );
                       tChildren[   4 ] = tNeighbor->get_child(  14 );
                       tChildren[   5 ] = tNeighbor->get_child(  15 );
                       tChildren[   6 ] = tNeighbor->get_child(  16 );
                       tChildren[   7 ] = tNeighbor->get_child(  17 );
                       tChildren[   8 ] = tNeighbor->get_child(  18 );
                       tChildren[   9 ] = tNeighbor->get_child(  19 );
                       tChildren[  10 ] = tNeighbor->get_child(  20 );
                       tChildren[  11 ] = tNeighbor->get_child(  21 );
                       tChildren[  12 ] = tNeighbor->get_child(  22 );
                       tChildren[  13 ] = tNeighbor->get_child(  23 );
                       tChildren[  14 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 1 exists
                if( tNeighbors[ 1 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 1 ];

                    // test if neighbor 1 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(   0 );
                       tChildren[   3 ] = tNeighbor->get_child(   1 );
                       tChildren[   4 ] = tNeighbor->get_child(   2 );
                       tChildren[   7 ] = tNeighbor->get_child(   5 );
                       tChildren[   8 ] = tNeighbor->get_child(   6 );
                       tChildren[   9 ] = tNeighbor->get_child(   7 );
                       tChildren[  12 ] = tNeighbor->get_child(  10 );
                       tChildren[  13 ] = tNeighbor->get_child(  11 );
                       tChildren[  14 ] = tNeighbor->get_child(  12 );
                       tChildren[  17 ] = tNeighbor->get_child(  15 );
                       tChildren[  18 ] = tNeighbor->get_child(  16 );
                       tChildren[  19 ] = tNeighbor->get_child(  17 );
                       tChildren[  22 ] = tNeighbor->get_child(  20 );
                       tChildren[  23 ] = tNeighbor->get_child(  21 );
                       tChildren[  24 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 2 exists
                if( tNeighbors[ 2 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 2 ];

                    // test if neighbor 2 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(   0 );
                       tChildren[  11 ] = tNeighbor->get_child(   1 );
                       tChildren[  12 ] = tNeighbor->get_child(   2 );
                       tChildren[  13 ] = tNeighbor->get_child(   3 );
                       tChildren[  14 ] = tNeighbor->get_child(   4 );
                       tChildren[  15 ] = tNeighbor->get_child(   5 );
                       tChildren[  16 ] = tNeighbor->get_child(   6 );
                       tChildren[  17 ] = tNeighbor->get_child(   7 );
                       tChildren[  18 ] = tNeighbor->get_child(   8 );
                       tChildren[  19 ] = tNeighbor->get_child(   9 );
                       tChildren[  20 ] = tNeighbor->get_child(  10 );
                       tChildren[  21 ] = tNeighbor->get_child(  11 );
                       tChildren[  22 ] = tNeighbor->get_child(  12 );
                       tChildren[  23 ] = tNeighbor->get_child(  13 );
                       tChildren[  24 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 3 exists
                if( tNeighbors[ 3 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 3 ];

                    // test if neighbor 3 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(   2 );
                       tChildren[   1 ] = tNeighbor->get_child(   3 );
                       tChildren[   2 ] = tNeighbor->get_child(   4 );
                       tChildren[   5 ] = tNeighbor->get_child(   7 );
                       tChildren[   6 ] = tNeighbor->get_child(   8 );
                       tChildren[   7 ] = tNeighbor->get_child(   9 );
                       tChildren[  10 ] = tNeighbor->get_child(  12 );
                       tChildren[  11 ] = tNeighbor->get_child(  13 );
                       tChildren[  12 ] = tNeighbor->get_child(  14 );
                       tChildren[  15 ] = tNeighbor->get_child(  17 );
                       tChildren[  16 ] = tNeighbor->get_child(  18 );
                       tChildren[  17 ] = tNeighbor->get_child(  19 );
                       tChildren[  20 ] = tNeighbor->get_child(  22 );
                       tChildren[  21 ] = tNeighbor->get_child(  23 );
                       tChildren[  22 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 4 exists
                if( tNeighbors[ 4 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 4 ];

                    // test if neighbor 4 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  12 );
                       tChildren[   1 ] = tNeighbor->get_child(  13 );
                       tChildren[   2 ] = tNeighbor->get_child(  14 );
                       tChildren[   5 ] = tNeighbor->get_child(  17 );
                       tChildren[   6 ] = tNeighbor->get_child(  18 );
                       tChildren[   7 ] = tNeighbor->get_child(  19 );
                       tChildren[  10 ] = tNeighbor->get_child(  22 );
                       tChildren[  11 ] = tNeighbor->get_child(  23 );
                       tChildren[  12 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 5 exists
                if( tNeighbors[ 5 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 5 ];

                    // test if neighbor 5 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  10 );
                       tChildren[   3 ] = tNeighbor->get_child(  11 );
                       tChildren[   4 ] = tNeighbor->get_child(  12 );
                       tChildren[   7 ] = tNeighbor->get_child(  15 );
                       tChildren[   8 ] = tNeighbor->get_child(  16 );
                       tChildren[   9 ] = tNeighbor->get_child(  17 );
                       tChildren[  12 ] = tNeighbor->get_child(  20 );
                       tChildren[  13 ] = tNeighbor->get_child(  21 );
                       tChildren[  14 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 6 exists
                if( tNeighbors[ 6 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 6 ];

                    // test if neighbor 6 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  12 ] = tNeighbor->get_child(   0 );
                       tChildren[  13 ] = tNeighbor->get_child(   1 );
                       tChildren[  14 ] = tNeighbor->get_child(   2 );
                       tChildren[  17 ] = tNeighbor->get_child(   5 );
                       tChildren[  18 ] = tNeighbor->get_child(   6 );
                       tChildren[  19 ] = tNeighbor->get_child(   7 );
                       tChildren[  22 ] = tNeighbor->get_child(  10 );
                       tChildren[  23 ] = tNeighbor->get_child(  11 );
                       tChildren[  24 ] = tNeighbor->get_child(  12 );
                    }
                }

                // test if neighbor 7 exists
                if( tNeighbors[ 7 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 7 ];

                    // test if neighbor 7 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(   2 );
                       tChildren[  11 ] = tNeighbor->get_child(   3 );
                       tChildren[  12 ] = tNeighbor->get_child(   4 );
                       tChildren[  15 ] = tNeighbor->get_child(   7 );
                       tChildren[  16 ] = tNeighbor->get_child(   8 );
                       tChildren[  17 ] = tNeighbor->get_child(   9 );
                       tChildren[  20 ] = tNeighbor->get_child(  12 );
                       tChildren[  21 ] = tNeighbor->get_child(  13 );
                       tChildren[  22 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 8 exists
                if( tNeighbors[ 8 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 8 ];

                    // test if neighbor 8 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 9 exists
                if( tNeighbors[ 9 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 9 ];

                    // test if neighbor 9 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  22 );
                       tChildren[   1 ] = tNeighbor->get_child(  23 );
                       tChildren[   2 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 10 exists
                if( tNeighbors[ 10 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 10 ];

                    // test if neighbor 10 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  20 );
                       tChildren[   1 ] = tNeighbor->get_child(  21 );
                       tChildren[   2 ] = tNeighbor->get_child(  22 );
                       tChildren[   3 ] = tNeighbor->get_child(  23 );
                       tChildren[   4 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 11 exists
                if( tNeighbors[ 11 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 11 ];

                    // test if neighbor 11 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  20 );
                       tChildren[   3 ] = tNeighbor->get_child(  21 );
                       tChildren[   4 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 12 exists
                if( tNeighbors[ 12 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 12 ];

                    // test if neighbor 12 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 13 exists
                if( tNeighbors[ 13 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 13 ];

                    // test if neighbor 13 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  14 );
                       tChildren[   5 ] = tNeighbor->get_child(  19 );
                       tChildren[  10 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 14 exists
                if( tNeighbors[ 14 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 14 ];

                    // test if neighbor 14 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  10 );
                       tChildren[   9 ] = tNeighbor->get_child(  15 );
                       tChildren[  14 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 15 exists
                if( tNeighbors[ 15 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 15 ];

                    // test if neighbor 15 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(   4 );
                       tChildren[   5 ] = tNeighbor->get_child(   9 );
                       tChildren[  10 ] = tNeighbor->get_child(  14 );
                       tChildren[  15 ] = tNeighbor->get_child(  19 );
                       tChildren[  20 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 16 exists
                if( tNeighbors[ 16 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 16 ];

                    // test if neighbor 16 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(   0 );
                       tChildren[   9 ] = tNeighbor->get_child(   5 );
                       tChildren[  14 ] = tNeighbor->get_child(  10 );
                       tChildren[  19 ] = tNeighbor->get_child(  15 );
                       tChildren[  24 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 17 exists
                if( tNeighbors[ 17 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 17 ];

                    // test if neighbor 17 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(   4 );
                       tChildren[  15 ] = tNeighbor->get_child(   9 );
                       tChildren[  20 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 18 exists
                if( tNeighbors[ 18 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 18 ];

                    // test if neighbor 18 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  14 ] = tNeighbor->get_child(   0 );
                       tChildren[  19 ] = tNeighbor->get_child(   5 );
                       tChildren[  24 ] = tNeighbor->get_child(  10 );
                    }
                }

                // test if neighbor 19 exists
                if( tNeighbors[ 19 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 19 ];

                    // test if neighbor 19 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 20 exists
                if( tNeighbors[ 20 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 20 ];

                    // test if neighbor 20 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   2 );
                       tChildren[  21 ] = tNeighbor->get_child(   3 );
                       tChildren[  22 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 21 exists
                if( tNeighbors[ 21 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 21 ];

                    // test if neighbor 21 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   0 );
                       tChildren[  21 ] = tNeighbor->get_child(   1 );
                       tChildren[  22 ] = tNeighbor->get_child(   2 );
                       tChildren[  23 ] = tNeighbor->get_child(   3 );
                       tChildren[  24 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 22 exists
                if( tNeighbors[ 22 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 22 ];

                    // test if neighbor 22 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  22 ] = tNeighbor->get_child(   0 );
                       tChildren[  23 ] = tNeighbor->get_child(   1 );
                       tChildren[  24 ] = tNeighbor->get_child(   2 );
                    }
                }

                // test if neighbor 23 exists
                if( tNeighbors[ 23 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis* tNeighbor = tNeighbors[ 23 ];

                    // test if neighbor 23 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  24 ] = tNeighbor->get_child(   0 );
                    }
                }

               // level of child basis
               uint tLevel = tBasis->get_level() + 1;

               // create container for children
               tBasis->init_children_container();

               // position of basis
               const luint* tParentIJ  = tBasis->get_ijk();

               // minumum i-position
               luint tIMin = 2*tParentIJ[ 0 ];

               // minumum j-position
               luint tJMin = 2*tParentIJ[ 1 ];

               // maximum i-position
               luint tIMax = tIMin + 5;

               // maximum j-position
               luint tJMax = tJMin + 5;

               // initialize counter
               uint tChildIndex = 0;

               // loop over all positions
               for( luint j=tJMin; j<tJMax; ++j )
               {
                   for( luint i=tIMin; i<tIMax; ++i )
                   {
                       // test if child exists
                       if( tChildren[ tChildIndex ] != nullptr )
                       {
                           // insert child
                           tBasis->insert_child( tChildIndex, tChildren[ tChildIndex ] );
                       }
                       else
                       {
                           // calculate i-j position of child
                           luint tIJ[ 2 ] = { i, j };

                           // create child
                           tBasis->insert_child( tChildIndex,
                               new BSpline< 3, 3, 0 >( tIJ, tLevel, gNoProcOwner ) );

                           // increment basis counter
                           tBasisCounter++;
                       }

                       // increment child index
                       tChildIndex++;
                   }
               }
            }
        }
        
        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------

    /**
     * Refines this element.
     *
     * @param aAllElementsOnProc Cell containing all B-spline elements including the aura
     * @return Number of created bases
     */
    template<>
    luint BSpline_Element< 3, 3, 0 >::refine( moris::Cell< Element* > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // refine basis if they have not been refined already
        for( uint k=0; k<16; ++k )
        {
            tBasisCounter += this->refine_basis( k );
        }

        // initialize temporary basis pattern
        Basis* tBasis[ 25 ] = { nullptr };

        // populate basis pattern
        if ( mBasis[   0 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   0 ]->get_child(  18 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   0 ]->get_child(  19 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   0 ]->get_child(  23 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   0 ]->get_child(  24 );
            }
        }

        if ( mBasis[   1 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   1 ]->get_child(  15 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   1 ]->get_child(  16 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   1 ]->get_child(  20 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   1 ]->get_child(  21 );
            }
        }

        if ( mBasis[   2 ] != nullptr )
        {
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   2 ]->get_child(   0 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   2 ]->get_child(   1 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   2 ]->get_child(   5 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   2 ]->get_child(   6 );
            }
        }

        if ( mBasis[   3 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   3 ]->get_child(   3 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   3 ]->get_child(   4 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   3 ]->get_child(   8 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   3 ]->get_child(   9 );
            }
        }

        if ( mBasis[   4 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   4 ]->get_child(  16 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   4 ]->get_child(  17 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   4 ]->get_child(  18 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   4 ]->get_child(  19 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   4 ]->get_child(  21 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   4 ]->get_child(  22 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   4 ]->get_child(  23 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   4 ]->get_child(  24 );
            }
        }

        if ( mBasis[   5 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   5 ]->get_child(  15 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   5 ]->get_child(  16 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   5 ]->get_child(  17 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   5 ]->get_child(  18 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   5 ]->get_child(  20 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   5 ]->get_child(  21 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   5 ]->get_child(  22 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   5 ]->get_child(  23 );
            }
        }

        if ( mBasis[   6 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   6 ]->get_child(   5 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   6 ]->get_child(   6 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   6 ]->get_child(  10 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   6 ]->get_child(  11 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   6 ]->get_child(  15 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   6 ]->get_child(  16 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   6 ]->get_child(  20 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   6 ]->get_child(  21 );
            }
        }

        if ( mBasis[   7 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   7 ]->get_child(   0 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   7 ]->get_child(   1 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   7 ]->get_child(   5 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   7 ]->get_child(   6 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   7 ]->get_child(  10 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   7 ]->get_child(  11 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   7 ]->get_child(  15 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   7 ]->get_child(  16 );
            }
        }

        if ( mBasis[   8 ] != nullptr )
        {
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   8 ]->get_child(   0 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   8 ]->get_child(   1 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   8 ]->get_child(   2 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   8 ]->get_child(   3 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   8 ]->get_child(   5 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   8 ]->get_child(   6 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   8 ]->get_child(   7 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   8 ]->get_child(   8 );
            }
        }

        if ( mBasis[   9 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   9 ]->get_child(   1 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   9 ]->get_child(   2 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   9 ]->get_child(   3 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   9 ]->get_child(   4 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   9 ]->get_child(   6 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   9 ]->get_child(   7 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   9 ]->get_child(   8 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   9 ]->get_child(   9 );
            }
        }

        if ( mBasis[  10 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  10 ]->get_child(   3 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  10 ]->get_child(   4 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  10 ]->get_child(   8 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  10 ]->get_child(   9 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  10 ]->get_child(  13 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  10 ]->get_child(  14 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  10 ]->get_child(  18 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  10 ]->get_child(  19 );
            }
        }

        if ( mBasis[  11 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  11 ]->get_child(   8 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  11 ]->get_child(   9 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  11 ]->get_child(  13 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  11 ]->get_child(  14 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  11 ]->get_child(  18 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  11 ]->get_child(  19 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  11 ]->get_child(  23 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  11 ]->get_child(  24 );
            }
        }

        if ( mBasis[  12 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  12 ]->get_child(   6 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  12 ]->get_child(   7 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  12 ]->get_child(   8 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  12 ]->get_child(   9 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  12 ]->get_child(  11 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  12 ]->get_child(  12 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  12 ]->get_child(  13 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  12 ]->get_child(  14 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  12 ]->get_child(  16 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  12 ]->get_child(  17 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  12 ]->get_child(  18 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  12 ]->get_child(  19 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  12 ]->get_child(  21 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  12 ]->get_child(  22 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  12 ]->get_child(  23 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  12 ]->get_child(  24 );
            }
        }

        if ( mBasis[  13 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  13 ]->get_child(   5 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  13 ]->get_child(   6 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  13 ]->get_child(   7 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  13 ]->get_child(   8 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  13 ]->get_child(  10 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  13 ]->get_child(  11 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  13 ]->get_child(  12 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  13 ]->get_child(  13 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  13 ]->get_child(  15 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  13 ]->get_child(  16 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  13 ]->get_child(  17 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  13 ]->get_child(  18 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  13 ]->get_child(  20 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  13 ]->get_child(  21 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  13 ]->get_child(  22 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  13 ]->get_child(  23 );
            }
        }

        if ( mBasis[  14 ] != nullptr )
        {
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  14 ]->get_child(   0 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  14 ]->get_child(   1 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  14 ]->get_child(   2 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  14 ]->get_child(   3 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  14 ]->get_child(   5 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  14 ]->get_child(   6 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  14 ]->get_child(   7 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  14 ]->get_child(   8 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  14 ]->get_child(  10 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  14 ]->get_child(  11 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  14 ]->get_child(  12 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  14 ]->get_child(  13 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  14 ]->get_child(  15 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  14 ]->get_child(  16 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  14 ]->get_child(  17 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  14 ]->get_child(  18 );
            }
        }

        if ( mBasis[  15 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  15 ]->get_child(   1 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  15 ]->get_child(   2 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  15 ]->get_child(   3 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  15 ]->get_child(   4 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  15 ]->get_child(   6 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  15 ]->get_child(   7 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  15 ]->get_child(   8 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  15 ]->get_child(   9 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  15 ]->get_child(  11 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  15 ]->get_child(  12 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  15 ]->get_child(  13 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  15 ]->get_child(  14 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  15 ]->get_child(  16 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  15 ]->get_child(  17 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  15 ]->get_child(  18 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  15 ]->get_child(  19 );
            }
        }

        // initialize child container
        Element * tChildren[ 4 ];

        // populate child container
        for( uint k=0; k<4; ++k)
        {
            tChildren[ k ] = aAllElementsOnProc(
                mElement->get_child( k )->get_memory_index() );
            tChildren[ k ]->init_basis_container();
        }

        // assign basis to child 1
        tChildren[ 0 ]->insert_basis(  0, tBasis[   0 ] );
        tChildren[ 0 ]->insert_basis(  1, tBasis[   3 ] );
        tChildren[ 0 ]->insert_basis(  2, tBasis[  18 ] );
        tChildren[ 0 ]->insert_basis(  3, tBasis[  15 ] );
        tChildren[ 0 ]->insert_basis(  4, tBasis[   1 ] );
        tChildren[ 0 ]->insert_basis(  5, tBasis[   2 ] );
        tChildren[ 0 ]->insert_basis(  6, tBasis[   8 ] );
        tChildren[ 0 ]->insert_basis(  7, tBasis[  13 ] );
        tChildren[ 0 ]->insert_basis(  8, tBasis[  17 ] );
        tChildren[ 0 ]->insert_basis(  9, tBasis[  16 ] );
        tChildren[ 0 ]->insert_basis( 10, tBasis[  10 ] );
        tChildren[ 0 ]->insert_basis( 11, tBasis[   5 ] );
        tChildren[ 0 ]->insert_basis( 12, tBasis[   6 ] );
        tChildren[ 0 ]->insert_basis( 13, tBasis[   7 ] );
        tChildren[ 0 ]->insert_basis( 14, tBasis[  12 ] );
        tChildren[ 0 ]->insert_basis( 15, tBasis[  11 ] );

        // assign basis to child 2
        tChildren[ 1 ]->insert_basis(  0, tBasis[   1 ] );
        tChildren[ 1 ]->insert_basis(  1, tBasis[   4 ] );
        tChildren[ 1 ]->insert_basis(  2, tBasis[  19 ] );
        tChildren[ 1 ]->insert_basis(  3, tBasis[  16 ] );
        tChildren[ 1 ]->insert_basis(  4, tBasis[   2 ] );
        tChildren[ 1 ]->insert_basis(  5, tBasis[   3 ] );
        tChildren[ 1 ]->insert_basis(  6, tBasis[   9 ] );
        tChildren[ 1 ]->insert_basis(  7, tBasis[  14 ] );
        tChildren[ 1 ]->insert_basis(  8, tBasis[  18 ] );
        tChildren[ 1 ]->insert_basis(  9, tBasis[  17 ] );
        tChildren[ 1 ]->insert_basis( 10, tBasis[  11 ] );
        tChildren[ 1 ]->insert_basis( 11, tBasis[   6 ] );
        tChildren[ 1 ]->insert_basis( 12, tBasis[   7 ] );
        tChildren[ 1 ]->insert_basis( 13, tBasis[   8 ] );
        tChildren[ 1 ]->insert_basis( 14, tBasis[  13 ] );
        tChildren[ 1 ]->insert_basis( 15, tBasis[  12 ] );

        // assign basis to child 3
        tChildren[ 2 ]->insert_basis(  0, tBasis[   5 ] );
        tChildren[ 2 ]->insert_basis(  1, tBasis[   8 ] );
        tChildren[ 2 ]->insert_basis(  2, tBasis[  23 ] );
        tChildren[ 2 ]->insert_basis(  3, tBasis[  20 ] );
        tChildren[ 2 ]->insert_basis(  4, tBasis[   6 ] );
        tChildren[ 2 ]->insert_basis(  5, tBasis[   7 ] );
        tChildren[ 2 ]->insert_basis(  6, tBasis[  13 ] );
        tChildren[ 2 ]->insert_basis(  7, tBasis[  18 ] );
        tChildren[ 2 ]->insert_basis(  8, tBasis[  22 ] );
        tChildren[ 2 ]->insert_basis(  9, tBasis[  21 ] );
        tChildren[ 2 ]->insert_basis( 10, tBasis[  15 ] );
        tChildren[ 2 ]->insert_basis( 11, tBasis[  10 ] );
        tChildren[ 2 ]->insert_basis( 12, tBasis[  11 ] );
        tChildren[ 2 ]->insert_basis( 13, tBasis[  12 ] );
        tChildren[ 2 ]->insert_basis( 14, tBasis[  17 ] );
        tChildren[ 2 ]->insert_basis( 15, tBasis[  16 ] );

        // assign basis to child 4
        tChildren[ 3 ]->insert_basis(  0, tBasis[   6 ] );
        tChildren[ 3 ]->insert_basis(  1, tBasis[   9 ] );
        tChildren[ 3 ]->insert_basis(  2, tBasis[  24 ] );
        tChildren[ 3 ]->insert_basis(  3, tBasis[  21 ] );
        tChildren[ 3 ]->insert_basis(  4, tBasis[   7 ] );
        tChildren[ 3 ]->insert_basis(  5, tBasis[   8 ] );
        tChildren[ 3 ]->insert_basis(  6, tBasis[  14 ] );
        tChildren[ 3 ]->insert_basis(  7, tBasis[  19 ] );
        tChildren[ 3 ]->insert_basis(  8, tBasis[  23 ] );
        tChildren[ 3 ]->insert_basis(  9, tBasis[  22 ] );
        tChildren[ 3 ]->insert_basis( 10, tBasis[  16 ] );
        tChildren[ 3 ]->insert_basis( 11, tBasis[  11 ] );
        tChildren[ 3 ]->insert_basis( 12, tBasis[  12 ] );
        tChildren[ 3 ]->insert_basis( 13, tBasis[  13 ] );
        tChildren[ 3 ]->insert_basis( 14, tBasis[  18 ] );
        tChildren[ 3 ]->insert_basis( 15, tBasis[  17 ] );

        // set basis flag of element
        mChildrenBasisFlag = true;
        
        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD16_HPP_ */

