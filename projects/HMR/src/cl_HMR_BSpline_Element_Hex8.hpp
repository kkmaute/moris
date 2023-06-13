/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element_Hex8.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX8_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX8_HPP_

#include "cl_HMR_BSpline_Element.hpp"
#include "fn_HMR_get_basis_neighbors_3d.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------

    /**
    * Returns the geometry type of this element
    *
    * @return mtk::Geometry_Type
    */
    template<>
    mtk::Geometry_Type
    BSpline_Element< 3, 8 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::HEX;
    }

// ----------------------------------------------------------------------------

    /**
     * node IDs needed for VTK output
     *
     * @param[out] Matrix< DDLUMat >
     *
     * @return void
     *
     */
    template<>
    void
    BSpline_Element< 3, 8 >::get_basis_indices_for_vtk(
        Matrix< DDLUMat > & aBasis )
    {
        // loop over all nodes
        for( uint k=0; k<8; ++k )
        {
            aBasis( k ) = mBasis[ k ]->get_memory_index();
        }
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
    BSpline_Element< 3, 8 >::get_ijk_of_basis(
        uint aBasisNumber,
        luint      * aIJK )
    {
        // get element local coordinate
        switch ( aBasisNumber )
        {
            case( 0 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 1 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 2 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 3 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 4 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 5 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 6 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 7 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
        }

        // get position of element on background mesh
        const luint * tElIJK = mElement->get_ijk();

        // add element offset
        aIJK[ 0 ] += tElIJK[ 0 ];
        aIJK[ 1 ] += tElIJK[ 1 ];
        aIJK[ 2 ] += tElIJK[ 2 ];
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
    BSpline_Element< 3, 8 >::link_basis_with_neighbors(
          moris::Cell< Element* > & aAllElementsOnProc )
    {
         // initialize frame of basis around basis from this element
         Basis* tBasis[ 56 ] = { nullptr };

         // get pointer to neighbor  0
         Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

         // test if neighbor  0 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  17 ] = tNeighbor->get_basis(   0 );
             tBasis[  18 ] = tNeighbor->get_basis(   1 );
             tBasis[  29 ] = tNeighbor->get_basis(   4 );
             tBasis[  30 ] = tNeighbor->get_basis(   5 );
         }

         // get pointer to neighbor  1
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

         // test if neighbor  1 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  21 ] = tNeighbor->get_basis(   1 );
             tBasis[  23 ] = tNeighbor->get_basis(   2 );
             tBasis[  33 ] = tNeighbor->get_basis(   5 );
             tBasis[  35 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  2
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

         // test if neighbor  2 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  25 ] = tNeighbor->get_basis(   3 );
             tBasis[  26 ] = tNeighbor->get_basis(   2 );
             tBasis[  37 ] = tNeighbor->get_basis(   7 );
             tBasis[  38 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  3
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

         // test if neighbor  3 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  20 ] = tNeighbor->get_basis(   0 );
             tBasis[  22 ] = tNeighbor->get_basis(   3 );
             tBasis[  32 ] = tNeighbor->get_basis(   4 );
             tBasis[  34 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  4
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

         // test if neighbor  4 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   5 ] = tNeighbor->get_basis(   0 );
             tBasis[   6 ] = tNeighbor->get_basis(   1 );
             tBasis[   9 ] = tNeighbor->get_basis(   3 );
             tBasis[  10 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  5
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

         // test if neighbor  5 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  45 ] = tNeighbor->get_basis(   4 );
             tBasis[  46 ] = tNeighbor->get_basis(   5 );
             tBasis[  49 ] = tNeighbor->get_basis(   7 );
             tBasis[  50 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  6
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

         // test if neighbor  6 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   1 ] = tNeighbor->get_basis(   0 );
             tBasis[   2 ] = tNeighbor->get_basis(   1 );
             tBasis[   5 ] = tNeighbor->get_basis(   3 );
             tBasis[   6 ] = tNeighbor->get_basis(   2 );
             tBasis[  17 ] = tNeighbor->get_basis(   4 );
             tBasis[  18 ] = tNeighbor->get_basis(   5 );
         }

         // get pointer to neighbor  7
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

         // test if neighbor  7 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   6 ] = tNeighbor->get_basis(   0 );
             tBasis[   7 ] = tNeighbor->get_basis(   1 );
             tBasis[  10 ] = tNeighbor->get_basis(   3 );
             tBasis[  11 ] = tNeighbor->get_basis(   2 );
             tBasis[  21 ] = tNeighbor->get_basis(   5 );
             tBasis[  23 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  8
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 8 );

         // test if neighbor  8 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   9 ] = tNeighbor->get_basis(   0 );
             tBasis[  10 ] = tNeighbor->get_basis(   1 );
             tBasis[  13 ] = tNeighbor->get_basis(   3 );
             tBasis[  14 ] = tNeighbor->get_basis(   2 );
             tBasis[  25 ] = tNeighbor->get_basis(   7 );
             tBasis[  26 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  9
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 9 );

         // test if neighbor  9 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   4 ] = tNeighbor->get_basis(   0 );
             tBasis[   5 ] = tNeighbor->get_basis(   1 );
             tBasis[   8 ] = tNeighbor->get_basis(   3 );
             tBasis[   9 ] = tNeighbor->get_basis(   2 );
             tBasis[  20 ] = tNeighbor->get_basis(   4 );
             tBasis[  22 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  10
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 10 );

         // test if neighbor  10 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  16 ] = tNeighbor->get_basis(   0 );
             tBasis[  17 ] = tNeighbor->get_basis(   1 );
             tBasis[  20 ] = tNeighbor->get_basis(   3 );
             tBasis[  28 ] = tNeighbor->get_basis(   4 );
             tBasis[  29 ] = tNeighbor->get_basis(   5 );
             tBasis[  32 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  11
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 11 );

         // test if neighbor  11 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  18 ] = tNeighbor->get_basis(   0 );
             tBasis[  19 ] = tNeighbor->get_basis(   1 );
             tBasis[  21 ] = tNeighbor->get_basis(   2 );
             tBasis[  30 ] = tNeighbor->get_basis(   4 );
             tBasis[  31 ] = tNeighbor->get_basis(   5 );
             tBasis[  33 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  12
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 12 );

         // test if neighbor  12 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  23 ] = tNeighbor->get_basis(   1 );
             tBasis[  26 ] = tNeighbor->get_basis(   3 );
             tBasis[  27 ] = tNeighbor->get_basis(   2 );
             tBasis[  35 ] = tNeighbor->get_basis(   5 );
             tBasis[  38 ] = tNeighbor->get_basis(   7 );
             tBasis[  39 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  13
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 13 );

         // test if neighbor  13 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  22 ] = tNeighbor->get_basis(   0 );
             tBasis[  24 ] = tNeighbor->get_basis(   3 );
             tBasis[  25 ] = tNeighbor->get_basis(   2 );
             tBasis[  34 ] = tNeighbor->get_basis(   4 );
             tBasis[  36 ] = tNeighbor->get_basis(   7 );
             tBasis[  37 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  14
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 14 );

         // test if neighbor  14 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  29 ] = tNeighbor->get_basis(   0 );
             tBasis[  30 ] = tNeighbor->get_basis(   1 );
             tBasis[  41 ] = tNeighbor->get_basis(   4 );
             tBasis[  42 ] = tNeighbor->get_basis(   5 );
             tBasis[  45 ] = tNeighbor->get_basis(   7 );
             tBasis[  46 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  15
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 15 );

         // test if neighbor  15 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  33 ] = tNeighbor->get_basis(   1 );
             tBasis[  35 ] = tNeighbor->get_basis(   2 );
             tBasis[  46 ] = tNeighbor->get_basis(   4 );
             tBasis[  47 ] = tNeighbor->get_basis(   5 );
             tBasis[  50 ] = tNeighbor->get_basis(   7 );
             tBasis[  51 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  16
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 16 );

         // test if neighbor  16 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  37 ] = tNeighbor->get_basis(   3 );
             tBasis[  38 ] = tNeighbor->get_basis(   2 );
             tBasis[  49 ] = tNeighbor->get_basis(   4 );
             tBasis[  50 ] = tNeighbor->get_basis(   5 );
             tBasis[  53 ] = tNeighbor->get_basis(   7 );
             tBasis[  54 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  17
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 17 );

         // test if neighbor  17 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  32 ] = tNeighbor->get_basis(   0 );
             tBasis[  34 ] = tNeighbor->get_basis(   3 );
             tBasis[  44 ] = tNeighbor->get_basis(   4 );
             tBasis[  45 ] = tNeighbor->get_basis(   5 );
             tBasis[  48 ] = tNeighbor->get_basis(   7 );
             tBasis[  49 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  18
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 18 );

         // test if neighbor  18 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   0 ] = tNeighbor->get_basis(   0 );
             tBasis[   1 ] = tNeighbor->get_basis(   1 );
             tBasis[   4 ] = tNeighbor->get_basis(   3 );
             tBasis[   5 ] = tNeighbor->get_basis(   2 );
             tBasis[  16 ] = tNeighbor->get_basis(   4 );
             tBasis[  17 ] = tNeighbor->get_basis(   5 );
             tBasis[  20 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  19
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 19 );

         // test if neighbor  19 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   2 ] = tNeighbor->get_basis(   0 );
             tBasis[   3 ] = tNeighbor->get_basis(   1 );
             tBasis[   6 ] = tNeighbor->get_basis(   3 );
             tBasis[   7 ] = tNeighbor->get_basis(   2 );
             tBasis[  18 ] = tNeighbor->get_basis(   4 );
             tBasis[  19 ] = tNeighbor->get_basis(   5 );
             tBasis[  21 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  20
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 20 );

         // test if neighbor  20 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  10 ] = tNeighbor->get_basis(   0 );
             tBasis[  11 ] = tNeighbor->get_basis(   1 );
             tBasis[  14 ] = tNeighbor->get_basis(   3 );
             tBasis[  15 ] = tNeighbor->get_basis(   2 );
             tBasis[  23 ] = tNeighbor->get_basis(   5 );
             tBasis[  26 ] = tNeighbor->get_basis(   7 );
             tBasis[  27 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  21
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 21 );

         // test if neighbor  21 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   8 ] = tNeighbor->get_basis(   0 );
             tBasis[   9 ] = tNeighbor->get_basis(   1 );
             tBasis[  12 ] = tNeighbor->get_basis(   3 );
             tBasis[  13 ] = tNeighbor->get_basis(   2 );
             tBasis[  22 ] = tNeighbor->get_basis(   4 );
             tBasis[  24 ] = tNeighbor->get_basis(   7 );
             tBasis[  25 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  22
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 22 );

         // test if neighbor  22 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  28 ] = tNeighbor->get_basis(   0 );
             tBasis[  29 ] = tNeighbor->get_basis(   1 );
             tBasis[  32 ] = tNeighbor->get_basis(   3 );
             tBasis[  40 ] = tNeighbor->get_basis(   4 );
             tBasis[  41 ] = tNeighbor->get_basis(   5 );
             tBasis[  44 ] = tNeighbor->get_basis(   7 );
             tBasis[  45 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  23
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 23 );

         // test if neighbor  23 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  30 ] = tNeighbor->get_basis(   0 );
             tBasis[  31 ] = tNeighbor->get_basis(   1 );
             tBasis[  33 ] = tNeighbor->get_basis(   2 );
             tBasis[  42 ] = tNeighbor->get_basis(   4 );
             tBasis[  43 ] = tNeighbor->get_basis(   5 );
             tBasis[  46 ] = tNeighbor->get_basis(   7 );
             tBasis[  47 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  24
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 24 );

         // test if neighbor  24 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  35 ] = tNeighbor->get_basis(   1 );
             tBasis[  38 ] = tNeighbor->get_basis(   3 );
             tBasis[  39 ] = tNeighbor->get_basis(   2 );
             tBasis[  50 ] = tNeighbor->get_basis(   4 );
             tBasis[  51 ] = tNeighbor->get_basis(   5 );
             tBasis[  54 ] = tNeighbor->get_basis(   7 );
             tBasis[  55 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  25
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 25 );

         // test if neighbor  25 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  34 ] = tNeighbor->get_basis(   0 );
             tBasis[  36 ] = tNeighbor->get_basis(   3 );
             tBasis[  37 ] = tNeighbor->get_basis(   2 );
             tBasis[  48 ] = tNeighbor->get_basis(   4 );
             tBasis[  49 ] = tNeighbor->get_basis(   5 );
             tBasis[  52 ] = tNeighbor->get_basis(   7 );
             tBasis[  53 ] = tNeighbor->get_basis(   6 );
         }

         // test if basis 0 exists
         if ( mBasis[   0 ] != nullptr )
         {
             // test if basis 0 has been processed
             if ( ! mBasis[   0 ]->is_flagged() )
             {
                 // link neighbors of basis 0
                 mBasis[   0 ]->insert_neighbor(  0, tBasis[  17 ] );
                 mBasis[   0 ]->insert_neighbor(  1, mBasis[   1 ] );
                 mBasis[   0 ]->insert_neighbor(  2, mBasis[   3 ] );
                 mBasis[   0 ]->insert_neighbor(  3, tBasis[  20 ] );
                 mBasis[   0 ]->insert_neighbor(  4, tBasis[   5 ] );
                 mBasis[   0 ]->insert_neighbor(  5, mBasis[   4 ] );
                 mBasis[   0 ]->insert_neighbor(  6, tBasis[   1 ] );
                 mBasis[   0 ]->insert_neighbor(  7, tBasis[   6 ] );
                 mBasis[   0 ]->insert_neighbor(  8, tBasis[   9 ] );
                 mBasis[   0 ]->insert_neighbor(  9, tBasis[   4 ] );
                 mBasis[   0 ]->insert_neighbor( 10, tBasis[  16 ] );
                 mBasis[   0 ]->insert_neighbor( 11, tBasis[  18 ] );
                 mBasis[   0 ]->insert_neighbor( 12, mBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor( 13, tBasis[  22 ] );
                 mBasis[   0 ]->insert_neighbor( 14, tBasis[  29 ] );
                 mBasis[   0 ]->insert_neighbor( 15, mBasis[   5 ] );
                 mBasis[   0 ]->insert_neighbor( 16, mBasis[   7 ] );
                 mBasis[   0 ]->insert_neighbor( 17, tBasis[  32 ] );
                 mBasis[   0 ]->insert_neighbor( 18, tBasis[   0 ] );
                 mBasis[   0 ]->insert_neighbor( 19, tBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor( 20, tBasis[  10 ] );
                 mBasis[   0 ]->insert_neighbor( 21, tBasis[   8 ] );
                 mBasis[   0 ]->insert_neighbor( 22, tBasis[  28 ] );
                 mBasis[   0 ]->insert_neighbor( 23, tBasis[  30 ] );
                 mBasis[   0 ]->insert_neighbor( 24, mBasis[   6 ] );
                 mBasis[   0 ]->insert_neighbor( 25, tBasis[  34 ] );

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
                 mBasis[   1 ]->insert_neighbor(  0, tBasis[  18 ] );
                 mBasis[   1 ]->insert_neighbor(  1, tBasis[  21 ] );
                 mBasis[   1 ]->insert_neighbor(  2, mBasis[   2 ] );
                 mBasis[   1 ]->insert_neighbor(  3, mBasis[   0 ] );
                 mBasis[   1 ]->insert_neighbor(  4, tBasis[   6 ] );
                 mBasis[   1 ]->insert_neighbor(  5, mBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor(  6, tBasis[   2 ] );
                 mBasis[   1 ]->insert_neighbor(  7, tBasis[   7 ] );
                 mBasis[   1 ]->insert_neighbor(  8, tBasis[  10 ] );
                 mBasis[   1 ]->insert_neighbor(  9, tBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor( 10, tBasis[  17 ] );
                 mBasis[   1 ]->insert_neighbor( 11, tBasis[  19 ] );
                 mBasis[   1 ]->insert_neighbor( 12, tBasis[  23 ] );
                 mBasis[   1 ]->insert_neighbor( 13, mBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor( 14, tBasis[  30 ] );
                 mBasis[   1 ]->insert_neighbor( 15, tBasis[  33 ] );
                 mBasis[   1 ]->insert_neighbor( 16, mBasis[   6 ] );
                 mBasis[   1 ]->insert_neighbor( 17, mBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor( 18, tBasis[   1 ] );
                 mBasis[   1 ]->insert_neighbor( 19, tBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor( 20, tBasis[  11 ] );
                 mBasis[   1 ]->insert_neighbor( 21, tBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor( 22, tBasis[  29 ] );
                 mBasis[   1 ]->insert_neighbor( 23, tBasis[  31 ] );
                 mBasis[   1 ]->insert_neighbor( 24, tBasis[  35 ] );
                 mBasis[   1 ]->insert_neighbor( 25, mBasis[   7 ] );

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
                 mBasis[   2 ]->insert_neighbor(  0, mBasis[   1 ] );
                 mBasis[   2 ]->insert_neighbor(  1, tBasis[  23 ] );
                 mBasis[   2 ]->insert_neighbor(  2, tBasis[  26 ] );
                 mBasis[   2 ]->insert_neighbor(  3, mBasis[   3 ] );
                 mBasis[   2 ]->insert_neighbor(  4, tBasis[  10 ] );
                 mBasis[   2 ]->insert_neighbor(  5, mBasis[   6 ] );
                 mBasis[   2 ]->insert_neighbor(  6, tBasis[   6 ] );
                 mBasis[   2 ]->insert_neighbor(  7, tBasis[  11 ] );
                 mBasis[   2 ]->insert_neighbor(  8, tBasis[  14 ] );
                 mBasis[   2 ]->insert_neighbor(  9, tBasis[   9 ] );
                 mBasis[   2 ]->insert_neighbor( 10, mBasis[   0 ] );
                 mBasis[   2 ]->insert_neighbor( 11, tBasis[  21 ] );
                 mBasis[   2 ]->insert_neighbor( 12, tBasis[  27 ] );
                 mBasis[   2 ]->insert_neighbor( 13, tBasis[  25 ] );
                 mBasis[   2 ]->insert_neighbor( 14, mBasis[   5 ] );
                 mBasis[   2 ]->insert_neighbor( 15, tBasis[  35 ] );
                 mBasis[   2 ]->insert_neighbor( 16, tBasis[  38 ] );
                 mBasis[   2 ]->insert_neighbor( 17, mBasis[   7 ] );
                 mBasis[   2 ]->insert_neighbor( 18, tBasis[   5 ] );
                 mBasis[   2 ]->insert_neighbor( 19, tBasis[   7 ] );
                 mBasis[   2 ]->insert_neighbor( 20, tBasis[  15 ] );
                 mBasis[   2 ]->insert_neighbor( 21, tBasis[  13 ] );
                 mBasis[   2 ]->insert_neighbor( 22, mBasis[   4 ] );
                 mBasis[   2 ]->insert_neighbor( 23, tBasis[  33 ] );
                 mBasis[   2 ]->insert_neighbor( 24, tBasis[  39 ] );
                 mBasis[   2 ]->insert_neighbor( 25, tBasis[  37 ] );

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
                 mBasis[   3 ]->insert_neighbor(  0, mBasis[   0 ] );
                 mBasis[   3 ]->insert_neighbor(  1, mBasis[   2 ] );
                 mBasis[   3 ]->insert_neighbor(  2, tBasis[  25 ] );
                 mBasis[   3 ]->insert_neighbor(  3, tBasis[  22 ] );
                 mBasis[   3 ]->insert_neighbor(  4, tBasis[   9 ] );
                 mBasis[   3 ]->insert_neighbor(  5, mBasis[   7 ] );
                 mBasis[   3 ]->insert_neighbor(  6, tBasis[   5 ] );
                 mBasis[   3 ]->insert_neighbor(  7, tBasis[  10 ] );
                 mBasis[   3 ]->insert_neighbor(  8, tBasis[  13 ] );
                 mBasis[   3 ]->insert_neighbor(  9, tBasis[   8 ] );
                 mBasis[   3 ]->insert_neighbor( 10, tBasis[  20 ] );
                 mBasis[   3 ]->insert_neighbor( 11, mBasis[   1 ] );
                 mBasis[   3 ]->insert_neighbor( 12, tBasis[  26 ] );
                 mBasis[   3 ]->insert_neighbor( 13, tBasis[  24 ] );
                 mBasis[   3 ]->insert_neighbor( 14, mBasis[   4 ] );
                 mBasis[   3 ]->insert_neighbor( 15, mBasis[   6 ] );
                 mBasis[   3 ]->insert_neighbor( 16, tBasis[  37 ] );
                 mBasis[   3 ]->insert_neighbor( 17, tBasis[  34 ] );
                 mBasis[   3 ]->insert_neighbor( 18, tBasis[   4 ] );
                 mBasis[   3 ]->insert_neighbor( 19, tBasis[   6 ] );
                 mBasis[   3 ]->insert_neighbor( 20, tBasis[  14 ] );
                 mBasis[   3 ]->insert_neighbor( 21, tBasis[  12 ] );
                 mBasis[   3 ]->insert_neighbor( 22, tBasis[  32 ] );
                 mBasis[   3 ]->insert_neighbor( 23, mBasis[   5 ] );
                 mBasis[   3 ]->insert_neighbor( 24, tBasis[  38 ] );
                 mBasis[   3 ]->insert_neighbor( 25, tBasis[  36 ] );

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
                 mBasis[   4 ]->insert_neighbor(  0, tBasis[  29 ] );
                 mBasis[   4 ]->insert_neighbor(  1, mBasis[   5 ] );
                 mBasis[   4 ]->insert_neighbor(  2, mBasis[   7 ] );
                 mBasis[   4 ]->insert_neighbor(  3, tBasis[  32 ] );
                 mBasis[   4 ]->insert_neighbor(  4, mBasis[   0 ] );
                 mBasis[   4 ]->insert_neighbor(  5, tBasis[  45 ] );
                 mBasis[   4 ]->insert_neighbor(  6, tBasis[  17 ] );
                 mBasis[   4 ]->insert_neighbor(  7, mBasis[   1 ] );
                 mBasis[   4 ]->insert_neighbor(  8, mBasis[   3 ] );
                 mBasis[   4 ]->insert_neighbor(  9, tBasis[  20 ] );
                 mBasis[   4 ]->insert_neighbor( 10, tBasis[  28 ] );
                 mBasis[   4 ]->insert_neighbor( 11, tBasis[  30 ] );
                 mBasis[   4 ]->insert_neighbor( 12, mBasis[   6 ] );
                 mBasis[   4 ]->insert_neighbor( 13, tBasis[  34 ] );
                 mBasis[   4 ]->insert_neighbor( 14, tBasis[  41 ] );
                 mBasis[   4 ]->insert_neighbor( 15, tBasis[  46 ] );
                 mBasis[   4 ]->insert_neighbor( 16, tBasis[  49 ] );
                 mBasis[   4 ]->insert_neighbor( 17, tBasis[  44 ] );
                 mBasis[   4 ]->insert_neighbor( 18, tBasis[  16 ] );
                 mBasis[   4 ]->insert_neighbor( 19, tBasis[  18 ] );
                 mBasis[   4 ]->insert_neighbor( 20, mBasis[   2 ] );
                 mBasis[   4 ]->insert_neighbor( 21, tBasis[  22 ] );
                 mBasis[   4 ]->insert_neighbor( 22, tBasis[  40 ] );
                 mBasis[   4 ]->insert_neighbor( 23, tBasis[  42 ] );
                 mBasis[   4 ]->insert_neighbor( 24, tBasis[  50 ] );
                 mBasis[   4 ]->insert_neighbor( 25, tBasis[  48 ] );

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
                 mBasis[   5 ]->insert_neighbor(  0, tBasis[  30 ] );
                 mBasis[   5 ]->insert_neighbor(  1, tBasis[  33 ] );
                 mBasis[   5 ]->insert_neighbor(  2, mBasis[   6 ] );
                 mBasis[   5 ]->insert_neighbor(  3, mBasis[   4 ] );
                 mBasis[   5 ]->insert_neighbor(  4, mBasis[   1 ] );
                 mBasis[   5 ]->insert_neighbor(  5, tBasis[  46 ] );
                 mBasis[   5 ]->insert_neighbor(  6, tBasis[  18 ] );
                 mBasis[   5 ]->insert_neighbor(  7, tBasis[  21 ] );
                 mBasis[   5 ]->insert_neighbor(  8, mBasis[   2 ] );
                 mBasis[   5 ]->insert_neighbor(  9, mBasis[   0 ] );
                 mBasis[   5 ]->insert_neighbor( 10, tBasis[  29 ] );
                 mBasis[   5 ]->insert_neighbor( 11, tBasis[  31 ] );
                 mBasis[   5 ]->insert_neighbor( 12, tBasis[  35 ] );
                 mBasis[   5 ]->insert_neighbor( 13, mBasis[   7 ] );
                 mBasis[   5 ]->insert_neighbor( 14, tBasis[  42 ] );
                 mBasis[   5 ]->insert_neighbor( 15, tBasis[  47 ] );
                 mBasis[   5 ]->insert_neighbor( 16, tBasis[  50 ] );
                 mBasis[   5 ]->insert_neighbor( 17, tBasis[  45 ] );
                 mBasis[   5 ]->insert_neighbor( 18, tBasis[  17 ] );
                 mBasis[   5 ]->insert_neighbor( 19, tBasis[  19 ] );
                 mBasis[   5 ]->insert_neighbor( 20, tBasis[  23 ] );
                 mBasis[   5 ]->insert_neighbor( 21, mBasis[   3 ] );
                 mBasis[   5 ]->insert_neighbor( 22, tBasis[  41 ] );
                 mBasis[   5 ]->insert_neighbor( 23, tBasis[  43 ] );
                 mBasis[   5 ]->insert_neighbor( 24, tBasis[  51 ] );
                 mBasis[   5 ]->insert_neighbor( 25, tBasis[  49 ] );

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
                 mBasis[   6 ]->insert_neighbor(  0, mBasis[   5 ] );
                 mBasis[   6 ]->insert_neighbor(  1, tBasis[  35 ] );
                 mBasis[   6 ]->insert_neighbor(  2, tBasis[  38 ] );
                 mBasis[   6 ]->insert_neighbor(  3, mBasis[   7 ] );
                 mBasis[   6 ]->insert_neighbor(  4, mBasis[   2 ] );
                 mBasis[   6 ]->insert_neighbor(  5, tBasis[  50 ] );
                 mBasis[   6 ]->insert_neighbor(  6, mBasis[   1 ] );
                 mBasis[   6 ]->insert_neighbor(  7, tBasis[  23 ] );
                 mBasis[   6 ]->insert_neighbor(  8, tBasis[  26 ] );
                 mBasis[   6 ]->insert_neighbor(  9, mBasis[   3 ] );
                 mBasis[   6 ]->insert_neighbor( 10, mBasis[   4 ] );
                 mBasis[   6 ]->insert_neighbor( 11, tBasis[  33 ] );
                 mBasis[   6 ]->insert_neighbor( 12, tBasis[  39 ] );
                 mBasis[   6 ]->insert_neighbor( 13, tBasis[  37 ] );
                 mBasis[   6 ]->insert_neighbor( 14, tBasis[  46 ] );
                 mBasis[   6 ]->insert_neighbor( 15, tBasis[  51 ] );
                 mBasis[   6 ]->insert_neighbor( 16, tBasis[  54 ] );
                 mBasis[   6 ]->insert_neighbor( 17, tBasis[  49 ] );
                 mBasis[   6 ]->insert_neighbor( 18, mBasis[   0 ] );
                 mBasis[   6 ]->insert_neighbor( 19, tBasis[  21 ] );
                 mBasis[   6 ]->insert_neighbor( 20, tBasis[  27 ] );
                 mBasis[   6 ]->insert_neighbor( 21, tBasis[  25 ] );
                 mBasis[   6 ]->insert_neighbor( 22, tBasis[  45 ] );
                 mBasis[   6 ]->insert_neighbor( 23, tBasis[  47 ] );
                 mBasis[   6 ]->insert_neighbor( 24, tBasis[  55 ] );
                 mBasis[   6 ]->insert_neighbor( 25, tBasis[  53 ] );

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
                 mBasis[   7 ]->insert_neighbor(  0, mBasis[   4 ] );
                 mBasis[   7 ]->insert_neighbor(  1, mBasis[   6 ] );
                 mBasis[   7 ]->insert_neighbor(  2, tBasis[  37 ] );
                 mBasis[   7 ]->insert_neighbor(  3, tBasis[  34 ] );
                 mBasis[   7 ]->insert_neighbor(  4, mBasis[   3 ] );
                 mBasis[   7 ]->insert_neighbor(  5, tBasis[  49 ] );
                 mBasis[   7 ]->insert_neighbor(  6, mBasis[   0 ] );
                 mBasis[   7 ]->insert_neighbor(  7, mBasis[   2 ] );
                 mBasis[   7 ]->insert_neighbor(  8, tBasis[  25 ] );
                 mBasis[   7 ]->insert_neighbor(  9, tBasis[  22 ] );
                 mBasis[   7 ]->insert_neighbor( 10, tBasis[  32 ] );
                 mBasis[   7 ]->insert_neighbor( 11, mBasis[   5 ] );
                 mBasis[   7 ]->insert_neighbor( 12, tBasis[  38 ] );
                 mBasis[   7 ]->insert_neighbor( 13, tBasis[  36 ] );
                 mBasis[   7 ]->insert_neighbor( 14, tBasis[  45 ] );
                 mBasis[   7 ]->insert_neighbor( 15, tBasis[  50 ] );
                 mBasis[   7 ]->insert_neighbor( 16, tBasis[  53 ] );
                 mBasis[   7 ]->insert_neighbor( 17, tBasis[  48 ] );
                 mBasis[   7 ]->insert_neighbor( 18, tBasis[  20 ] );
                 mBasis[   7 ]->insert_neighbor( 19, mBasis[   1 ] );
                 mBasis[   7 ]->insert_neighbor( 20, tBasis[  26 ] );
                 mBasis[   7 ]->insert_neighbor( 21, tBasis[  24 ] );
                 mBasis[   7 ]->insert_neighbor( 22, tBasis[  44 ] );
                 mBasis[   7 ]->insert_neighbor( 23, tBasis[  46 ] );
                 mBasis[   7 ]->insert_neighbor( 24, tBasis[  54 ] );
                 mBasis[   7 ]->insert_neighbor( 25, tBasis[  52 ] );

                 // flag this basis
                 mBasis[   7 ]->flag();
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
    luint BSpline_Element< 3, 8 >::refine_basis( uint aBasisNumber )
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
                // create temporary container for children
                tBasis->init_children_container();

                // pointer to basis neighbor
                Basis* tNeighbor;

                // get neighbor 0
                tNeighbor = tBasis->get_neighbor( 0 );

                // test if neighbor 0 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 1
                tNeighbor = tBasis->get_neighbor( 1 );

                // test if neighbor 1 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  2, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 24 ) );
                    }
                }

                // get neighbor 2
                tNeighbor = tBasis->get_neighbor( 2 );

                // test if neighbor 2 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  6, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 20 ) );
                    }
                }

                // get neighbor 3
                tNeighbor = tBasis->get_neighbor( 3 );

                // test if neighbor 3 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 4
                tNeighbor = tBasis->get_neighbor( 4 );

                // test if neighbor 4 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 5
                tNeighbor = tBasis->get_neighbor( 5 );

                // test if neighbor 5 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 18, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child(  8 ) );
                    }
                }

                // get neighbor 6
                tNeighbor = tBasis->get_neighbor( 6 );

                // test if neighbor 6 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 7
                tNeighbor = tBasis->get_neighbor( 7 );

                // test if neighbor 7 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  2, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 24 ) );
                    }
                }

                // get neighbor 8
                tNeighbor = tBasis->get_neighbor( 8 );

                // test if neighbor 8 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  6, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 20 ) );
                    }
                }

                // get neighbor 9
                tNeighbor = tBasis->get_neighbor( 9 );

                // test if neighbor 9 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 10
                tNeighbor = tBasis->get_neighbor( 10 );

                // test if neighbor 10 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 11
                tNeighbor = tBasis->get_neighbor( 11 );

                // test if neighbor 11 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  2, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 24 ) );
                    }
                }

                // get neighbor 12
                tNeighbor = tBasis->get_neighbor( 12 );

                // test if neighbor 12 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  8, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 18 ) );
                    }
                }

                // get neighbor 13
                tNeighbor = tBasis->get_neighbor( 13 );

                // test if neighbor 13 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  6, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 20 ) );
                    }
                }

                // get neighbor 14
                tNeighbor = tBasis->get_neighbor( 14 );

                // test if neighbor 14 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 18, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child(  8 ) );
                    }
                }

                // get neighbor 15
                tNeighbor = tBasis->get_neighbor( 15 );

                // test if neighbor 15 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 20, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child(  6 ) );
                    }
                }

                // get neighbor 16
                tNeighbor = tBasis->get_neighbor( 16 );

                // test if neighbor 16 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 24, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child(  2 ) );
                    }
                }

                // get neighbor 17
                tNeighbor = tBasis->get_neighbor( 17 );

                // test if neighbor 17 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 18, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child(  8 ) );
                    }
                }

                // get neighbor 18
                tNeighbor = tBasis->get_neighbor( 18 );

                // test if neighbor 18 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  0, tNeighbor->get_child( 26 ) );
                    }
                }

                // get neighbor 19
                tNeighbor = tBasis->get_neighbor( 19 );

                // test if neighbor 19 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  2, tNeighbor->get_child( 24 ) );
                    }
                }

                // get neighbor 20
                tNeighbor = tBasis->get_neighbor( 20 );

                // test if neighbor 20 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  8, tNeighbor->get_child( 18 ) );
                    }
                }

                // get neighbor 21
                tNeighbor = tBasis->get_neighbor( 21 );

                // test if neighbor 21 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child(  6, tNeighbor->get_child( 20 ) );
                    }
                }

                // get neighbor 22
                tNeighbor = tBasis->get_neighbor( 22 );

                // test if neighbor 22 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 18, tNeighbor->get_child(  8 ) );
                    }
                }

                // get neighbor 23
                tNeighbor = tBasis->get_neighbor( 23 );

                // test if neighbor 23 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 20, tNeighbor->get_child(  6 ) );
                    }
                }

                // get neighbor 24
                tNeighbor = tBasis->get_neighbor( 24 );

                // test if neighbor 24 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 26, tNeighbor->get_child(  0 ) );
                    }
                }

                // get neighbor 25
                tNeighbor = tBasis->get_neighbor( 25 );

                // test if neighbor 25 exists
                if ( tNeighbor != nullptr )
                {
                    // test if neighbor has children
                    if ( tNeighbor->has_children() )
                    {
                        // copy children of neighbor
                        tBasis->insert_child( 24, tNeighbor->get_child(  2 ) );
                    }
                }

                // level of child basis
                uint tLevel = tBasis->get_level() + 1;

                // position of basis
                const luint* tParentIJK  = tBasis->get_ijk();

                // minumum i-position
                luint tIMin = 2*tParentIJK[ 0 ];

                // minumum j-position
                luint tJMin = 2*tParentIJK[ 1 ];

                // minumum k-position
                luint tKMin = 2*tParentIJK[ 2 ];

                // maximum i-position
                luint tIMax = tIMin + 3;

                // maximum j-position
                luint tJMax = tJMin + 3;

                // maximum K-position
                luint tKMax = tKMin + 3;

                // initialize counter
                uint tChildIndex = 0;

                // loop over all positions
                for( luint k=tKMin; k<tKMax; ++k )
                {
                    for( luint j=tJMin; j<tJMax; ++j )
                    {
                        for( luint i=tIMin; i<tIMax; ++i )
                        {
                            // test if child does not exist
                            if( tBasis->get_child( tChildIndex ) == nullptr )
                            {
                                 // calculate i-j-k position of child
                                 luint tIJK[ 3 ] = { i, j, k };

                                 // create child
                                 tBasis->insert_child( tChildIndex,
                                     new BSpline< 3, 27, 26 >( tIJK, tLevel, gNoProcOwner ) );

                                // increment basis counter
                                tBasisCounter++;
                            }

                            // increment child index
                            tChildIndex++;
                        }
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
    luint BSpline_Element< 3, 8 >::refine( moris::Cell< Element* > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;

        // refine basis if they have not been refined already
        for( uint k=0; k<8; ++k )
        {
            tBasisCounter += this->refine_basis( k );
        }

        // initialize temporary basis pattern
        Basis* tBasis[ 27 ] = { nullptr };

        // populate basis pattern
        if ( mBasis[   0 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   0 ]->get_child(  13 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   0 ]->get_child(  14 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   0 ]->get_child(  16 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   0 ]->get_child(  17 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   0 ]->get_child(  22 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   0 ]->get_child(  23 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   0 ]->get_child(  25 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   0 ]->get_child(  26 );
            }
        }

        if ( mBasis[   1 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   1 ]->get_child(  12 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   1 ]->get_child(  13 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   1 ]->get_child(  15 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   1 ]->get_child(  16 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   1 ]->get_child(  21 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   1 ]->get_child(  22 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   1 ]->get_child(  24 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   1 ]->get_child(  25 );
            }
        }

        if ( mBasis[   2 ] != nullptr )
        {
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   2 ]->get_child(   9 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   2 ]->get_child(  10 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   2 ]->get_child(  12 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   2 ]->get_child(  13 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   2 ]->get_child(  18 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   2 ]->get_child(  19 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   2 ]->get_child(  21 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   2 ]->get_child(  22 );
            }
        }

        if ( mBasis[   3 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   3 ]->get_child(  10 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   3 ]->get_child(  11 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   3 ]->get_child(  13 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   3 ]->get_child(  14 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   3 ]->get_child(  19 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   3 ]->get_child(  20 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   3 ]->get_child(  22 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   3 ]->get_child(  23 );
            }
        }

        if ( mBasis[   4 ] != nullptr )
        {
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   4 ]->get_child(   4 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   4 ]->get_child(   5 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   4 ]->get_child(   7 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   4 ]->get_child(   8 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   4 ]->get_child(  13 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   4 ]->get_child(  14 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   4 ]->get_child(  16 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   4 ]->get_child(  17 );
            }
        }

        if ( mBasis[   5 ] != nullptr )
        {
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   5 ]->get_child(   3 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   5 ]->get_child(   4 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   5 ]->get_child(   6 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   5 ]->get_child(   7 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   5 ]->get_child(  12 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   5 ]->get_child(  13 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   5 ]->get_child(  15 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   5 ]->get_child(  16 );
            }
        }

        if ( mBasis[   6 ] != nullptr )
        {
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   6 ]->get_child(   0 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   6 ]->get_child(   1 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   6 ]->get_child(   3 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   6 ]->get_child(   4 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   6 ]->get_child(   9 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   6 ]->get_child(  10 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[   6 ]->get_child(  12 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   6 ]->get_child(  13 );
            }
        }

        if ( mBasis[   7 ] != nullptr )
        {
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   7 ]->get_child(   1 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   7 ]->get_child(   2 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   7 ]->get_child(   4 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   7 ]->get_child(   5 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   7 ]->get_child(  10 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   7 ]->get_child(  11 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   7 ]->get_child(  13 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[   7 ]->get_child(  14 );
            }
        }

        // initialize child container
        Element * tChildren[ 8 ];

        // populate child container
        for( uint k=0; k<8; ++k)
        {
            tChildren[ k ] = aAllElementsOnProc(
                mElement->get_child( k )->get_memory_index() );
            tChildren[ k ]->init_basis_container();
        }

        // assign basis to child 1
        tChildren[ 0 ]->insert_basis(  0, tBasis[   0 ] );
        tChildren[ 0 ]->insert_basis(  1, tBasis[   1 ] );
        tChildren[ 0 ]->insert_basis(  2, tBasis[   4 ] );
        tChildren[ 0 ]->insert_basis(  3, tBasis[   3 ] );
        tChildren[ 0 ]->insert_basis(  4, tBasis[   9 ] );
        tChildren[ 0 ]->insert_basis(  5, tBasis[  10 ] );
        tChildren[ 0 ]->insert_basis(  6, tBasis[  13 ] );
        tChildren[ 0 ]->insert_basis(  7, tBasis[  12 ] );

        // assign basis to child 2
        tChildren[ 1 ]->insert_basis(  0, tBasis[   1 ] );
        tChildren[ 1 ]->insert_basis(  1, tBasis[   2 ] );
        tChildren[ 1 ]->insert_basis(  2, tBasis[   5 ] );
        tChildren[ 1 ]->insert_basis(  3, tBasis[   4 ] );
        tChildren[ 1 ]->insert_basis(  4, tBasis[  10 ] );
        tChildren[ 1 ]->insert_basis(  5, tBasis[  11 ] );
        tChildren[ 1 ]->insert_basis(  6, tBasis[  14 ] );
        tChildren[ 1 ]->insert_basis(  7, tBasis[  13 ] );

        // assign basis to child 3
        tChildren[ 2 ]->insert_basis(  0, tBasis[   3 ] );
        tChildren[ 2 ]->insert_basis(  1, tBasis[   4 ] );
        tChildren[ 2 ]->insert_basis(  2, tBasis[   7 ] );
        tChildren[ 2 ]->insert_basis(  3, tBasis[   6 ] );
        tChildren[ 2 ]->insert_basis(  4, tBasis[  12 ] );
        tChildren[ 2 ]->insert_basis(  5, tBasis[  13 ] );
        tChildren[ 2 ]->insert_basis(  6, tBasis[  16 ] );
        tChildren[ 2 ]->insert_basis(  7, tBasis[  15 ] );

        // assign basis to child 4
        tChildren[ 3 ]->insert_basis(  0, tBasis[   4 ] );
        tChildren[ 3 ]->insert_basis(  1, tBasis[   5 ] );
        tChildren[ 3 ]->insert_basis(  2, tBasis[   8 ] );
        tChildren[ 3 ]->insert_basis(  3, tBasis[   7 ] );
        tChildren[ 3 ]->insert_basis(  4, tBasis[  13 ] );
        tChildren[ 3 ]->insert_basis(  5, tBasis[  14 ] );
        tChildren[ 3 ]->insert_basis(  6, tBasis[  17 ] );
        tChildren[ 3 ]->insert_basis(  7, tBasis[  16 ] );

        // assign basis to child 5
        tChildren[ 4 ]->insert_basis(  0, tBasis[   9 ] );
        tChildren[ 4 ]->insert_basis(  1, tBasis[  10 ] );
        tChildren[ 4 ]->insert_basis(  2, tBasis[  13 ] );
        tChildren[ 4 ]->insert_basis(  3, tBasis[  12 ] );
        tChildren[ 4 ]->insert_basis(  4, tBasis[  18 ] );
        tChildren[ 4 ]->insert_basis(  5, tBasis[  19 ] );
        tChildren[ 4 ]->insert_basis(  6, tBasis[  22 ] );
        tChildren[ 4 ]->insert_basis(  7, tBasis[  21 ] );

        // assign basis to child 6
        tChildren[ 5 ]->insert_basis(  0, tBasis[  10 ] );
        tChildren[ 5 ]->insert_basis(  1, tBasis[  11 ] );
        tChildren[ 5 ]->insert_basis(  2, tBasis[  14 ] );
        tChildren[ 5 ]->insert_basis(  3, tBasis[  13 ] );
        tChildren[ 5 ]->insert_basis(  4, tBasis[  19 ] );
        tChildren[ 5 ]->insert_basis(  5, tBasis[  20 ] );
        tChildren[ 5 ]->insert_basis(  6, tBasis[  23 ] );
        tChildren[ 5 ]->insert_basis(  7, tBasis[  22 ] );

        // assign basis to child 7
        tChildren[ 6 ]->insert_basis(  0, tBasis[  12 ] );
        tChildren[ 6 ]->insert_basis(  1, tBasis[  13 ] );
        tChildren[ 6 ]->insert_basis(  2, tBasis[  16 ] );
        tChildren[ 6 ]->insert_basis(  3, tBasis[  15 ] );
        tChildren[ 6 ]->insert_basis(  4, tBasis[  21 ] );
        tChildren[ 6 ]->insert_basis(  5, tBasis[  22 ] );
        tChildren[ 6 ]->insert_basis(  6, tBasis[  25 ] );
        tChildren[ 6 ]->insert_basis(  7, tBasis[  24 ] );

        // assign basis to child 8
        tChildren[ 7 ]->insert_basis(  0, tBasis[  13 ] );
        tChildren[ 7 ]->insert_basis(  1, tBasis[  14 ] );
        tChildren[ 7 ]->insert_basis(  2, tBasis[  17 ] );
        tChildren[ 7 ]->insert_basis(  3, tBasis[  16 ] );
        tChildren[ 7 ]->insert_basis(  4, tBasis[  22 ] );
        tChildren[ 7 ]->insert_basis(  5, tBasis[  23 ] );
        tChildren[ 7 ]->insert_basis(  6, tBasis[  26 ] );
        tChildren[ 7 ]->insert_basis(  7, tBasis[  25 ] );

        // set basis flag of element
        mChildrenBasisFlag = true;
        
        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX8_HPP_ */

