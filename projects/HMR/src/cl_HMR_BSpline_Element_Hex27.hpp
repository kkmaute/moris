/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element_Hex27.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX27_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX27_HPP_

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
    BSpline_Element< 2, 2, 2 >::get_geometry_type() const
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
    BSpline_Element< 2, 2, 2 >::get_basis_indices_for_vtk(
        Matrix< DDLUMat > & aBasis )
    {
        // assemble nodes in correct order
       aBasis(  0 ) =  mBasis[  0 ]->get_memory_index();
       aBasis(  1 ) =  mBasis[  1 ]->get_memory_index();
       aBasis(  2 ) =  mBasis[  2 ]->get_memory_index();
       aBasis(  3 ) =  mBasis[  3 ]->get_memory_index();
       aBasis(  4 ) =  mBasis[  4 ]->get_memory_index();
       aBasis(  5 ) =  mBasis[  5 ]->get_memory_index();
       aBasis(  6 ) =  mBasis[  6 ]->get_memory_index();
       aBasis(  7 ) =  mBasis[  7 ]->get_memory_index();
       aBasis(  8 ) =  mBasis[  8 ]->get_memory_index();
       aBasis(  9 ) =  mBasis[ 11 ]->get_memory_index();
       aBasis( 10 ) =  mBasis[ 12 ]->get_memory_index();
       aBasis( 11 ) =  mBasis[  9 ]->get_memory_index();
       aBasis( 12 ) =  mBasis[ 13 ]->get_memory_index();
       aBasis( 13 ) =  mBasis[ 10 ]->get_memory_index();
       aBasis( 14 ) =  mBasis[ 14 ]->get_memory_index();
       aBasis( 15 ) =  mBasis[ 15 ]->get_memory_index();
       aBasis( 16 ) =  mBasis[ 16 ]->get_memory_index();
       aBasis( 17 ) =  mBasis[ 19 ]->get_memory_index();
       aBasis( 18 ) =  mBasis[ 17 ]->get_memory_index();
       aBasis( 19 ) =  mBasis[ 18 ]->get_memory_index();
       aBasis( 20 ) =  mBasis[ 21 ]->get_memory_index();
       aBasis( 21 ) =  mBasis[ 25 ]->get_memory_index();
       aBasis( 22 ) =  mBasis[ 23 ]->get_memory_index();
       aBasis( 23 ) =  mBasis[ 24 ]->get_memory_index();
       aBasis( 24 ) =  mBasis[ 26 ]->get_memory_index();
       aBasis( 25 ) =  mBasis[ 22 ]->get_memory_index();
       aBasis( 26 ) =  mBasis[ 20 ]->get_memory_index();
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
    BSpline_Element< 2, 2, 2 >::get_ijk_of_basis(
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
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  1 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  2 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  3 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  4 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case(  5 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case(  6 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case(  7 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case(  8 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  9 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 10 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 11 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 12 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 13 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 14 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 15 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 16 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 17 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 18 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 19 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 20 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 21 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 22 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 23 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 24 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 25 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 26 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
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
    BSpline_Element< 2, 2, 2 >::link_basis_with_neighbors(
          moris::Cell< Element* > & aAllElementsOnProc )
    {
         // initialize frame of basis around basis from this element
         Basis* tBasis[ 98 ] = { nullptr };

         // get pointer to neighbor  0
         Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

         // test if neighbor  0 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  26 ] = tNeighbor->get_basis(   0 );
             tBasis[  27 ] = tNeighbor->get_basis(   8 );
             tBasis[  28 ] = tNeighbor->get_basis(   1 );
             tBasis[  42 ] = tNeighbor->get_basis(  12 );
             tBasis[  43 ] = tNeighbor->get_basis(  25 );
             tBasis[  44 ] = tNeighbor->get_basis(  13 );
             tBasis[  58 ] = tNeighbor->get_basis(   4 );
             tBasis[  59 ] = tNeighbor->get_basis(  16 );
             tBasis[  60 ] = tNeighbor->get_basis(   5 );
         }

         // get pointer to neighbor  1
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

         // test if neighbor  1 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  31 ] = tNeighbor->get_basis(   1 );
             tBasis[  33 ] = tNeighbor->get_basis(   9 );
             tBasis[  35 ] = tNeighbor->get_basis(   2 );
             tBasis[  47 ] = tNeighbor->get_basis(  13 );
             tBasis[  49 ] = tNeighbor->get_basis(  24 );
             tBasis[  51 ] = tNeighbor->get_basis(  14 );
             tBasis[  63 ] = tNeighbor->get_basis(   5 );
             tBasis[  65 ] = tNeighbor->get_basis(  17 );
             tBasis[  67 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  2
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

         // test if neighbor  2 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  37 ] = tNeighbor->get_basis(   3 );
             tBasis[  38 ] = tNeighbor->get_basis(  10 );
             tBasis[  39 ] = tNeighbor->get_basis(   2 );
             tBasis[  53 ] = tNeighbor->get_basis(  15 );
             tBasis[  54 ] = tNeighbor->get_basis(  26 );
             tBasis[  55 ] = tNeighbor->get_basis(  14 );
             tBasis[  69 ] = tNeighbor->get_basis(   7 );
             tBasis[  70 ] = tNeighbor->get_basis(  18 );
             tBasis[  71 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  3
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

         // test if neighbor  3 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  30 ] = tNeighbor->get_basis(   0 );
             tBasis[  32 ] = tNeighbor->get_basis(  11 );
             tBasis[  34 ] = tNeighbor->get_basis(   3 );
             tBasis[  46 ] = tNeighbor->get_basis(  12 );
             tBasis[  48 ] = tNeighbor->get_basis(  23 );
             tBasis[  50 ] = tNeighbor->get_basis(  15 );
             tBasis[  62 ] = tNeighbor->get_basis(   4 );
             tBasis[  64 ] = tNeighbor->get_basis(  19 );
             tBasis[  66 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  4
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

         // test if neighbor  4 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   6 ] = tNeighbor->get_basis(   0 );
             tBasis[   7 ] = tNeighbor->get_basis(   8 );
             tBasis[   8 ] = tNeighbor->get_basis(   1 );
             tBasis[  11 ] = tNeighbor->get_basis(  11 );
             tBasis[  12 ] = tNeighbor->get_basis(  21 );
             tBasis[  13 ] = tNeighbor->get_basis(   9 );
             tBasis[  16 ] = tNeighbor->get_basis(   3 );
             tBasis[  17 ] = tNeighbor->get_basis(  10 );
             tBasis[  18 ] = tNeighbor->get_basis(   2 );
         }

         // get pointer to neighbor  5
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

         // test if neighbor  5 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  79 ] = tNeighbor->get_basis(   4 );
             tBasis[  80 ] = tNeighbor->get_basis(  16 );
             tBasis[  81 ] = tNeighbor->get_basis(   5 );
             tBasis[  84 ] = tNeighbor->get_basis(  19 );
             tBasis[  85 ] = tNeighbor->get_basis(  22 );
             tBasis[  86 ] = tNeighbor->get_basis(  17 );
             tBasis[  89 ] = tNeighbor->get_basis(   7 );
             tBasis[  90 ] = tNeighbor->get_basis(  18 );
             tBasis[  91 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  6
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

         // test if neighbor  6 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   1 ] = tNeighbor->get_basis(   0 );
             tBasis[   2 ] = tNeighbor->get_basis(   8 );
             tBasis[   3 ] = tNeighbor->get_basis(   1 );
             tBasis[   6 ] = tNeighbor->get_basis(  11 );
             tBasis[   7 ] = tNeighbor->get_basis(  21 );
             tBasis[   8 ] = tNeighbor->get_basis(   9 );
             tBasis[  11 ] = tNeighbor->get_basis(   3 );
             tBasis[  12 ] = tNeighbor->get_basis(  10 );
             tBasis[  13 ] = tNeighbor->get_basis(   2 );
             tBasis[  26 ] = tNeighbor->get_basis(  12 );
             tBasis[  27 ] = tNeighbor->get_basis(  25 );
             tBasis[  28 ] = tNeighbor->get_basis(  13 );
             tBasis[  42 ] = tNeighbor->get_basis(   4 );
             tBasis[  43 ] = tNeighbor->get_basis(  16 );
             tBasis[  44 ] = tNeighbor->get_basis(   5 );
         }

         // get pointer to neighbor  7
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

         // test if neighbor  7 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   7 ] = tNeighbor->get_basis(   0 );
             tBasis[   8 ] = tNeighbor->get_basis(   8 );
             tBasis[   9 ] = tNeighbor->get_basis(   1 );
             tBasis[  12 ] = tNeighbor->get_basis(  11 );
             tBasis[  13 ] = tNeighbor->get_basis(  21 );
             tBasis[  14 ] = tNeighbor->get_basis(   9 );
             tBasis[  17 ] = tNeighbor->get_basis(   3 );
             tBasis[  18 ] = tNeighbor->get_basis(  10 );
             tBasis[  19 ] = tNeighbor->get_basis(   2 );
             tBasis[  31 ] = tNeighbor->get_basis(  13 );
             tBasis[  33 ] = tNeighbor->get_basis(  24 );
             tBasis[  35 ] = tNeighbor->get_basis(  14 );
             tBasis[  47 ] = tNeighbor->get_basis(   5 );
             tBasis[  49 ] = tNeighbor->get_basis(  17 );
             tBasis[  51 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  8
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 8 );

         // test if neighbor  8 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  11 ] = tNeighbor->get_basis(   0 );
             tBasis[  12 ] = tNeighbor->get_basis(   8 );
             tBasis[  13 ] = tNeighbor->get_basis(   1 );
             tBasis[  16 ] = tNeighbor->get_basis(  11 );
             tBasis[  17 ] = tNeighbor->get_basis(  21 );
             tBasis[  18 ] = tNeighbor->get_basis(   9 );
             tBasis[  21 ] = tNeighbor->get_basis(   3 );
             tBasis[  22 ] = tNeighbor->get_basis(  10 );
             tBasis[  23 ] = tNeighbor->get_basis(   2 );
             tBasis[  37 ] = tNeighbor->get_basis(  15 );
             tBasis[  38 ] = tNeighbor->get_basis(  26 );
             tBasis[  39 ] = tNeighbor->get_basis(  14 );
             tBasis[  53 ] = tNeighbor->get_basis(   7 );
             tBasis[  54 ] = tNeighbor->get_basis(  18 );
             tBasis[  55 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  9
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 9 );

         // test if neighbor  9 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   5 ] = tNeighbor->get_basis(   0 );
             tBasis[   6 ] = tNeighbor->get_basis(   8 );
             tBasis[   7 ] = tNeighbor->get_basis(   1 );
             tBasis[  10 ] = tNeighbor->get_basis(  11 );
             tBasis[  11 ] = tNeighbor->get_basis(  21 );
             tBasis[  12 ] = tNeighbor->get_basis(   9 );
             tBasis[  15 ] = tNeighbor->get_basis(   3 );
             tBasis[  16 ] = tNeighbor->get_basis(  10 );
             tBasis[  17 ] = tNeighbor->get_basis(   2 );
             tBasis[  30 ] = tNeighbor->get_basis(  12 );
             tBasis[  32 ] = tNeighbor->get_basis(  23 );
             tBasis[  34 ] = tNeighbor->get_basis(  15 );
             tBasis[  46 ] = tNeighbor->get_basis(   4 );
             tBasis[  48 ] = tNeighbor->get_basis(  19 );
             tBasis[  50 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  10
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 10 );

         // test if neighbor  10 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  25 ] = tNeighbor->get_basis(   0 );
             tBasis[  26 ] = tNeighbor->get_basis(   8 );
             tBasis[  27 ] = tNeighbor->get_basis(   1 );
             tBasis[  30 ] = tNeighbor->get_basis(  11 );
             tBasis[  32 ] = tNeighbor->get_basis(   3 );
             tBasis[  41 ] = tNeighbor->get_basis(  12 );
             tBasis[  42 ] = tNeighbor->get_basis(  25 );
             tBasis[  43 ] = tNeighbor->get_basis(  13 );
             tBasis[  46 ] = tNeighbor->get_basis(  23 );
             tBasis[  48 ] = tNeighbor->get_basis(  15 );
             tBasis[  57 ] = tNeighbor->get_basis(   4 );
             tBasis[  58 ] = tNeighbor->get_basis(  16 );
             tBasis[  59 ] = tNeighbor->get_basis(   5 );
             tBasis[  62 ] = tNeighbor->get_basis(  19 );
             tBasis[  64 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  11
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 11 );

         // test if neighbor  11 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  27 ] = tNeighbor->get_basis(   0 );
             tBasis[  28 ] = tNeighbor->get_basis(   8 );
             tBasis[  29 ] = tNeighbor->get_basis(   1 );
             tBasis[  31 ] = tNeighbor->get_basis(   9 );
             tBasis[  33 ] = tNeighbor->get_basis(   2 );
             tBasis[  43 ] = tNeighbor->get_basis(  12 );
             tBasis[  44 ] = tNeighbor->get_basis(  25 );
             tBasis[  45 ] = tNeighbor->get_basis(  13 );
             tBasis[  47 ] = tNeighbor->get_basis(  24 );
             tBasis[  49 ] = tNeighbor->get_basis(  14 );
             tBasis[  59 ] = tNeighbor->get_basis(   4 );
             tBasis[  60 ] = tNeighbor->get_basis(  16 );
             tBasis[  61 ] = tNeighbor->get_basis(   5 );
             tBasis[  63 ] = tNeighbor->get_basis(  17 );
             tBasis[  65 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  12
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 12 );

         // test if neighbor  12 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  33 ] = tNeighbor->get_basis(   1 );
             tBasis[  35 ] = tNeighbor->get_basis(   9 );
             tBasis[  38 ] = tNeighbor->get_basis(   3 );
             tBasis[  39 ] = tNeighbor->get_basis(  10 );
             tBasis[  40 ] = tNeighbor->get_basis(   2 );
             tBasis[  49 ] = tNeighbor->get_basis(  13 );
             tBasis[  51 ] = tNeighbor->get_basis(  24 );
             tBasis[  54 ] = tNeighbor->get_basis(  15 );
             tBasis[  55 ] = tNeighbor->get_basis(  26 );
             tBasis[  56 ] = tNeighbor->get_basis(  14 );
             tBasis[  65 ] = tNeighbor->get_basis(   5 );
             tBasis[  67 ] = tNeighbor->get_basis(  17 );
             tBasis[  70 ] = tNeighbor->get_basis(   7 );
             tBasis[  71 ] = tNeighbor->get_basis(  18 );
             tBasis[  72 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  13
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 13 );

         // test if neighbor  13 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  32 ] = tNeighbor->get_basis(   0 );
             tBasis[  34 ] = tNeighbor->get_basis(  11 );
             tBasis[  36 ] = tNeighbor->get_basis(   3 );
             tBasis[  37 ] = tNeighbor->get_basis(  10 );
             tBasis[  38 ] = tNeighbor->get_basis(   2 );
             tBasis[  48 ] = tNeighbor->get_basis(  12 );
             tBasis[  50 ] = tNeighbor->get_basis(  23 );
             tBasis[  52 ] = tNeighbor->get_basis(  15 );
             tBasis[  53 ] = tNeighbor->get_basis(  26 );
             tBasis[  54 ] = tNeighbor->get_basis(  14 );
             tBasis[  64 ] = tNeighbor->get_basis(   4 );
             tBasis[  66 ] = tNeighbor->get_basis(  19 );
             tBasis[  68 ] = tNeighbor->get_basis(   7 );
             tBasis[  69 ] = tNeighbor->get_basis(  18 );
             tBasis[  70 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  14
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 14 );

         // test if neighbor  14 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  42 ] = tNeighbor->get_basis(   0 );
             tBasis[  43 ] = tNeighbor->get_basis(   8 );
             tBasis[  44 ] = tNeighbor->get_basis(   1 );
             tBasis[  58 ] = tNeighbor->get_basis(  12 );
             tBasis[  59 ] = tNeighbor->get_basis(  25 );
             tBasis[  60 ] = tNeighbor->get_basis(  13 );
             tBasis[  74 ] = tNeighbor->get_basis(   4 );
             tBasis[  75 ] = tNeighbor->get_basis(  16 );
             tBasis[  76 ] = tNeighbor->get_basis(   5 );
             tBasis[  79 ] = tNeighbor->get_basis(  19 );
             tBasis[  80 ] = tNeighbor->get_basis(  22 );
             tBasis[  81 ] = tNeighbor->get_basis(  17 );
             tBasis[  84 ] = tNeighbor->get_basis(   7 );
             tBasis[  85 ] = tNeighbor->get_basis(  18 );
             tBasis[  86 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  15
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 15 );

         // test if neighbor  15 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  47 ] = tNeighbor->get_basis(   1 );
             tBasis[  49 ] = tNeighbor->get_basis(   9 );
             tBasis[  51 ] = tNeighbor->get_basis(   2 );
             tBasis[  63 ] = tNeighbor->get_basis(  13 );
             tBasis[  65 ] = tNeighbor->get_basis(  24 );
             tBasis[  67 ] = tNeighbor->get_basis(  14 );
             tBasis[  80 ] = tNeighbor->get_basis(   4 );
             tBasis[  81 ] = tNeighbor->get_basis(  16 );
             tBasis[  82 ] = tNeighbor->get_basis(   5 );
             tBasis[  85 ] = tNeighbor->get_basis(  19 );
             tBasis[  86 ] = tNeighbor->get_basis(  22 );
             tBasis[  87 ] = tNeighbor->get_basis(  17 );
             tBasis[  90 ] = tNeighbor->get_basis(   7 );
             tBasis[  91 ] = tNeighbor->get_basis(  18 );
             tBasis[  92 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  16
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 16 );

         // test if neighbor  16 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  53 ] = tNeighbor->get_basis(   3 );
             tBasis[  54 ] = tNeighbor->get_basis(  10 );
             tBasis[  55 ] = tNeighbor->get_basis(   2 );
             tBasis[  69 ] = tNeighbor->get_basis(  15 );
             tBasis[  70 ] = tNeighbor->get_basis(  26 );
             tBasis[  71 ] = tNeighbor->get_basis(  14 );
             tBasis[  84 ] = tNeighbor->get_basis(   4 );
             tBasis[  85 ] = tNeighbor->get_basis(  16 );
             tBasis[  86 ] = tNeighbor->get_basis(   5 );
             tBasis[  89 ] = tNeighbor->get_basis(  19 );
             tBasis[  90 ] = tNeighbor->get_basis(  22 );
             tBasis[  91 ] = tNeighbor->get_basis(  17 );
             tBasis[  94 ] = tNeighbor->get_basis(   7 );
             tBasis[  95 ] = tNeighbor->get_basis(  18 );
             tBasis[  96 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  17
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 17 );

         // test if neighbor  17 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  46 ] = tNeighbor->get_basis(   0 );
             tBasis[  48 ] = tNeighbor->get_basis(  11 );
             tBasis[  50 ] = tNeighbor->get_basis(   3 );
             tBasis[  62 ] = tNeighbor->get_basis(  12 );
             tBasis[  64 ] = tNeighbor->get_basis(  23 );
             tBasis[  66 ] = tNeighbor->get_basis(  15 );
             tBasis[  78 ] = tNeighbor->get_basis(   4 );
             tBasis[  79 ] = tNeighbor->get_basis(  16 );
             tBasis[  80 ] = tNeighbor->get_basis(   5 );
             tBasis[  83 ] = tNeighbor->get_basis(  19 );
             tBasis[  84 ] = tNeighbor->get_basis(  22 );
             tBasis[  85 ] = tNeighbor->get_basis(  17 );
             tBasis[  88 ] = tNeighbor->get_basis(   7 );
             tBasis[  89 ] = tNeighbor->get_basis(  18 );
             tBasis[  90 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  18
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 18 );

         // test if neighbor  18 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   0 ] = tNeighbor->get_basis(   0 );
             tBasis[   1 ] = tNeighbor->get_basis(   8 );
             tBasis[   2 ] = tNeighbor->get_basis(   1 );
             tBasis[   5 ] = tNeighbor->get_basis(  11 );
             tBasis[   6 ] = tNeighbor->get_basis(  21 );
             tBasis[   7 ] = tNeighbor->get_basis(   9 );
             tBasis[  10 ] = tNeighbor->get_basis(   3 );
             tBasis[  11 ] = tNeighbor->get_basis(  10 );
             tBasis[  12 ] = tNeighbor->get_basis(   2 );
             tBasis[  25 ] = tNeighbor->get_basis(  12 );
             tBasis[  26 ] = tNeighbor->get_basis(  25 );
             tBasis[  27 ] = tNeighbor->get_basis(  13 );
             tBasis[  30 ] = tNeighbor->get_basis(  23 );
             tBasis[  32 ] = tNeighbor->get_basis(  15 );
             tBasis[  41 ] = tNeighbor->get_basis(   4 );
             tBasis[  42 ] = tNeighbor->get_basis(  16 );
             tBasis[  43 ] = tNeighbor->get_basis(   5 );
             tBasis[  46 ] = tNeighbor->get_basis(  19 );
             tBasis[  48 ] = tNeighbor->get_basis(   7 );
         }

         // get pointer to neighbor  19
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 19 );

         // test if neighbor  19 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   2 ] = tNeighbor->get_basis(   0 );
             tBasis[   3 ] = tNeighbor->get_basis(   8 );
             tBasis[   4 ] = tNeighbor->get_basis(   1 );
             tBasis[   7 ] = tNeighbor->get_basis(  11 );
             tBasis[   8 ] = tNeighbor->get_basis(  21 );
             tBasis[   9 ] = tNeighbor->get_basis(   9 );
             tBasis[  12 ] = tNeighbor->get_basis(   3 );
             tBasis[  13 ] = tNeighbor->get_basis(  10 );
             tBasis[  14 ] = tNeighbor->get_basis(   2 );
             tBasis[  27 ] = tNeighbor->get_basis(  12 );
             tBasis[  28 ] = tNeighbor->get_basis(  25 );
             tBasis[  29 ] = tNeighbor->get_basis(  13 );
             tBasis[  31 ] = tNeighbor->get_basis(  24 );
             tBasis[  33 ] = tNeighbor->get_basis(  14 );
             tBasis[  43 ] = tNeighbor->get_basis(   4 );
             tBasis[  44 ] = tNeighbor->get_basis(  16 );
             tBasis[  45 ] = tNeighbor->get_basis(   5 );
             tBasis[  47 ] = tNeighbor->get_basis(  17 );
             tBasis[  49 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  20
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 20 );

         // test if neighbor  20 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  12 ] = tNeighbor->get_basis(   0 );
             tBasis[  13 ] = tNeighbor->get_basis(   8 );
             tBasis[  14 ] = tNeighbor->get_basis(   1 );
             tBasis[  17 ] = tNeighbor->get_basis(  11 );
             tBasis[  18 ] = tNeighbor->get_basis(  21 );
             tBasis[  19 ] = tNeighbor->get_basis(   9 );
             tBasis[  22 ] = tNeighbor->get_basis(   3 );
             tBasis[  23 ] = tNeighbor->get_basis(  10 );
             tBasis[  24 ] = tNeighbor->get_basis(   2 );
             tBasis[  33 ] = tNeighbor->get_basis(  13 );
             tBasis[  35 ] = tNeighbor->get_basis(  24 );
             tBasis[  38 ] = tNeighbor->get_basis(  15 );
             tBasis[  39 ] = tNeighbor->get_basis(  26 );
             tBasis[  40 ] = tNeighbor->get_basis(  14 );
             tBasis[  49 ] = tNeighbor->get_basis(   5 );
             tBasis[  51 ] = tNeighbor->get_basis(  17 );
             tBasis[  54 ] = tNeighbor->get_basis(   7 );
             tBasis[  55 ] = tNeighbor->get_basis(  18 );
             tBasis[  56 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  21
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 21 );

         // test if neighbor  21 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  10 ] = tNeighbor->get_basis(   0 );
             tBasis[  11 ] = tNeighbor->get_basis(   8 );
             tBasis[  12 ] = tNeighbor->get_basis(   1 );
             tBasis[  15 ] = tNeighbor->get_basis(  11 );
             tBasis[  16 ] = tNeighbor->get_basis(  21 );
             tBasis[  17 ] = tNeighbor->get_basis(   9 );
             tBasis[  20 ] = tNeighbor->get_basis(   3 );
             tBasis[  21 ] = tNeighbor->get_basis(  10 );
             tBasis[  22 ] = tNeighbor->get_basis(   2 );
             tBasis[  32 ] = tNeighbor->get_basis(  12 );
             tBasis[  34 ] = tNeighbor->get_basis(  23 );
             tBasis[  36 ] = tNeighbor->get_basis(  15 );
             tBasis[  37 ] = tNeighbor->get_basis(  26 );
             tBasis[  38 ] = tNeighbor->get_basis(  14 );
             tBasis[  48 ] = tNeighbor->get_basis(   4 );
             tBasis[  50 ] = tNeighbor->get_basis(  19 );
             tBasis[  52 ] = tNeighbor->get_basis(   7 );
             tBasis[  53 ] = tNeighbor->get_basis(  18 );
             tBasis[  54 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  22
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 22 );

         // test if neighbor  22 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  41 ] = tNeighbor->get_basis(   0 );
             tBasis[  42 ] = tNeighbor->get_basis(   8 );
             tBasis[  43 ] = tNeighbor->get_basis(   1 );
             tBasis[  46 ] = tNeighbor->get_basis(  11 );
             tBasis[  48 ] = tNeighbor->get_basis(   3 );
             tBasis[  57 ] = tNeighbor->get_basis(  12 );
             tBasis[  58 ] = tNeighbor->get_basis(  25 );
             tBasis[  59 ] = tNeighbor->get_basis(  13 );
             tBasis[  62 ] = tNeighbor->get_basis(  23 );
             tBasis[  64 ] = tNeighbor->get_basis(  15 );
             tBasis[  73 ] = tNeighbor->get_basis(   4 );
             tBasis[  74 ] = tNeighbor->get_basis(  16 );
             tBasis[  75 ] = tNeighbor->get_basis(   5 );
             tBasis[  78 ] = tNeighbor->get_basis(  19 );
             tBasis[  79 ] = tNeighbor->get_basis(  22 );
             tBasis[  80 ] = tNeighbor->get_basis(  17 );
             tBasis[  83 ] = tNeighbor->get_basis(   7 );
             tBasis[  84 ] = tNeighbor->get_basis(  18 );
             tBasis[  85 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  23
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 23 );

         // test if neighbor  23 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  43 ] = tNeighbor->get_basis(   0 );
             tBasis[  44 ] = tNeighbor->get_basis(   8 );
             tBasis[  45 ] = tNeighbor->get_basis(   1 );
             tBasis[  47 ] = tNeighbor->get_basis(   9 );
             tBasis[  49 ] = tNeighbor->get_basis(   2 );
             tBasis[  59 ] = tNeighbor->get_basis(  12 );
             tBasis[  60 ] = tNeighbor->get_basis(  25 );
             tBasis[  61 ] = tNeighbor->get_basis(  13 );
             tBasis[  63 ] = tNeighbor->get_basis(  24 );
             tBasis[  65 ] = tNeighbor->get_basis(  14 );
             tBasis[  75 ] = tNeighbor->get_basis(   4 );
             tBasis[  76 ] = tNeighbor->get_basis(  16 );
             tBasis[  77 ] = tNeighbor->get_basis(   5 );
             tBasis[  80 ] = tNeighbor->get_basis(  19 );
             tBasis[  81 ] = tNeighbor->get_basis(  22 );
             tBasis[  82 ] = tNeighbor->get_basis(  17 );
             tBasis[  85 ] = tNeighbor->get_basis(   7 );
             tBasis[  86 ] = tNeighbor->get_basis(  18 );
             tBasis[  87 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  24
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 24 );

         // test if neighbor  24 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  49 ] = tNeighbor->get_basis(   1 );
             tBasis[  51 ] = tNeighbor->get_basis(   9 );
             tBasis[  54 ] = tNeighbor->get_basis(   3 );
             tBasis[  55 ] = tNeighbor->get_basis(  10 );
             tBasis[  56 ] = tNeighbor->get_basis(   2 );
             tBasis[  65 ] = tNeighbor->get_basis(  13 );
             tBasis[  67 ] = tNeighbor->get_basis(  24 );
             tBasis[  70 ] = tNeighbor->get_basis(  15 );
             tBasis[  71 ] = tNeighbor->get_basis(  26 );
             tBasis[  72 ] = tNeighbor->get_basis(  14 );
             tBasis[  85 ] = tNeighbor->get_basis(   4 );
             tBasis[  86 ] = tNeighbor->get_basis(  16 );
             tBasis[  87 ] = tNeighbor->get_basis(   5 );
             tBasis[  90 ] = tNeighbor->get_basis(  19 );
             tBasis[  91 ] = tNeighbor->get_basis(  22 );
             tBasis[  92 ] = tNeighbor->get_basis(  17 );
             tBasis[  95 ] = tNeighbor->get_basis(   7 );
             tBasis[  96 ] = tNeighbor->get_basis(  18 );
             tBasis[  97 ] = tNeighbor->get_basis(   6 );
         }

         // get pointer to neighbor  25
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 25 );

         // test if neighbor  25 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  48 ] = tNeighbor->get_basis(   0 );
             tBasis[  50 ] = tNeighbor->get_basis(  11 );
             tBasis[  52 ] = tNeighbor->get_basis(   3 );
             tBasis[  53 ] = tNeighbor->get_basis(  10 );
             tBasis[  54 ] = tNeighbor->get_basis(   2 );
             tBasis[  64 ] = tNeighbor->get_basis(  12 );
             tBasis[  66 ] = tNeighbor->get_basis(  23 );
             tBasis[  68 ] = tNeighbor->get_basis(  15 );
             tBasis[  69 ] = tNeighbor->get_basis(  26 );
             tBasis[  70 ] = tNeighbor->get_basis(  14 );
             tBasis[  83 ] = tNeighbor->get_basis(   4 );
             tBasis[  84 ] = tNeighbor->get_basis(  16 );
             tBasis[  85 ] = tNeighbor->get_basis(   5 );
             tBasis[  88 ] = tNeighbor->get_basis(  19 );
             tBasis[  89 ] = tNeighbor->get_basis(  22 );
             tBasis[  90 ] = tNeighbor->get_basis(  17 );
             tBasis[  93 ] = tNeighbor->get_basis(   7 );
             tBasis[  94 ] = tNeighbor->get_basis(  18 );
             tBasis[  95 ] = tNeighbor->get_basis(   6 );
         }

         // test if basis 0 exists
         if ( mBasis[   0 ] != nullptr )
         {
             // test if basis 0 has been processed
             if ( ! mBasis[   0 ]->is_flagged() )
             {
                 // link neighbors of basis 0
                 mBasis[   0 ]->insert_neighbor(  0, tBasis[  26 ] );
                 mBasis[   0 ]->insert_neighbor(  1, mBasis[   8 ] );
                 mBasis[   0 ]->insert_neighbor(  2, mBasis[  11 ] );
                 mBasis[   0 ]->insert_neighbor(  3, tBasis[  30 ] );
                 mBasis[   0 ]->insert_neighbor(  4, tBasis[   6 ] );
                 mBasis[   0 ]->insert_neighbor(  5, mBasis[  12 ] );
                 mBasis[   0 ]->insert_neighbor(  6, tBasis[   1 ] );
                 mBasis[   0 ]->insert_neighbor(  7, tBasis[   7 ] );
                 mBasis[   0 ]->insert_neighbor(  8, tBasis[  11 ] );
                 mBasis[   0 ]->insert_neighbor(  9, tBasis[   5 ] );
                 mBasis[   0 ]->insert_neighbor( 10, tBasis[  25 ] );
                 mBasis[   0 ]->insert_neighbor( 11, tBasis[  27 ] );
                 mBasis[   0 ]->insert_neighbor( 12, mBasis[  21 ] );
                 mBasis[   0 ]->insert_neighbor( 13, tBasis[  32 ] );
                 mBasis[   0 ]->insert_neighbor( 14, tBasis[  42 ] );
                 mBasis[   0 ]->insert_neighbor( 15, mBasis[  25 ] );
                 mBasis[   0 ]->insert_neighbor( 16, mBasis[  23 ] );
                 mBasis[   0 ]->insert_neighbor( 17, tBasis[  46 ] );
                 mBasis[   0 ]->insert_neighbor( 18, tBasis[   0 ] );
                 mBasis[   0 ]->insert_neighbor( 19, tBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor( 20, tBasis[  12 ] );
                 mBasis[   0 ]->insert_neighbor( 21, tBasis[  10 ] );
                 mBasis[   0 ]->insert_neighbor( 22, tBasis[  41 ] );
                 mBasis[   0 ]->insert_neighbor( 23, tBasis[  43 ] );
                 mBasis[   0 ]->insert_neighbor( 24, mBasis[  20 ] );
                 mBasis[   0 ]->insert_neighbor( 25, tBasis[  48 ] );

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
                 mBasis[   1 ]->insert_neighbor(  0, tBasis[  28 ] );
                 mBasis[   1 ]->insert_neighbor(  1, tBasis[  31 ] );
                 mBasis[   1 ]->insert_neighbor(  2, mBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor(  3, mBasis[   8 ] );
                 mBasis[   1 ]->insert_neighbor(  4, tBasis[   8 ] );
                 mBasis[   1 ]->insert_neighbor(  5, mBasis[  13 ] );
                 mBasis[   1 ]->insert_neighbor(  6, tBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor(  7, tBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor(  8, tBasis[  13 ] );
                 mBasis[   1 ]->insert_neighbor(  9, tBasis[   7 ] );
                 mBasis[   1 ]->insert_neighbor( 10, tBasis[  27 ] );
                 mBasis[   1 ]->insert_neighbor( 11, tBasis[  29 ] );
                 mBasis[   1 ]->insert_neighbor( 12, tBasis[  33 ] );
                 mBasis[   1 ]->insert_neighbor( 13, mBasis[  21 ] );
                 mBasis[   1 ]->insert_neighbor( 14, tBasis[  44 ] );
                 mBasis[   1 ]->insert_neighbor( 15, tBasis[  47 ] );
                 mBasis[   1 ]->insert_neighbor( 16, mBasis[  24 ] );
                 mBasis[   1 ]->insert_neighbor( 17, mBasis[  25 ] );
                 mBasis[   1 ]->insert_neighbor( 18, tBasis[   2 ] );
                 mBasis[   1 ]->insert_neighbor( 19, tBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor( 20, tBasis[  14 ] );
                 mBasis[   1 ]->insert_neighbor( 21, tBasis[  12 ] );
                 mBasis[   1 ]->insert_neighbor( 22, tBasis[  43 ] );
                 mBasis[   1 ]->insert_neighbor( 23, tBasis[  45 ] );
                 mBasis[   1 ]->insert_neighbor( 24, tBasis[  49 ] );
                 mBasis[   1 ]->insert_neighbor( 25, mBasis[  20 ] );

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
                 mBasis[   2 ]->insert_neighbor(  0, mBasis[   9 ] );
                 mBasis[   2 ]->insert_neighbor(  1, tBasis[  35 ] );
                 mBasis[   2 ]->insert_neighbor(  2, tBasis[  39 ] );
                 mBasis[   2 ]->insert_neighbor(  3, mBasis[  10 ] );
                 mBasis[   2 ]->insert_neighbor(  4, tBasis[  18 ] );
                 mBasis[   2 ]->insert_neighbor(  5, mBasis[  14 ] );
                 mBasis[   2 ]->insert_neighbor(  6, tBasis[  13 ] );
                 mBasis[   2 ]->insert_neighbor(  7, tBasis[  19 ] );
                 mBasis[   2 ]->insert_neighbor(  8, tBasis[  23 ] );
                 mBasis[   2 ]->insert_neighbor(  9, tBasis[  17 ] );
                 mBasis[   2 ]->insert_neighbor( 10, mBasis[  21 ] );
                 mBasis[   2 ]->insert_neighbor( 11, tBasis[  33 ] );
                 mBasis[   2 ]->insert_neighbor( 12, tBasis[  40 ] );
                 mBasis[   2 ]->insert_neighbor( 13, tBasis[  38 ] );
                 mBasis[   2 ]->insert_neighbor( 14, mBasis[  24 ] );
                 mBasis[   2 ]->insert_neighbor( 15, tBasis[  51 ] );
                 mBasis[   2 ]->insert_neighbor( 16, tBasis[  55 ] );
                 mBasis[   2 ]->insert_neighbor( 17, mBasis[  26 ] );
                 mBasis[   2 ]->insert_neighbor( 18, tBasis[  12 ] );
                 mBasis[   2 ]->insert_neighbor( 19, tBasis[  14 ] );
                 mBasis[   2 ]->insert_neighbor( 20, tBasis[  24 ] );
                 mBasis[   2 ]->insert_neighbor( 21, tBasis[  22 ] );
                 mBasis[   2 ]->insert_neighbor( 22, mBasis[  20 ] );
                 mBasis[   2 ]->insert_neighbor( 23, tBasis[  49 ] );
                 mBasis[   2 ]->insert_neighbor( 24, tBasis[  56 ] );
                 mBasis[   2 ]->insert_neighbor( 25, tBasis[  54 ] );

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
                 mBasis[   3 ]->insert_neighbor(  0, mBasis[  11 ] );
                 mBasis[   3 ]->insert_neighbor(  1, mBasis[  10 ] );
                 mBasis[   3 ]->insert_neighbor(  2, tBasis[  37 ] );
                 mBasis[   3 ]->insert_neighbor(  3, tBasis[  34 ] );
                 mBasis[   3 ]->insert_neighbor(  4, tBasis[  16 ] );
                 mBasis[   3 ]->insert_neighbor(  5, mBasis[  15 ] );
                 mBasis[   3 ]->insert_neighbor(  6, tBasis[  11 ] );
                 mBasis[   3 ]->insert_neighbor(  7, tBasis[  17 ] );
                 mBasis[   3 ]->insert_neighbor(  8, tBasis[  21 ] );
                 mBasis[   3 ]->insert_neighbor(  9, tBasis[  15 ] );
                 mBasis[   3 ]->insert_neighbor( 10, tBasis[  32 ] );
                 mBasis[   3 ]->insert_neighbor( 11, mBasis[  21 ] );
                 mBasis[   3 ]->insert_neighbor( 12, tBasis[  38 ] );
                 mBasis[   3 ]->insert_neighbor( 13, tBasis[  36 ] );
                 mBasis[   3 ]->insert_neighbor( 14, mBasis[  23 ] );
                 mBasis[   3 ]->insert_neighbor( 15, mBasis[  26 ] );
                 mBasis[   3 ]->insert_neighbor( 16, tBasis[  53 ] );
                 mBasis[   3 ]->insert_neighbor( 17, tBasis[  50 ] );
                 mBasis[   3 ]->insert_neighbor( 18, tBasis[  10 ] );
                 mBasis[   3 ]->insert_neighbor( 19, tBasis[  12 ] );
                 mBasis[   3 ]->insert_neighbor( 20, tBasis[  22 ] );
                 mBasis[   3 ]->insert_neighbor( 21, tBasis[  20 ] );
                 mBasis[   3 ]->insert_neighbor( 22, tBasis[  48 ] );
                 mBasis[   3 ]->insert_neighbor( 23, mBasis[  20 ] );
                 mBasis[   3 ]->insert_neighbor( 24, tBasis[  54 ] );
                 mBasis[   3 ]->insert_neighbor( 25, tBasis[  52 ] );

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
                 mBasis[   4 ]->insert_neighbor(  0, tBasis[  58 ] );
                 mBasis[   4 ]->insert_neighbor(  1, mBasis[  16 ] );
                 mBasis[   4 ]->insert_neighbor(  2, mBasis[  19 ] );
                 mBasis[   4 ]->insert_neighbor(  3, tBasis[  62 ] );
                 mBasis[   4 ]->insert_neighbor(  4, mBasis[  12 ] );
                 mBasis[   4 ]->insert_neighbor(  5, tBasis[  79 ] );
                 mBasis[   4 ]->insert_neighbor(  6, tBasis[  42 ] );
                 mBasis[   4 ]->insert_neighbor(  7, mBasis[  25 ] );
                 mBasis[   4 ]->insert_neighbor(  8, mBasis[  23 ] );
                 mBasis[   4 ]->insert_neighbor(  9, tBasis[  46 ] );
                 mBasis[   4 ]->insert_neighbor( 10, tBasis[  57 ] );
                 mBasis[   4 ]->insert_neighbor( 11, tBasis[  59 ] );
                 mBasis[   4 ]->insert_neighbor( 12, mBasis[  22 ] );
                 mBasis[   4 ]->insert_neighbor( 13, tBasis[  64 ] );
                 mBasis[   4 ]->insert_neighbor( 14, tBasis[  74 ] );
                 mBasis[   4 ]->insert_neighbor( 15, tBasis[  80 ] );
                 mBasis[   4 ]->insert_neighbor( 16, tBasis[  84 ] );
                 mBasis[   4 ]->insert_neighbor( 17, tBasis[  78 ] );
                 mBasis[   4 ]->insert_neighbor( 18, tBasis[  41 ] );
                 mBasis[   4 ]->insert_neighbor( 19, tBasis[  43 ] );
                 mBasis[   4 ]->insert_neighbor( 20, mBasis[  20 ] );
                 mBasis[   4 ]->insert_neighbor( 21, tBasis[  48 ] );
                 mBasis[   4 ]->insert_neighbor( 22, tBasis[  73 ] );
                 mBasis[   4 ]->insert_neighbor( 23, tBasis[  75 ] );
                 mBasis[   4 ]->insert_neighbor( 24, tBasis[  85 ] );
                 mBasis[   4 ]->insert_neighbor( 25, tBasis[  83 ] );

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
                 mBasis[   5 ]->insert_neighbor(  0, tBasis[  60 ] );
                 mBasis[   5 ]->insert_neighbor(  1, tBasis[  63 ] );
                 mBasis[   5 ]->insert_neighbor(  2, mBasis[  17 ] );
                 mBasis[   5 ]->insert_neighbor(  3, mBasis[  16 ] );
                 mBasis[   5 ]->insert_neighbor(  4, mBasis[  13 ] );
                 mBasis[   5 ]->insert_neighbor(  5, tBasis[  81 ] );
                 mBasis[   5 ]->insert_neighbor(  6, tBasis[  44 ] );
                 mBasis[   5 ]->insert_neighbor(  7, tBasis[  47 ] );
                 mBasis[   5 ]->insert_neighbor(  8, mBasis[  24 ] );
                 mBasis[   5 ]->insert_neighbor(  9, mBasis[  25 ] );
                 mBasis[   5 ]->insert_neighbor( 10, tBasis[  59 ] );
                 mBasis[   5 ]->insert_neighbor( 11, tBasis[  61 ] );
                 mBasis[   5 ]->insert_neighbor( 12, tBasis[  65 ] );
                 mBasis[   5 ]->insert_neighbor( 13, mBasis[  22 ] );
                 mBasis[   5 ]->insert_neighbor( 14, tBasis[  76 ] );
                 mBasis[   5 ]->insert_neighbor( 15, tBasis[  82 ] );
                 mBasis[   5 ]->insert_neighbor( 16, tBasis[  86 ] );
                 mBasis[   5 ]->insert_neighbor( 17, tBasis[  80 ] );
                 mBasis[   5 ]->insert_neighbor( 18, tBasis[  43 ] );
                 mBasis[   5 ]->insert_neighbor( 19, tBasis[  45 ] );
                 mBasis[   5 ]->insert_neighbor( 20, tBasis[  49 ] );
                 mBasis[   5 ]->insert_neighbor( 21, mBasis[  20 ] );
                 mBasis[   5 ]->insert_neighbor( 22, tBasis[  75 ] );
                 mBasis[   5 ]->insert_neighbor( 23, tBasis[  77 ] );
                 mBasis[   5 ]->insert_neighbor( 24, tBasis[  87 ] );
                 mBasis[   5 ]->insert_neighbor( 25, tBasis[  85 ] );

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
                 mBasis[   6 ]->insert_neighbor(  0, mBasis[  17 ] );
                 mBasis[   6 ]->insert_neighbor(  1, tBasis[  67 ] );
                 mBasis[   6 ]->insert_neighbor(  2, tBasis[  71 ] );
                 mBasis[   6 ]->insert_neighbor(  3, mBasis[  18 ] );
                 mBasis[   6 ]->insert_neighbor(  4, mBasis[  14 ] );
                 mBasis[   6 ]->insert_neighbor(  5, tBasis[  91 ] );
                 mBasis[   6 ]->insert_neighbor(  6, mBasis[  24 ] );
                 mBasis[   6 ]->insert_neighbor(  7, tBasis[  51 ] );
                 mBasis[   6 ]->insert_neighbor(  8, tBasis[  55 ] );
                 mBasis[   6 ]->insert_neighbor(  9, mBasis[  26 ] );
                 mBasis[   6 ]->insert_neighbor( 10, mBasis[  22 ] );
                 mBasis[   6 ]->insert_neighbor( 11, tBasis[  65 ] );
                 mBasis[   6 ]->insert_neighbor( 12, tBasis[  72 ] );
                 mBasis[   6 ]->insert_neighbor( 13, tBasis[  70 ] );
                 mBasis[   6 ]->insert_neighbor( 14, tBasis[  86 ] );
                 mBasis[   6 ]->insert_neighbor( 15, tBasis[  92 ] );
                 mBasis[   6 ]->insert_neighbor( 16, tBasis[  96 ] );
                 mBasis[   6 ]->insert_neighbor( 17, tBasis[  90 ] );
                 mBasis[   6 ]->insert_neighbor( 18, mBasis[  20 ] );
                 mBasis[   6 ]->insert_neighbor( 19, tBasis[  49 ] );
                 mBasis[   6 ]->insert_neighbor( 20, tBasis[  56 ] );
                 mBasis[   6 ]->insert_neighbor( 21, tBasis[  54 ] );
                 mBasis[   6 ]->insert_neighbor( 22, tBasis[  85 ] );
                 mBasis[   6 ]->insert_neighbor( 23, tBasis[  87 ] );
                 mBasis[   6 ]->insert_neighbor( 24, tBasis[  97 ] );
                 mBasis[   6 ]->insert_neighbor( 25, tBasis[  95 ] );

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
                 mBasis[   7 ]->insert_neighbor(  0, mBasis[  19 ] );
                 mBasis[   7 ]->insert_neighbor(  1, mBasis[  18 ] );
                 mBasis[   7 ]->insert_neighbor(  2, tBasis[  69 ] );
                 mBasis[   7 ]->insert_neighbor(  3, tBasis[  66 ] );
                 mBasis[   7 ]->insert_neighbor(  4, mBasis[  15 ] );
                 mBasis[   7 ]->insert_neighbor(  5, tBasis[  89 ] );
                 mBasis[   7 ]->insert_neighbor(  6, mBasis[  23 ] );
                 mBasis[   7 ]->insert_neighbor(  7, mBasis[  26 ] );
                 mBasis[   7 ]->insert_neighbor(  8, tBasis[  53 ] );
                 mBasis[   7 ]->insert_neighbor(  9, tBasis[  50 ] );
                 mBasis[   7 ]->insert_neighbor( 10, tBasis[  64 ] );
                 mBasis[   7 ]->insert_neighbor( 11, mBasis[  22 ] );
                 mBasis[   7 ]->insert_neighbor( 12, tBasis[  70 ] );
                 mBasis[   7 ]->insert_neighbor( 13, tBasis[  68 ] );
                 mBasis[   7 ]->insert_neighbor( 14, tBasis[  84 ] );
                 mBasis[   7 ]->insert_neighbor( 15, tBasis[  90 ] );
                 mBasis[   7 ]->insert_neighbor( 16, tBasis[  94 ] );
                 mBasis[   7 ]->insert_neighbor( 17, tBasis[  88 ] );
                 mBasis[   7 ]->insert_neighbor( 18, tBasis[  48 ] );
                 mBasis[   7 ]->insert_neighbor( 19, mBasis[  20 ] );
                 mBasis[   7 ]->insert_neighbor( 20, tBasis[  54 ] );
                 mBasis[   7 ]->insert_neighbor( 21, tBasis[  52 ] );
                 mBasis[   7 ]->insert_neighbor( 22, tBasis[  83 ] );
                 mBasis[   7 ]->insert_neighbor( 23, tBasis[  85 ] );
                 mBasis[   7 ]->insert_neighbor( 24, tBasis[  95 ] );
                 mBasis[   7 ]->insert_neighbor( 25, tBasis[  93 ] );

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
                 mBasis[   8 ]->insert_neighbor(  0, tBasis[  27 ] );
                 mBasis[   8 ]->insert_neighbor(  1, mBasis[   1 ] );
                 mBasis[   8 ]->insert_neighbor(  2, mBasis[  21 ] );
                 mBasis[   8 ]->insert_neighbor(  3, mBasis[   0 ] );
                 mBasis[   8 ]->insert_neighbor(  4, tBasis[   7 ] );
                 mBasis[   8 ]->insert_neighbor(  5, mBasis[  25 ] );
                 mBasis[   8 ]->insert_neighbor(  6, tBasis[   2 ] );
                 mBasis[   8 ]->insert_neighbor(  7, tBasis[   8 ] );
                 mBasis[   8 ]->insert_neighbor(  8, tBasis[  12 ] );
                 mBasis[   8 ]->insert_neighbor(  9, tBasis[   6 ] );
                 mBasis[   8 ]->insert_neighbor( 10, tBasis[  26 ] );
                 mBasis[   8 ]->insert_neighbor( 11, tBasis[  28 ] );
                 mBasis[   8 ]->insert_neighbor( 12, mBasis[   9 ] );
                 mBasis[   8 ]->insert_neighbor( 13, mBasis[  11 ] );
                 mBasis[   8 ]->insert_neighbor( 14, tBasis[  43 ] );
                 mBasis[   8 ]->insert_neighbor( 15, mBasis[  13 ] );
                 mBasis[   8 ]->insert_neighbor( 16, mBasis[  20 ] );
                 mBasis[   8 ]->insert_neighbor( 17, mBasis[  12 ] );
                 mBasis[   8 ]->insert_neighbor( 18, tBasis[   1 ] );
                 mBasis[   8 ]->insert_neighbor( 19, tBasis[   3 ] );
                 mBasis[   8 ]->insert_neighbor( 20, tBasis[  13 ] );
                 mBasis[   8 ]->insert_neighbor( 21, tBasis[  11 ] );
                 mBasis[   8 ]->insert_neighbor( 22, tBasis[  42 ] );
                 mBasis[   8 ]->insert_neighbor( 23, tBasis[  44 ] );
                 mBasis[   8 ]->insert_neighbor( 24, mBasis[  24 ] );
                 mBasis[   8 ]->insert_neighbor( 25, mBasis[  23 ] );

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
                 mBasis[   9 ]->insert_neighbor(  0, mBasis[   1 ] );
                 mBasis[   9 ]->insert_neighbor(  1, tBasis[  33 ] );
                 mBasis[   9 ]->insert_neighbor(  2, mBasis[   2 ] );
                 mBasis[   9 ]->insert_neighbor(  3, mBasis[  21 ] );
                 mBasis[   9 ]->insert_neighbor(  4, tBasis[  13 ] );
                 mBasis[   9 ]->insert_neighbor(  5, mBasis[  24 ] );
                 mBasis[   9 ]->insert_neighbor(  6, tBasis[   8 ] );
                 mBasis[   9 ]->insert_neighbor(  7, tBasis[  14 ] );
                 mBasis[   9 ]->insert_neighbor(  8, tBasis[  18 ] );
                 mBasis[   9 ]->insert_neighbor(  9, tBasis[  12 ] );
                 mBasis[   9 ]->insert_neighbor( 10, mBasis[   8 ] );
                 mBasis[   9 ]->insert_neighbor( 11, tBasis[  31 ] );
                 mBasis[   9 ]->insert_neighbor( 12, tBasis[  35 ] );
                 mBasis[   9 ]->insert_neighbor( 13, mBasis[  10 ] );
                 mBasis[   9 ]->insert_neighbor( 14, mBasis[  13 ] );
                 mBasis[   9 ]->insert_neighbor( 15, tBasis[  49 ] );
                 mBasis[   9 ]->insert_neighbor( 16, mBasis[  14 ] );
                 mBasis[   9 ]->insert_neighbor( 17, mBasis[  20 ] );
                 mBasis[   9 ]->insert_neighbor( 18, tBasis[   7 ] );
                 mBasis[   9 ]->insert_neighbor( 19, tBasis[   9 ] );
                 mBasis[   9 ]->insert_neighbor( 20, tBasis[  19 ] );
                 mBasis[   9 ]->insert_neighbor( 21, tBasis[  17 ] );
                 mBasis[   9 ]->insert_neighbor( 22, mBasis[  25 ] );
                 mBasis[   9 ]->insert_neighbor( 23, tBasis[  47 ] );
                 mBasis[   9 ]->insert_neighbor( 24, tBasis[  51 ] );
                 mBasis[   9 ]->insert_neighbor( 25, mBasis[  26 ] );

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
                 mBasis[  10 ]->insert_neighbor(  0, mBasis[  21 ] );
                 mBasis[  10 ]->insert_neighbor(  1, mBasis[   2 ] );
                 mBasis[  10 ]->insert_neighbor(  2, tBasis[  38 ] );
                 mBasis[  10 ]->insert_neighbor(  3, mBasis[   3 ] );
                 mBasis[  10 ]->insert_neighbor(  4, tBasis[  17 ] );
                 mBasis[  10 ]->insert_neighbor(  5, mBasis[  26 ] );
                 mBasis[  10 ]->insert_neighbor(  6, tBasis[  12 ] );
                 mBasis[  10 ]->insert_neighbor(  7, tBasis[  18 ] );
                 mBasis[  10 ]->insert_neighbor(  8, tBasis[  22 ] );
                 mBasis[  10 ]->insert_neighbor(  9, tBasis[  16 ] );
                 mBasis[  10 ]->insert_neighbor( 10, mBasis[  11 ] );
                 mBasis[  10 ]->insert_neighbor( 11, mBasis[   9 ] );
                 mBasis[  10 ]->insert_neighbor( 12, tBasis[  39 ] );
                 mBasis[  10 ]->insert_neighbor( 13, tBasis[  37 ] );
                 mBasis[  10 ]->insert_neighbor( 14, mBasis[  20 ] );
                 mBasis[  10 ]->insert_neighbor( 15, mBasis[  14 ] );
                 mBasis[  10 ]->insert_neighbor( 16, tBasis[  54 ] );
                 mBasis[  10 ]->insert_neighbor( 17, mBasis[  15 ] );
                 mBasis[  10 ]->insert_neighbor( 18, tBasis[  11 ] );
                 mBasis[  10 ]->insert_neighbor( 19, tBasis[  13 ] );
                 mBasis[  10 ]->insert_neighbor( 20, tBasis[  23 ] );
                 mBasis[  10 ]->insert_neighbor( 21, tBasis[  21 ] );
                 mBasis[  10 ]->insert_neighbor( 22, mBasis[  23 ] );
                 mBasis[  10 ]->insert_neighbor( 23, mBasis[  24 ] );
                 mBasis[  10 ]->insert_neighbor( 24, tBasis[  55 ] );
                 mBasis[  10 ]->insert_neighbor( 25, tBasis[  53 ] );

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
                 mBasis[  11 ]->insert_neighbor(  1, mBasis[  21 ] );
                 mBasis[  11 ]->insert_neighbor(  2, mBasis[   3 ] );
                 mBasis[  11 ]->insert_neighbor(  3, tBasis[  32 ] );
                 mBasis[  11 ]->insert_neighbor(  4, tBasis[  11 ] );
                 mBasis[  11 ]->insert_neighbor(  5, mBasis[  23 ] );
                 mBasis[  11 ]->insert_neighbor(  6, tBasis[   6 ] );
                 mBasis[  11 ]->insert_neighbor(  7, tBasis[  12 ] );
                 mBasis[  11 ]->insert_neighbor(  8, tBasis[  16 ] );
                 mBasis[  11 ]->insert_neighbor(  9, tBasis[  10 ] );
                 mBasis[  11 ]->insert_neighbor( 10, tBasis[  30 ] );
                 mBasis[  11 ]->insert_neighbor( 11, mBasis[   8 ] );
                 mBasis[  11 ]->insert_neighbor( 12, mBasis[  10 ] );
                 mBasis[  11 ]->insert_neighbor( 13, tBasis[  34 ] );
                 mBasis[  11 ]->insert_neighbor( 14, mBasis[  12 ] );
                 mBasis[  11 ]->insert_neighbor( 15, mBasis[  20 ] );
                 mBasis[  11 ]->insert_neighbor( 16, mBasis[  15 ] );
                 mBasis[  11 ]->insert_neighbor( 17, tBasis[  48 ] );
                 mBasis[  11 ]->insert_neighbor( 18, tBasis[   5 ] );
                 mBasis[  11 ]->insert_neighbor( 19, tBasis[   7 ] );
                 mBasis[  11 ]->insert_neighbor( 20, tBasis[  17 ] );
                 mBasis[  11 ]->insert_neighbor( 21, tBasis[  15 ] );
                 mBasis[  11 ]->insert_neighbor( 22, tBasis[  46 ] );
                 mBasis[  11 ]->insert_neighbor( 23, mBasis[  25 ] );
                 mBasis[  11 ]->insert_neighbor( 24, mBasis[  26 ] );
                 mBasis[  11 ]->insert_neighbor( 25, tBasis[  50 ] );

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
                 mBasis[  12 ]->insert_neighbor(  0, tBasis[  42 ] );
                 mBasis[  12 ]->insert_neighbor(  1, mBasis[  25 ] );
                 mBasis[  12 ]->insert_neighbor(  2, mBasis[  23 ] );
                 mBasis[  12 ]->insert_neighbor(  3, tBasis[  46 ] );
                 mBasis[  12 ]->insert_neighbor(  4, mBasis[   0 ] );
                 mBasis[  12 ]->insert_neighbor(  5, mBasis[   4 ] );
                 mBasis[  12 ]->insert_neighbor(  6, tBasis[  26 ] );
                 mBasis[  12 ]->insert_neighbor(  7, mBasis[   8 ] );
                 mBasis[  12 ]->insert_neighbor(  8, mBasis[  11 ] );
                 mBasis[  12 ]->insert_neighbor(  9, tBasis[  30 ] );
                 mBasis[  12 ]->insert_neighbor( 10, tBasis[  41 ] );
                 mBasis[  12 ]->insert_neighbor( 11, tBasis[  43 ] );
                 mBasis[  12 ]->insert_neighbor( 12, mBasis[  20 ] );
                 mBasis[  12 ]->insert_neighbor( 13, tBasis[  48 ] );
                 mBasis[  12 ]->insert_neighbor( 14, tBasis[  58 ] );
                 mBasis[  12 ]->insert_neighbor( 15, mBasis[  16 ] );
                 mBasis[  12 ]->insert_neighbor( 16, mBasis[  19 ] );
                 mBasis[  12 ]->insert_neighbor( 17, tBasis[  62 ] );
                 mBasis[  12 ]->insert_neighbor( 18, tBasis[  25 ] );
                 mBasis[  12 ]->insert_neighbor( 19, tBasis[  27 ] );
                 mBasis[  12 ]->insert_neighbor( 20, mBasis[  21 ] );
                 mBasis[  12 ]->insert_neighbor( 21, tBasis[  32 ] );
                 mBasis[  12 ]->insert_neighbor( 22, tBasis[  57 ] );
                 mBasis[  12 ]->insert_neighbor( 23, tBasis[  59 ] );
                 mBasis[  12 ]->insert_neighbor( 24, mBasis[  22 ] );
                 mBasis[  12 ]->insert_neighbor( 25, tBasis[  64 ] );

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
                 mBasis[  13 ]->insert_neighbor(  0, tBasis[  44 ] );
                 mBasis[  13 ]->insert_neighbor(  1, tBasis[  47 ] );
                 mBasis[  13 ]->insert_neighbor(  2, mBasis[  24 ] );
                 mBasis[  13 ]->insert_neighbor(  3, mBasis[  25 ] );
                 mBasis[  13 ]->insert_neighbor(  4, mBasis[   1 ] );
                 mBasis[  13 ]->insert_neighbor(  5, mBasis[   5 ] );
                 mBasis[  13 ]->insert_neighbor(  6, tBasis[  28 ] );
                 mBasis[  13 ]->insert_neighbor(  7, tBasis[  31 ] );
                 mBasis[  13 ]->insert_neighbor(  8, mBasis[   9 ] );
                 mBasis[  13 ]->insert_neighbor(  9, mBasis[   8 ] );
                 mBasis[  13 ]->insert_neighbor( 10, tBasis[  43 ] );
                 mBasis[  13 ]->insert_neighbor( 11, tBasis[  45 ] );
                 mBasis[  13 ]->insert_neighbor( 12, tBasis[  49 ] );
                 mBasis[  13 ]->insert_neighbor( 13, mBasis[  20 ] );
                 mBasis[  13 ]->insert_neighbor( 14, tBasis[  60 ] );
                 mBasis[  13 ]->insert_neighbor( 15, tBasis[  63 ] );
                 mBasis[  13 ]->insert_neighbor( 16, mBasis[  17 ] );
                 mBasis[  13 ]->insert_neighbor( 17, mBasis[  16 ] );
                 mBasis[  13 ]->insert_neighbor( 18, tBasis[  27 ] );
                 mBasis[  13 ]->insert_neighbor( 19, tBasis[  29 ] );
                 mBasis[  13 ]->insert_neighbor( 20, tBasis[  33 ] );
                 mBasis[  13 ]->insert_neighbor( 21, mBasis[  21 ] );
                 mBasis[  13 ]->insert_neighbor( 22, tBasis[  59 ] );
                 mBasis[  13 ]->insert_neighbor( 23, tBasis[  61 ] );
                 mBasis[  13 ]->insert_neighbor( 24, tBasis[  65 ] );
                 mBasis[  13 ]->insert_neighbor( 25, mBasis[  22 ] );

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
                 mBasis[  14 ]->insert_neighbor(  0, mBasis[  24 ] );
                 mBasis[  14 ]->insert_neighbor(  1, tBasis[  51 ] );
                 mBasis[  14 ]->insert_neighbor(  2, tBasis[  55 ] );
                 mBasis[  14 ]->insert_neighbor(  3, mBasis[  26 ] );
                 mBasis[  14 ]->insert_neighbor(  4, mBasis[   2 ] );
                 mBasis[  14 ]->insert_neighbor(  5, mBasis[   6 ] );
                 mBasis[  14 ]->insert_neighbor(  6, mBasis[   9 ] );
                 mBasis[  14 ]->insert_neighbor(  7, tBasis[  35 ] );
                 mBasis[  14 ]->insert_neighbor(  8, tBasis[  39 ] );
                 mBasis[  14 ]->insert_neighbor(  9, mBasis[  10 ] );
                 mBasis[  14 ]->insert_neighbor( 10, mBasis[  20 ] );
                 mBasis[  14 ]->insert_neighbor( 11, tBasis[  49 ] );
                 mBasis[  14 ]->insert_neighbor( 12, tBasis[  56 ] );
                 mBasis[  14 ]->insert_neighbor( 13, tBasis[  54 ] );
                 mBasis[  14 ]->insert_neighbor( 14, mBasis[  17 ] );
                 mBasis[  14 ]->insert_neighbor( 15, tBasis[  67 ] );
                 mBasis[  14 ]->insert_neighbor( 16, tBasis[  71 ] );
                 mBasis[  14 ]->insert_neighbor( 17, mBasis[  18 ] );
                 mBasis[  14 ]->insert_neighbor( 18, mBasis[  21 ] );
                 mBasis[  14 ]->insert_neighbor( 19, tBasis[  33 ] );
                 mBasis[  14 ]->insert_neighbor( 20, tBasis[  40 ] );
                 mBasis[  14 ]->insert_neighbor( 21, tBasis[  38 ] );
                 mBasis[  14 ]->insert_neighbor( 22, mBasis[  22 ] );
                 mBasis[  14 ]->insert_neighbor( 23, tBasis[  65 ] );
                 mBasis[  14 ]->insert_neighbor( 24, tBasis[  72 ] );
                 mBasis[  14 ]->insert_neighbor( 25, tBasis[  70 ] );

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
                 mBasis[  15 ]->insert_neighbor(  0, mBasis[  23 ] );
                 mBasis[  15 ]->insert_neighbor(  1, mBasis[  26 ] );
                 mBasis[  15 ]->insert_neighbor(  2, tBasis[  53 ] );
                 mBasis[  15 ]->insert_neighbor(  3, tBasis[  50 ] );
                 mBasis[  15 ]->insert_neighbor(  4, mBasis[   3 ] );
                 mBasis[  15 ]->insert_neighbor(  5, mBasis[   7 ] );
                 mBasis[  15 ]->insert_neighbor(  6, mBasis[  11 ] );
                 mBasis[  15 ]->insert_neighbor(  7, mBasis[  10 ] );
                 mBasis[  15 ]->insert_neighbor(  8, tBasis[  37 ] );
                 mBasis[  15 ]->insert_neighbor(  9, tBasis[  34 ] );
                 mBasis[  15 ]->insert_neighbor( 10, tBasis[  48 ] );
                 mBasis[  15 ]->insert_neighbor( 11, mBasis[  20 ] );
                 mBasis[  15 ]->insert_neighbor( 12, tBasis[  54 ] );
                 mBasis[  15 ]->insert_neighbor( 13, tBasis[  52 ] );
                 mBasis[  15 ]->insert_neighbor( 14, mBasis[  19 ] );
                 mBasis[  15 ]->insert_neighbor( 15, mBasis[  18 ] );
                 mBasis[  15 ]->insert_neighbor( 16, tBasis[  69 ] );
                 mBasis[  15 ]->insert_neighbor( 17, tBasis[  66 ] );
                 mBasis[  15 ]->insert_neighbor( 18, tBasis[  32 ] );
                 mBasis[  15 ]->insert_neighbor( 19, mBasis[  21 ] );
                 mBasis[  15 ]->insert_neighbor( 20, tBasis[  38 ] );
                 mBasis[  15 ]->insert_neighbor( 21, tBasis[  36 ] );
                 mBasis[  15 ]->insert_neighbor( 22, tBasis[  64 ] );
                 mBasis[  15 ]->insert_neighbor( 23, mBasis[  22 ] );
                 mBasis[  15 ]->insert_neighbor( 24, tBasis[  70 ] );
                 mBasis[  15 ]->insert_neighbor( 25, tBasis[  68 ] );

                 // flag this basis
                 mBasis[  15 ]->flag();
             }

         }

         // test if basis 16 exists
         if ( mBasis[  16 ] != nullptr )
         {
             // test if basis 16 has been processed
             if ( ! mBasis[  16 ]->is_flagged() )
             {
                 // link neighbors of basis 16
                 mBasis[  16 ]->insert_neighbor(  0, tBasis[  59 ] );
                 mBasis[  16 ]->insert_neighbor(  1, mBasis[   5 ] );
                 mBasis[  16 ]->insert_neighbor(  2, mBasis[  22 ] );
                 mBasis[  16 ]->insert_neighbor(  3, mBasis[   4 ] );
                 mBasis[  16 ]->insert_neighbor(  4, mBasis[  25 ] );
                 mBasis[  16 ]->insert_neighbor(  5, tBasis[  80 ] );
                 mBasis[  16 ]->insert_neighbor(  6, tBasis[  43 ] );
                 mBasis[  16 ]->insert_neighbor(  7, mBasis[  13 ] );
                 mBasis[  16 ]->insert_neighbor(  8, mBasis[  20 ] );
                 mBasis[  16 ]->insert_neighbor(  9, mBasis[  12 ] );
                 mBasis[  16 ]->insert_neighbor( 10, tBasis[  58 ] );
                 mBasis[  16 ]->insert_neighbor( 11, tBasis[  60 ] );
                 mBasis[  16 ]->insert_neighbor( 12, mBasis[  17 ] );
                 mBasis[  16 ]->insert_neighbor( 13, mBasis[  19 ] );
                 mBasis[  16 ]->insert_neighbor( 14, tBasis[  75 ] );
                 mBasis[  16 ]->insert_neighbor( 15, tBasis[  81 ] );
                 mBasis[  16 ]->insert_neighbor( 16, tBasis[  85 ] );
                 mBasis[  16 ]->insert_neighbor( 17, tBasis[  79 ] );
                 mBasis[  16 ]->insert_neighbor( 18, tBasis[  42 ] );
                 mBasis[  16 ]->insert_neighbor( 19, tBasis[  44 ] );
                 mBasis[  16 ]->insert_neighbor( 20, mBasis[  24 ] );
                 mBasis[  16 ]->insert_neighbor( 21, mBasis[  23 ] );
                 mBasis[  16 ]->insert_neighbor( 22, tBasis[  74 ] );
                 mBasis[  16 ]->insert_neighbor( 23, tBasis[  76 ] );
                 mBasis[  16 ]->insert_neighbor( 24, tBasis[  86 ] );
                 mBasis[  16 ]->insert_neighbor( 25, tBasis[  84 ] );

                 // flag this basis
                 mBasis[  16 ]->flag();
             }

         }

         // test if basis 17 exists
         if ( mBasis[  17 ] != nullptr )
         {
             // test if basis 17 has been processed
             if ( ! mBasis[  17 ]->is_flagged() )
             {
                 // link neighbors of basis 17
                 mBasis[  17 ]->insert_neighbor(  0, mBasis[   5 ] );
                 mBasis[  17 ]->insert_neighbor(  1, tBasis[  65 ] );
                 mBasis[  17 ]->insert_neighbor(  2, mBasis[   6 ] );
                 mBasis[  17 ]->insert_neighbor(  3, mBasis[  22 ] );
                 mBasis[  17 ]->insert_neighbor(  4, mBasis[  24 ] );
                 mBasis[  17 ]->insert_neighbor(  5, tBasis[  86 ] );
                 mBasis[  17 ]->insert_neighbor(  6, mBasis[  13 ] );
                 mBasis[  17 ]->insert_neighbor(  7, tBasis[  49 ] );
                 mBasis[  17 ]->insert_neighbor(  8, mBasis[  14 ] );
                 mBasis[  17 ]->insert_neighbor(  9, mBasis[  20 ] );
                 mBasis[  17 ]->insert_neighbor( 10, mBasis[  16 ] );
                 mBasis[  17 ]->insert_neighbor( 11, tBasis[  63 ] );
                 mBasis[  17 ]->insert_neighbor( 12, tBasis[  67 ] );
                 mBasis[  17 ]->insert_neighbor( 13, mBasis[  18 ] );
                 mBasis[  17 ]->insert_neighbor( 14, tBasis[  81 ] );
                 mBasis[  17 ]->insert_neighbor( 15, tBasis[  87 ] );
                 mBasis[  17 ]->insert_neighbor( 16, tBasis[  91 ] );
                 mBasis[  17 ]->insert_neighbor( 17, tBasis[  85 ] );
                 mBasis[  17 ]->insert_neighbor( 18, mBasis[  25 ] );
                 mBasis[  17 ]->insert_neighbor( 19, tBasis[  47 ] );
                 mBasis[  17 ]->insert_neighbor( 20, tBasis[  51 ] );
                 mBasis[  17 ]->insert_neighbor( 21, mBasis[  26 ] );
                 mBasis[  17 ]->insert_neighbor( 22, tBasis[  80 ] );
                 mBasis[  17 ]->insert_neighbor( 23, tBasis[  82 ] );
                 mBasis[  17 ]->insert_neighbor( 24, tBasis[  92 ] );
                 mBasis[  17 ]->insert_neighbor( 25, tBasis[  90 ] );

                 // flag this basis
                 mBasis[  17 ]->flag();
             }

         }

         // test if basis 18 exists
         if ( mBasis[  18 ] != nullptr )
         {
             // test if basis 18 has been processed
             if ( ! mBasis[  18 ]->is_flagged() )
             {
                 // link neighbors of basis 18
                 mBasis[  18 ]->insert_neighbor(  0, mBasis[  22 ] );
                 mBasis[  18 ]->insert_neighbor(  1, mBasis[   6 ] );
                 mBasis[  18 ]->insert_neighbor(  2, tBasis[  70 ] );
                 mBasis[  18 ]->insert_neighbor(  3, mBasis[   7 ] );
                 mBasis[  18 ]->insert_neighbor(  4, mBasis[  26 ] );
                 mBasis[  18 ]->insert_neighbor(  5, tBasis[  90 ] );
                 mBasis[  18 ]->insert_neighbor(  6, mBasis[  20 ] );
                 mBasis[  18 ]->insert_neighbor(  7, mBasis[  14 ] );
                 mBasis[  18 ]->insert_neighbor(  8, tBasis[  54 ] );
                 mBasis[  18 ]->insert_neighbor(  9, mBasis[  15 ] );
                 mBasis[  18 ]->insert_neighbor( 10, mBasis[  19 ] );
                 mBasis[  18 ]->insert_neighbor( 11, mBasis[  17 ] );
                 mBasis[  18 ]->insert_neighbor( 12, tBasis[  71 ] );
                 mBasis[  18 ]->insert_neighbor( 13, tBasis[  69 ] );
                 mBasis[  18 ]->insert_neighbor( 14, tBasis[  85 ] );
                 mBasis[  18 ]->insert_neighbor( 15, tBasis[  91 ] );
                 mBasis[  18 ]->insert_neighbor( 16, tBasis[  95 ] );
                 mBasis[  18 ]->insert_neighbor( 17, tBasis[  89 ] );
                 mBasis[  18 ]->insert_neighbor( 18, mBasis[  23 ] );
                 mBasis[  18 ]->insert_neighbor( 19, mBasis[  24 ] );
                 mBasis[  18 ]->insert_neighbor( 20, tBasis[  55 ] );
                 mBasis[  18 ]->insert_neighbor( 21, tBasis[  53 ] );
                 mBasis[  18 ]->insert_neighbor( 22, tBasis[  84 ] );
                 mBasis[  18 ]->insert_neighbor( 23, tBasis[  86 ] );
                 mBasis[  18 ]->insert_neighbor( 24, tBasis[  96 ] );
                 mBasis[  18 ]->insert_neighbor( 25, tBasis[  94 ] );

                 // flag this basis
                 mBasis[  18 ]->flag();
             }

         }

         // test if basis 19 exists
         if ( mBasis[  19 ] != nullptr )
         {
             // test if basis 19 has been processed
             if ( ! mBasis[  19 ]->is_flagged() )
             {
                 // link neighbors of basis 19
                 mBasis[  19 ]->insert_neighbor(  0, mBasis[   4 ] );
                 mBasis[  19 ]->insert_neighbor(  1, mBasis[  22 ] );
                 mBasis[  19 ]->insert_neighbor(  2, mBasis[   7 ] );
                 mBasis[  19 ]->insert_neighbor(  3, tBasis[  64 ] );
                 mBasis[  19 ]->insert_neighbor(  4, mBasis[  23 ] );
                 mBasis[  19 ]->insert_neighbor(  5, tBasis[  84 ] );
                 mBasis[  19 ]->insert_neighbor(  6, mBasis[  12 ] );
                 mBasis[  19 ]->insert_neighbor(  7, mBasis[  20 ] );
                 mBasis[  19 ]->insert_neighbor(  8, mBasis[  15 ] );
                 mBasis[  19 ]->insert_neighbor(  9, tBasis[  48 ] );
                 mBasis[  19 ]->insert_neighbor( 10, tBasis[  62 ] );
                 mBasis[  19 ]->insert_neighbor( 11, mBasis[  16 ] );
                 mBasis[  19 ]->insert_neighbor( 12, mBasis[  18 ] );
                 mBasis[  19 ]->insert_neighbor( 13, tBasis[  66 ] );
                 mBasis[  19 ]->insert_neighbor( 14, tBasis[  79 ] );
                 mBasis[  19 ]->insert_neighbor( 15, tBasis[  85 ] );
                 mBasis[  19 ]->insert_neighbor( 16, tBasis[  89 ] );
                 mBasis[  19 ]->insert_neighbor( 17, tBasis[  83 ] );
                 mBasis[  19 ]->insert_neighbor( 18, tBasis[  46 ] );
                 mBasis[  19 ]->insert_neighbor( 19, mBasis[  25 ] );
                 mBasis[  19 ]->insert_neighbor( 20, mBasis[  26 ] );
                 mBasis[  19 ]->insert_neighbor( 21, tBasis[  50 ] );
                 mBasis[  19 ]->insert_neighbor( 22, tBasis[  78 ] );
                 mBasis[  19 ]->insert_neighbor( 23, tBasis[  80 ] );
                 mBasis[  19 ]->insert_neighbor( 24, tBasis[  90 ] );
                 mBasis[  19 ]->insert_neighbor( 25, tBasis[  88 ] );

                 // flag this basis
                 mBasis[  19 ]->flag();
             }

         }

         // test if basis 20 exists
         if ( mBasis[  20 ] != nullptr )
         {
             // test if basis 20 has been processed
             if ( ! mBasis[  20 ]->is_flagged() )
             {
                 // link neighbors of basis 20
                 mBasis[  20 ]->insert_neighbor(  0, mBasis[  25 ] );
                 mBasis[  20 ]->insert_neighbor(  1, mBasis[  24 ] );
                 mBasis[  20 ]->insert_neighbor(  2, mBasis[  26 ] );
                 mBasis[  20 ]->insert_neighbor(  3, mBasis[  23 ] );
                 mBasis[  20 ]->insert_neighbor(  4, mBasis[  21 ] );
                 mBasis[  20 ]->insert_neighbor(  5, mBasis[  22 ] );
                 mBasis[  20 ]->insert_neighbor(  6, mBasis[   8 ] );
                 mBasis[  20 ]->insert_neighbor(  7, mBasis[   9 ] );
                 mBasis[  20 ]->insert_neighbor(  8, mBasis[  10 ] );
                 mBasis[  20 ]->insert_neighbor(  9, mBasis[  11 ] );
                 mBasis[  20 ]->insert_neighbor( 10, mBasis[  12 ] );
                 mBasis[  20 ]->insert_neighbor( 11, mBasis[  13 ] );
                 mBasis[  20 ]->insert_neighbor( 12, mBasis[  14 ] );
                 mBasis[  20 ]->insert_neighbor( 13, mBasis[  15 ] );
                 mBasis[  20 ]->insert_neighbor( 14, mBasis[  16 ] );
                 mBasis[  20 ]->insert_neighbor( 15, mBasis[  17 ] );
                 mBasis[  20 ]->insert_neighbor( 16, mBasis[  18 ] );
                 mBasis[  20 ]->insert_neighbor( 17, mBasis[  19 ] );
                 mBasis[  20 ]->insert_neighbor( 18, mBasis[   0 ] );
                 mBasis[  20 ]->insert_neighbor( 19, mBasis[   1 ] );
                 mBasis[  20 ]->insert_neighbor( 20, mBasis[   2 ] );
                 mBasis[  20 ]->insert_neighbor( 21, mBasis[   3 ] );
                 mBasis[  20 ]->insert_neighbor( 22, mBasis[   4 ] );
                 mBasis[  20 ]->insert_neighbor( 23, mBasis[   5 ] );
                 mBasis[  20 ]->insert_neighbor( 24, mBasis[   6 ] );
                 mBasis[  20 ]->insert_neighbor( 25, mBasis[   7 ] );

                 // flag this basis
                 mBasis[  20 ]->flag();
             }

         }

         // test if basis 21 exists
         if ( mBasis[  21 ] != nullptr )
         {
             // test if basis 21 has been processed
             if ( ! mBasis[  21 ]->is_flagged() )
             {
                 // link neighbors of basis 21
                 mBasis[  21 ]->insert_neighbor(  0, mBasis[   8 ] );
                 mBasis[  21 ]->insert_neighbor(  1, mBasis[   9 ] );
                 mBasis[  21 ]->insert_neighbor(  2, mBasis[  10 ] );
                 mBasis[  21 ]->insert_neighbor(  3, mBasis[  11 ] );
                 mBasis[  21 ]->insert_neighbor(  4, tBasis[  12 ] );
                 mBasis[  21 ]->insert_neighbor(  5, mBasis[  20 ] );
                 mBasis[  21 ]->insert_neighbor(  6, tBasis[   7 ] );
                 mBasis[  21 ]->insert_neighbor(  7, tBasis[  13 ] );
                 mBasis[  21 ]->insert_neighbor(  8, tBasis[  17 ] );
                 mBasis[  21 ]->insert_neighbor(  9, tBasis[  11 ] );
                 mBasis[  21 ]->insert_neighbor( 10, mBasis[   0 ] );
                 mBasis[  21 ]->insert_neighbor( 11, mBasis[   1 ] );
                 mBasis[  21 ]->insert_neighbor( 12, mBasis[   2 ] );
                 mBasis[  21 ]->insert_neighbor( 13, mBasis[   3 ] );
                 mBasis[  21 ]->insert_neighbor( 14, mBasis[  25 ] );
                 mBasis[  21 ]->insert_neighbor( 15, mBasis[  24 ] );
                 mBasis[  21 ]->insert_neighbor( 16, mBasis[  26 ] );
                 mBasis[  21 ]->insert_neighbor( 17, mBasis[  23 ] );
                 mBasis[  21 ]->insert_neighbor( 18, tBasis[   6 ] );
                 mBasis[  21 ]->insert_neighbor( 19, tBasis[   8 ] );
                 mBasis[  21 ]->insert_neighbor( 20, tBasis[  18 ] );
                 mBasis[  21 ]->insert_neighbor( 21, tBasis[  16 ] );
                 mBasis[  21 ]->insert_neighbor( 22, mBasis[  12 ] );
                 mBasis[  21 ]->insert_neighbor( 23, mBasis[  13 ] );
                 mBasis[  21 ]->insert_neighbor( 24, mBasis[  14 ] );
                 mBasis[  21 ]->insert_neighbor( 25, mBasis[  15 ] );

                 // flag this basis
                 mBasis[  21 ]->flag();
             }

         }

         // test if basis 22 exists
         if ( mBasis[  22 ] != nullptr )
         {
             // test if basis 22 has been processed
             if ( ! mBasis[  22 ]->is_flagged() )
             {
                 // link neighbors of basis 22
                 mBasis[  22 ]->insert_neighbor(  0, mBasis[  16 ] );
                 mBasis[  22 ]->insert_neighbor(  1, mBasis[  17 ] );
                 mBasis[  22 ]->insert_neighbor(  2, mBasis[  18 ] );
                 mBasis[  22 ]->insert_neighbor(  3, mBasis[  19 ] );
                 mBasis[  22 ]->insert_neighbor(  4, mBasis[  20 ] );
                 mBasis[  22 ]->insert_neighbor(  5, tBasis[  85 ] );
                 mBasis[  22 ]->insert_neighbor(  6, mBasis[  25 ] );
                 mBasis[  22 ]->insert_neighbor(  7, mBasis[  24 ] );
                 mBasis[  22 ]->insert_neighbor(  8, mBasis[  26 ] );
                 mBasis[  22 ]->insert_neighbor(  9, mBasis[  23 ] );
                 mBasis[  22 ]->insert_neighbor( 10, mBasis[   4 ] );
                 mBasis[  22 ]->insert_neighbor( 11, mBasis[   5 ] );
                 mBasis[  22 ]->insert_neighbor( 12, mBasis[   6 ] );
                 mBasis[  22 ]->insert_neighbor( 13, mBasis[   7 ] );
                 mBasis[  22 ]->insert_neighbor( 14, tBasis[  80 ] );
                 mBasis[  22 ]->insert_neighbor( 15, tBasis[  86 ] );
                 mBasis[  22 ]->insert_neighbor( 16, tBasis[  90 ] );
                 mBasis[  22 ]->insert_neighbor( 17, tBasis[  84 ] );
                 mBasis[  22 ]->insert_neighbor( 18, mBasis[  12 ] );
                 mBasis[  22 ]->insert_neighbor( 19, mBasis[  13 ] );
                 mBasis[  22 ]->insert_neighbor( 20, mBasis[  14 ] );
                 mBasis[  22 ]->insert_neighbor( 21, mBasis[  15 ] );
                 mBasis[  22 ]->insert_neighbor( 22, tBasis[  79 ] );
                 mBasis[  22 ]->insert_neighbor( 23, tBasis[  81 ] );
                 mBasis[  22 ]->insert_neighbor( 24, tBasis[  91 ] );
                 mBasis[  22 ]->insert_neighbor( 25, tBasis[  89 ] );

                 // flag this basis
                 mBasis[  22 ]->flag();
             }

         }

         // test if basis 23 exists
         if ( mBasis[  23 ] != nullptr )
         {
             // test if basis 23 has been processed
             if ( ! mBasis[  23 ]->is_flagged() )
             {
                 // link neighbors of basis 23
                 mBasis[  23 ]->insert_neighbor(  0, mBasis[  12 ] );
                 mBasis[  23 ]->insert_neighbor(  1, mBasis[  20 ] );
                 mBasis[  23 ]->insert_neighbor(  2, mBasis[  15 ] );
                 mBasis[  23 ]->insert_neighbor(  3, tBasis[  48 ] );
                 mBasis[  23 ]->insert_neighbor(  4, mBasis[  11 ] );
                 mBasis[  23 ]->insert_neighbor(  5, mBasis[  19 ] );
                 mBasis[  23 ]->insert_neighbor(  6, mBasis[   0 ] );
                 mBasis[  23 ]->insert_neighbor(  7, mBasis[  21 ] );
                 mBasis[  23 ]->insert_neighbor(  8, mBasis[   3 ] );
                 mBasis[  23 ]->insert_neighbor(  9, tBasis[  32 ] );
                 mBasis[  23 ]->insert_neighbor( 10, tBasis[  46 ] );
                 mBasis[  23 ]->insert_neighbor( 11, mBasis[  25 ] );
                 mBasis[  23 ]->insert_neighbor( 12, mBasis[  26 ] );
                 mBasis[  23 ]->insert_neighbor( 13, tBasis[  50 ] );
                 mBasis[  23 ]->insert_neighbor( 14, mBasis[   4 ] );
                 mBasis[  23 ]->insert_neighbor( 15, mBasis[  22 ] );
                 mBasis[  23 ]->insert_neighbor( 16, mBasis[   7 ] );
                 mBasis[  23 ]->insert_neighbor( 17, tBasis[  64 ] );
                 mBasis[  23 ]->insert_neighbor( 18, tBasis[  30 ] );
                 mBasis[  23 ]->insert_neighbor( 19, mBasis[   8 ] );
                 mBasis[  23 ]->insert_neighbor( 20, mBasis[  10 ] );
                 mBasis[  23 ]->insert_neighbor( 21, tBasis[  34 ] );
                 mBasis[  23 ]->insert_neighbor( 22, tBasis[  62 ] );
                 mBasis[  23 ]->insert_neighbor( 23, mBasis[  16 ] );
                 mBasis[  23 ]->insert_neighbor( 24, mBasis[  18 ] );
                 mBasis[  23 ]->insert_neighbor( 25, tBasis[  66 ] );

                 // flag this basis
                 mBasis[  23 ]->flag();
             }

         }

         // test if basis 24 exists
         if ( mBasis[  24 ] != nullptr )
         {
             // test if basis 24 has been processed
             if ( ! mBasis[  24 ]->is_flagged() )
             {
                 // link neighbors of basis 24
                 mBasis[  24 ]->insert_neighbor(  0, mBasis[  13 ] );
                 mBasis[  24 ]->insert_neighbor(  1, tBasis[  49 ] );
                 mBasis[  24 ]->insert_neighbor(  2, mBasis[  14 ] );
                 mBasis[  24 ]->insert_neighbor(  3, mBasis[  20 ] );
                 mBasis[  24 ]->insert_neighbor(  4, mBasis[   9 ] );
                 mBasis[  24 ]->insert_neighbor(  5, mBasis[  17 ] );
                 mBasis[  24 ]->insert_neighbor(  6, mBasis[   1 ] );
                 mBasis[  24 ]->insert_neighbor(  7, tBasis[  33 ] );
                 mBasis[  24 ]->insert_neighbor(  8, mBasis[   2 ] );
                 mBasis[  24 ]->insert_neighbor(  9, mBasis[  21 ] );
                 mBasis[  24 ]->insert_neighbor( 10, mBasis[  25 ] );
                 mBasis[  24 ]->insert_neighbor( 11, tBasis[  47 ] );
                 mBasis[  24 ]->insert_neighbor( 12, tBasis[  51 ] );
                 mBasis[  24 ]->insert_neighbor( 13, mBasis[  26 ] );
                 mBasis[  24 ]->insert_neighbor( 14, mBasis[   5 ] );
                 mBasis[  24 ]->insert_neighbor( 15, tBasis[  65 ] );
                 mBasis[  24 ]->insert_neighbor( 16, mBasis[   6 ] );
                 mBasis[  24 ]->insert_neighbor( 17, mBasis[  22 ] );
                 mBasis[  24 ]->insert_neighbor( 18, mBasis[   8 ] );
                 mBasis[  24 ]->insert_neighbor( 19, tBasis[  31 ] );
                 mBasis[  24 ]->insert_neighbor( 20, tBasis[  35 ] );
                 mBasis[  24 ]->insert_neighbor( 21, mBasis[  10 ] );
                 mBasis[  24 ]->insert_neighbor( 22, mBasis[  16 ] );
                 mBasis[  24 ]->insert_neighbor( 23, tBasis[  63 ] );
                 mBasis[  24 ]->insert_neighbor( 24, tBasis[  67 ] );
                 mBasis[  24 ]->insert_neighbor( 25, mBasis[  18 ] );

                 // flag this basis
                 mBasis[  24 ]->flag();
             }

         }

         // test if basis 25 exists
         if ( mBasis[  25 ] != nullptr )
         {
             // test if basis 25 has been processed
             if ( ! mBasis[  25 ]->is_flagged() )
             {
                 // link neighbors of basis 25
                 mBasis[  25 ]->insert_neighbor(  0, tBasis[  43 ] );
                 mBasis[  25 ]->insert_neighbor(  1, mBasis[  13 ] );
                 mBasis[  25 ]->insert_neighbor(  2, mBasis[  20 ] );
                 mBasis[  25 ]->insert_neighbor(  3, mBasis[  12 ] );
                 mBasis[  25 ]->insert_neighbor(  4, mBasis[   8 ] );
                 mBasis[  25 ]->insert_neighbor(  5, mBasis[  16 ] );
                 mBasis[  25 ]->insert_neighbor(  6, tBasis[  27 ] );
                 mBasis[  25 ]->insert_neighbor(  7, mBasis[   1 ] );
                 mBasis[  25 ]->insert_neighbor(  8, mBasis[  21 ] );
                 mBasis[  25 ]->insert_neighbor(  9, mBasis[   0 ] );
                 mBasis[  25 ]->insert_neighbor( 10, tBasis[  42 ] );
                 mBasis[  25 ]->insert_neighbor( 11, tBasis[  44 ] );
                 mBasis[  25 ]->insert_neighbor( 12, mBasis[  24 ] );
                 mBasis[  25 ]->insert_neighbor( 13, mBasis[  23 ] );
                 mBasis[  25 ]->insert_neighbor( 14, tBasis[  59 ] );
                 mBasis[  25 ]->insert_neighbor( 15, mBasis[   5 ] );
                 mBasis[  25 ]->insert_neighbor( 16, mBasis[  22 ] );
                 mBasis[  25 ]->insert_neighbor( 17, mBasis[   4 ] );
                 mBasis[  25 ]->insert_neighbor( 18, tBasis[  26 ] );
                 mBasis[  25 ]->insert_neighbor( 19, tBasis[  28 ] );
                 mBasis[  25 ]->insert_neighbor( 20, mBasis[   9 ] );
                 mBasis[  25 ]->insert_neighbor( 21, mBasis[  11 ] );
                 mBasis[  25 ]->insert_neighbor( 22, tBasis[  58 ] );
                 mBasis[  25 ]->insert_neighbor( 23, tBasis[  60 ] );
                 mBasis[  25 ]->insert_neighbor( 24, mBasis[  17 ] );
                 mBasis[  25 ]->insert_neighbor( 25, mBasis[  19 ] );

                 // flag this basis
                 mBasis[  25 ]->flag();
             }

         }

         // test if basis 26 exists
         if ( mBasis[  26 ] != nullptr )
         {
             // test if basis 26 has been processed
             if ( ! mBasis[  26 ]->is_flagged() )
             {
                 // link neighbors of basis 26
                 mBasis[  26 ]->insert_neighbor(  0, mBasis[  20 ] );
                 mBasis[  26 ]->insert_neighbor(  1, mBasis[  14 ] );
                 mBasis[  26 ]->insert_neighbor(  2, tBasis[  54 ] );
                 mBasis[  26 ]->insert_neighbor(  3, mBasis[  15 ] );
                 mBasis[  26 ]->insert_neighbor(  4, mBasis[  10 ] );
                 mBasis[  26 ]->insert_neighbor(  5, mBasis[  18 ] );
                 mBasis[  26 ]->insert_neighbor(  6, mBasis[  21 ] );
                 mBasis[  26 ]->insert_neighbor(  7, mBasis[   2 ] );
                 mBasis[  26 ]->insert_neighbor(  8, tBasis[  38 ] );
                 mBasis[  26 ]->insert_neighbor(  9, mBasis[   3 ] );
                 mBasis[  26 ]->insert_neighbor( 10, mBasis[  23 ] );
                 mBasis[  26 ]->insert_neighbor( 11, mBasis[  24 ] );
                 mBasis[  26 ]->insert_neighbor( 12, tBasis[  55 ] );
                 mBasis[  26 ]->insert_neighbor( 13, tBasis[  53 ] );
                 mBasis[  26 ]->insert_neighbor( 14, mBasis[  22 ] );
                 mBasis[  26 ]->insert_neighbor( 15, mBasis[   6 ] );
                 mBasis[  26 ]->insert_neighbor( 16, tBasis[  70 ] );
                 mBasis[  26 ]->insert_neighbor( 17, mBasis[   7 ] );
                 mBasis[  26 ]->insert_neighbor( 18, mBasis[  11 ] );
                 mBasis[  26 ]->insert_neighbor( 19, mBasis[   9 ] );
                 mBasis[  26 ]->insert_neighbor( 20, tBasis[  39 ] );
                 mBasis[  26 ]->insert_neighbor( 21, tBasis[  37 ] );
                 mBasis[  26 ]->insert_neighbor( 22, mBasis[  19 ] );
                 mBasis[  26 ]->insert_neighbor( 23, mBasis[  17 ] );
                 mBasis[  26 ]->insert_neighbor( 24, tBasis[  71 ] );
                 mBasis[  26 ]->insert_neighbor( 25, tBasis[  69 ] );

                 // flag this basis
                 mBasis[  26 ]->flag();
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
    luint BSpline_Element< 2, 2, 2 >::refine_basis( uint aBasisNumber )
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
                        tBasis->insert_child(  0, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 29 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 31 ) );
                        tBasis->insert_child( 32, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child( 34, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 61 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  3, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 29 ) );
                        tBasis->insert_child( 34, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 53 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 61 ) );
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
                        tBasis->insert_child(  8, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 40, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 53 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 55 ) );
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
                        tBasis->insert_child(  1, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 31 ) );
                        tBasis->insert_child( 32, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child( 40, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 55 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  0, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 53 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 55 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 61 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child( 32, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 34, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 40, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 29 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 31 ) );
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
                        tBasis->insert_child(  0, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child(  2, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 61 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  2, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 53 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 61 ) );
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
                        tBasis->insert_child(  8, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child( 10, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 53 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 55 ) );
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
                        tBasis->insert_child(  0, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child(  8, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 55 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  0, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 31 ) );
                        tBasis->insert_child( 32, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  2, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 29 ) );
                        tBasis->insert_child( 34, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 61 ) );
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
                        tBasis->insert_child( 10, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 53 ) );
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
                        tBasis->insert_child(  8, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 40, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 55 ) );
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
                        tBasis->insert_child( 32, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 34, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 29 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 31 ) );
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
                        tBasis->insert_child( 34, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 29 ) );
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
                        tBasis->insert_child( 40, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 42, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 21 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 23 ) );
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
                        tBasis->insert_child( 32, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 40, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 23 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 31 ) );
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
                        tBasis->insert_child(  0, tNeighbor->get_child( 42 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 43 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 46 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 47 ) );
                        tBasis->insert_child( 16, tNeighbor->get_child( 58 ) );
                        tBasis->insert_child( 17, tNeighbor->get_child( 59 ) );
                        tBasis->insert_child( 20, tNeighbor->get_child( 62 ) );
                        tBasis->insert_child( 21, tNeighbor->get_child( 63 ) );
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
                        tBasis->insert_child(  2, tNeighbor->get_child( 40 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child( 41 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 44 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 45 ) );
                        tBasis->insert_child( 18, tNeighbor->get_child( 56 ) );
                        tBasis->insert_child( 19, tNeighbor->get_child( 57 ) );
                        tBasis->insert_child( 22, tNeighbor->get_child( 60 ) );
                        tBasis->insert_child( 23, tNeighbor->get_child( 61 ) );
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
                        tBasis->insert_child( 10, tNeighbor->get_child( 32 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child( 33 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child( 36 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child( 37 ) );
                        tBasis->insert_child( 26, tNeighbor->get_child( 48 ) );
                        tBasis->insert_child( 27, tNeighbor->get_child( 49 ) );
                        tBasis->insert_child( 30, tNeighbor->get_child( 52 ) );
                        tBasis->insert_child( 31, tNeighbor->get_child( 53 ) );
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
                        tBasis->insert_child(  8, tNeighbor->get_child( 34 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child( 35 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child( 38 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child( 39 ) );
                        tBasis->insert_child( 24, tNeighbor->get_child( 50 ) );
                        tBasis->insert_child( 25, tNeighbor->get_child( 51 ) );
                        tBasis->insert_child( 28, tNeighbor->get_child( 54 ) );
                        tBasis->insert_child( 29, tNeighbor->get_child( 55 ) );
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
                        tBasis->insert_child( 32, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child( 33, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child( 36, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child( 37, tNeighbor->get_child( 15 ) );
                        tBasis->insert_child( 48, tNeighbor->get_child( 26 ) );
                        tBasis->insert_child( 49, tNeighbor->get_child( 27 ) );
                        tBasis->insert_child( 52, tNeighbor->get_child( 30 ) );
                        tBasis->insert_child( 53, tNeighbor->get_child( 31 ) );
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
                        tBasis->insert_child( 34, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child( 35, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child( 38, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child( 39, tNeighbor->get_child( 13 ) );
                        tBasis->insert_child( 50, tNeighbor->get_child( 24 ) );
                        tBasis->insert_child( 51, tNeighbor->get_child( 25 ) );
                        tBasis->insert_child( 54, tNeighbor->get_child( 28 ) );
                        tBasis->insert_child( 55, tNeighbor->get_child( 29 ) );
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
                        tBasis->insert_child( 42, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 43, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 46, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 47, tNeighbor->get_child(  5 ) );
                        tBasis->insert_child( 58, tNeighbor->get_child( 16 ) );
                        tBasis->insert_child( 59, tNeighbor->get_child( 17 ) );
                        tBasis->insert_child( 62, tNeighbor->get_child( 20 ) );
                        tBasis->insert_child( 63, tNeighbor->get_child( 21 ) );
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
                        tBasis->insert_child( 40, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child( 41, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 44, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 45, tNeighbor->get_child(  7 ) );
                        tBasis->insert_child( 56, tNeighbor->get_child( 18 ) );
                        tBasis->insert_child( 57, tNeighbor->get_child( 19 ) );
                        tBasis->insert_child( 60, tNeighbor->get_child( 22 ) );
                        tBasis->insert_child( 61, tNeighbor->get_child( 23 ) );
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
                luint tIMax = tIMin + 4;

                // maximum j-position
                luint tJMax = tJMin + 4;

                // maximum K-position
                luint tKMax = tKMin + 4;

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
                                     new BSpline< 3, 64, 26 >( tIJK, tLevel, gNoProcOwner ) );

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
    luint BSpline_Element< 2, 2, 2 >::refine( moris::Cell< Element* > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;

        // refine basis if they have not been refined already
        for( uint k=0; k<27; ++k )
        {
            tBasisCounter += this->refine_basis( k );
        }

        // initialize temporary basis pattern
        Basis* tBasis[ 64 ] = { nullptr };

        // populate basis pattern
        if ( mBasis[   0 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   0 ]->get_child(  42 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   0 ]->get_child(  43 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   0 ]->get_child(  46 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   0 ]->get_child(  47 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   0 ]->get_child(  58 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   0 ]->get_child(  59 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   0 ]->get_child(  62 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   0 ]->get_child(  63 );
            }
        }

        if ( mBasis[   1 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   1 ]->get_child(  40 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   1 ]->get_child(  41 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   1 ]->get_child(  44 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   1 ]->get_child(  45 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   1 ]->get_child(  56 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   1 ]->get_child(  57 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   1 ]->get_child(  60 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   1 ]->get_child(  61 );
            }
        }

        if ( mBasis[   2 ] != nullptr )
        {
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   2 ]->get_child(  32 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   2 ]->get_child(  33 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   2 ]->get_child(  36 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   2 ]->get_child(  37 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   2 ]->get_child(  48 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[   2 ]->get_child(  49 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[   2 ]->get_child(  52 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[   2 ]->get_child(  53 );
            }
        }

        if ( mBasis[   3 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   3 ]->get_child(  34 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   3 ]->get_child(  35 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   3 ]->get_child(  38 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   3 ]->get_child(  39 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   3 ]->get_child(  50 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[   3 ]->get_child(  51 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[   3 ]->get_child(  54 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[   3 ]->get_child(  55 );
            }
        }

        if ( mBasis[   4 ] != nullptr )
        {
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[   4 ]->get_child(  10 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[   4 ]->get_child(  11 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[   4 ]->get_child(  14 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[   4 ]->get_child(  15 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[   4 ]->get_child(  26 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[   4 ]->get_child(  27 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[   4 ]->get_child(  30 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[   4 ]->get_child(  31 );
            }
        }

        if ( mBasis[   5 ] != nullptr )
        {
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[   5 ]->get_child(   8 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[   5 ]->get_child(   9 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[   5 ]->get_child(  12 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[   5 ]->get_child(  13 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[   5 ]->get_child(  24 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[   5 ]->get_child(  25 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[   5 ]->get_child(  28 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[   5 ]->get_child(  29 );
            }
        }

        if ( mBasis[   6 ] != nullptr )
        {
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[   6 ]->get_child(   0 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[   6 ]->get_child(   1 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[   6 ]->get_child(   4 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[   6 ]->get_child(   5 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[   6 ]->get_child(  16 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[   6 ]->get_child(  17 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[   6 ]->get_child(  20 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[   6 ]->get_child(  21 );
            }
        }

        if ( mBasis[   7 ] != nullptr )
        {
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[   7 ]->get_child(   2 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[   7 ]->get_child(   3 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[   7 ]->get_child(   6 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[   7 ]->get_child(   7 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[   7 ]->get_child(  18 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[   7 ]->get_child(  19 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[   7 ]->get_child(  22 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[   7 ]->get_child(  23 );
            }
        }

        if ( mBasis[   8 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   8 ]->get_child(  40 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   8 ]->get_child(  41 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   8 ]->get_child(  42 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   8 ]->get_child(  43 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   8 ]->get_child(  44 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   8 ]->get_child(  45 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   8 ]->get_child(  46 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   8 ]->get_child(  47 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   8 ]->get_child(  56 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[   8 ]->get_child(  57 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   8 ]->get_child(  58 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   8 ]->get_child(  59 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   8 ]->get_child(  60 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   8 ]->get_child(  61 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   8 ]->get_child(  62 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   8 ]->get_child(  63 );
            }
        }

        if ( mBasis[   9 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   9 ]->get_child(  32 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   9 ]->get_child(  33 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   9 ]->get_child(  36 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   9 ]->get_child(  37 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   9 ]->get_child(  40 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   9 ]->get_child(  41 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   9 ]->get_child(  44 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   9 ]->get_child(  45 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   9 ]->get_child(  48 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   9 ]->get_child(  49 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[   9 ]->get_child(  52 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   9 ]->get_child(  53 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   9 ]->get_child(  56 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[   9 ]->get_child(  57 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[   9 ]->get_child(  60 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[   9 ]->get_child(  61 );
            }
        }

        if ( mBasis[  10 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  10 ]->get_child(  32 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  10 ]->get_child(  33 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  10 ]->get_child(  34 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  10 ]->get_child(  35 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  10 ]->get_child(  36 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  10 ]->get_child(  37 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  10 ]->get_child(  38 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  10 ]->get_child(  39 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  10 ]->get_child(  48 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  10 ]->get_child(  49 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  10 ]->get_child(  50 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  10 ]->get_child(  51 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  10 ]->get_child(  52 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  10 ]->get_child(  53 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  10 ]->get_child(  54 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  10 ]->get_child(  55 );
            }
        }

        if ( mBasis[  11 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  11 ]->get_child(  34 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  11 ]->get_child(  35 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  11 ]->get_child(  38 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  11 ]->get_child(  39 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  11 ]->get_child(  42 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  11 ]->get_child(  43 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  11 ]->get_child(  46 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  11 ]->get_child(  47 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  11 ]->get_child(  50 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  11 ]->get_child(  51 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  11 ]->get_child(  54 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  11 ]->get_child(  55 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  11 ]->get_child(  58 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  11 ]->get_child(  59 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  11 ]->get_child(  62 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  11 ]->get_child(  63 );
            }
        }

        if ( mBasis[  12 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  12 ]->get_child(  10 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  12 ]->get_child(  11 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  12 ]->get_child(  14 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  12 ]->get_child(  15 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  12 ]->get_child(  26 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  12 ]->get_child(  27 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  12 ]->get_child(  30 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  12 ]->get_child(  31 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  12 ]->get_child(  42 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  12 ]->get_child(  43 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  12 ]->get_child(  46 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  12 ]->get_child(  47 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  12 ]->get_child(  58 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  12 ]->get_child(  59 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  12 ]->get_child(  62 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  12 ]->get_child(  63 );
            }
        }

        if ( mBasis[  13 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  13 ]->get_child(   8 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  13 ]->get_child(   9 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  13 ]->get_child(  12 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  13 ]->get_child(  13 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  13 ]->get_child(  24 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  13 ]->get_child(  25 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  13 ]->get_child(  28 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  13 ]->get_child(  29 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  13 ]->get_child(  40 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  13 ]->get_child(  41 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  13 ]->get_child(  44 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  13 ]->get_child(  45 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  13 ]->get_child(  56 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  13 ]->get_child(  57 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  13 ]->get_child(  60 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  13 ]->get_child(  61 );
            }
        }

        if ( mBasis[  14 ] != nullptr )
        {
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  14 ]->get_child(   0 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  14 ]->get_child(   1 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  14 ]->get_child(   4 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  14 ]->get_child(   5 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  14 ]->get_child(  16 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  14 ]->get_child(  17 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  14 ]->get_child(  20 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  14 ]->get_child(  21 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  14 ]->get_child(  32 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  14 ]->get_child(  33 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  14 ]->get_child(  36 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  14 ]->get_child(  37 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  14 ]->get_child(  48 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  14 ]->get_child(  49 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  14 ]->get_child(  52 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  14 ]->get_child(  53 );
            }
        }

        if ( mBasis[  15 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  15 ]->get_child(   2 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  15 ]->get_child(   3 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  15 ]->get_child(   6 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  15 ]->get_child(   7 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  15 ]->get_child(  18 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  15 ]->get_child(  19 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  15 ]->get_child(  22 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  15 ]->get_child(  23 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  15 ]->get_child(  34 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  15 ]->get_child(  35 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  15 ]->get_child(  38 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  15 ]->get_child(  39 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  15 ]->get_child(  50 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  15 ]->get_child(  51 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  15 ]->get_child(  54 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  15 ]->get_child(  55 );
            }
        }

        if ( mBasis[  16 ] != nullptr )
        {
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  16 ]->get_child(   8 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  16 ]->get_child(   9 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  16 ]->get_child(  10 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  16 ]->get_child(  11 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  16 ]->get_child(  12 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  16 ]->get_child(  13 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  16 ]->get_child(  14 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  16 ]->get_child(  15 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  16 ]->get_child(  24 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  16 ]->get_child(  25 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  16 ]->get_child(  26 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  16 ]->get_child(  27 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  16 ]->get_child(  28 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  16 ]->get_child(  29 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  16 ]->get_child(  30 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  16 ]->get_child(  31 );
            }
        }

        if ( mBasis[  17 ] != nullptr )
        {
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  17 ]->get_child(   0 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  17 ]->get_child(   1 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  17 ]->get_child(   4 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  17 ]->get_child(   5 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  17 ]->get_child(   8 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  17 ]->get_child(   9 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  17 ]->get_child(  12 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  17 ]->get_child(  13 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  17 ]->get_child(  16 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  17 ]->get_child(  17 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  17 ]->get_child(  20 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  17 ]->get_child(  21 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  17 ]->get_child(  24 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  17 ]->get_child(  25 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  17 ]->get_child(  28 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  17 ]->get_child(  29 );
            }
        }

        if ( mBasis[  18 ] != nullptr )
        {
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  18 ]->get_child(   0 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  18 ]->get_child(   1 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  18 ]->get_child(   2 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  18 ]->get_child(   3 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  18 ]->get_child(   4 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  18 ]->get_child(   5 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  18 ]->get_child(   6 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  18 ]->get_child(   7 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  18 ]->get_child(  16 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  18 ]->get_child(  17 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  18 ]->get_child(  18 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  18 ]->get_child(  19 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  18 ]->get_child(  20 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  18 ]->get_child(  21 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  18 ]->get_child(  22 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  18 ]->get_child(  23 );
            }
        }

        if ( mBasis[  19 ] != nullptr )
        {
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  19 ]->get_child(   2 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  19 ]->get_child(   3 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  19 ]->get_child(   6 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  19 ]->get_child(   7 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  19 ]->get_child(  10 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  19 ]->get_child(  11 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  19 ]->get_child(  14 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  19 ]->get_child(  15 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  19 ]->get_child(  18 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  19 ]->get_child(  19 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  19 ]->get_child(  22 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  19 ]->get_child(  23 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  19 ]->get_child(  26 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  19 ]->get_child(  27 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  19 ]->get_child(  30 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  19 ]->get_child(  31 );
            }
        }

        if ( mBasis[  20 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  20 ]->get_child(   0 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  20 ]->get_child(   1 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  20 ]->get_child(   2 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  20 ]->get_child(   3 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  20 ]->get_child(   4 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  20 ]->get_child(   5 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  20 ]->get_child(   6 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  20 ]->get_child(   7 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  20 ]->get_child(   8 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  20 ]->get_child(   9 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  20 ]->get_child(  10 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  20 ]->get_child(  11 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  20 ]->get_child(  12 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  20 ]->get_child(  13 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  20 ]->get_child(  14 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  20 ]->get_child(  15 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  20 ]->get_child(  16 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  20 ]->get_child(  17 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  20 ]->get_child(  18 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  20 ]->get_child(  19 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  20 ]->get_child(  20 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  20 ]->get_child(  21 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  20 ]->get_child(  22 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  20 ]->get_child(  23 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  20 ]->get_child(  24 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  20 ]->get_child(  25 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  20 ]->get_child(  26 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  20 ]->get_child(  27 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  20 ]->get_child(  28 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  20 ]->get_child(  29 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  20 ]->get_child(  30 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  20 ]->get_child(  31 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  20 ]->get_child(  32 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  20 ]->get_child(  33 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  20 ]->get_child(  34 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  20 ]->get_child(  35 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  20 ]->get_child(  36 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  20 ]->get_child(  37 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  20 ]->get_child(  38 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  20 ]->get_child(  39 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  20 ]->get_child(  40 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  20 ]->get_child(  41 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  20 ]->get_child(  42 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  20 ]->get_child(  43 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  20 ]->get_child(  44 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  20 ]->get_child(  45 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  20 ]->get_child(  46 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  20 ]->get_child(  47 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  20 ]->get_child(  48 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  20 ]->get_child(  49 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  20 ]->get_child(  50 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  20 ]->get_child(  51 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  20 ]->get_child(  52 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  20 ]->get_child(  53 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  20 ]->get_child(  54 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  20 ]->get_child(  55 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  20 ]->get_child(  56 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  20 ]->get_child(  57 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  20 ]->get_child(  58 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  20 ]->get_child(  59 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  20 ]->get_child(  60 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  20 ]->get_child(  61 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  20 ]->get_child(  62 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  20 ]->get_child(  63 );
            }
        }

        if ( mBasis[  21 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  21 ]->get_child(  32 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  21 ]->get_child(  33 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  21 ]->get_child(  34 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  21 ]->get_child(  35 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  21 ]->get_child(  36 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  21 ]->get_child(  37 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  21 ]->get_child(  38 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  21 ]->get_child(  39 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  21 ]->get_child(  40 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  21 ]->get_child(  41 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  21 ]->get_child(  42 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  21 ]->get_child(  43 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  21 ]->get_child(  44 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  21 ]->get_child(  45 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  21 ]->get_child(  46 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  21 ]->get_child(  47 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  21 ]->get_child(  48 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  21 ]->get_child(  49 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  21 ]->get_child(  50 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  21 ]->get_child(  51 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  21 ]->get_child(  52 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  21 ]->get_child(  53 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  21 ]->get_child(  54 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  21 ]->get_child(  55 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  21 ]->get_child(  56 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  21 ]->get_child(  57 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  21 ]->get_child(  58 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  21 ]->get_child(  59 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  21 ]->get_child(  60 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  21 ]->get_child(  61 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  21 ]->get_child(  62 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  21 ]->get_child(  63 );
            }
        }

        if ( mBasis[  22 ] != nullptr )
        {
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  22 ]->get_child(   0 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  22 ]->get_child(   1 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  22 ]->get_child(   2 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  22 ]->get_child(   3 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  22 ]->get_child(   4 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  22 ]->get_child(   5 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  22 ]->get_child(   6 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  22 ]->get_child(   7 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  22 ]->get_child(   8 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  22 ]->get_child(   9 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  22 ]->get_child(  10 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  22 ]->get_child(  11 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  22 ]->get_child(  12 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  22 ]->get_child(  13 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  22 ]->get_child(  14 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  22 ]->get_child(  15 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  22 ]->get_child(  16 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  22 ]->get_child(  17 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  22 ]->get_child(  18 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  22 ]->get_child(  19 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  22 ]->get_child(  20 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  22 ]->get_child(  21 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  22 ]->get_child(  22 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  22 ]->get_child(  23 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  22 ]->get_child(  24 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  22 ]->get_child(  25 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  22 ]->get_child(  26 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  22 ]->get_child(  27 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  22 ]->get_child(  28 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  22 ]->get_child(  29 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  22 ]->get_child(  30 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  22 ]->get_child(  31 );
            }
        }

        if ( mBasis[  23 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  23 ]->get_child(   2 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  23 ]->get_child(   3 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  23 ]->get_child(   6 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  23 ]->get_child(   7 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  23 ]->get_child(  10 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  23 ]->get_child(  11 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  23 ]->get_child(  14 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  23 ]->get_child(  15 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  23 ]->get_child(  18 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  23 ]->get_child(  19 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  23 ]->get_child(  22 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  23 ]->get_child(  23 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  23 ]->get_child(  26 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  23 ]->get_child(  27 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  23 ]->get_child(  30 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  23 ]->get_child(  31 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  23 ]->get_child(  34 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  23 ]->get_child(  35 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  23 ]->get_child(  38 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  23 ]->get_child(  39 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  23 ]->get_child(  42 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  23 ]->get_child(  43 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  23 ]->get_child(  46 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  23 ]->get_child(  47 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  23 ]->get_child(  50 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  23 ]->get_child(  51 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  23 ]->get_child(  54 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  23 ]->get_child(  55 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  23 ]->get_child(  58 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  23 ]->get_child(  59 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  23 ]->get_child(  62 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  23 ]->get_child(  63 );
            }
        }

        if ( mBasis[  24 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  24 ]->get_child(   0 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  24 ]->get_child(   1 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  24 ]->get_child(   4 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  24 ]->get_child(   5 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  24 ]->get_child(   8 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  24 ]->get_child(   9 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  24 ]->get_child(  12 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  24 ]->get_child(  13 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  24 ]->get_child(  16 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  24 ]->get_child(  17 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  24 ]->get_child(  20 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  24 ]->get_child(  21 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  24 ]->get_child(  24 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  24 ]->get_child(  25 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  24 ]->get_child(  28 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  24 ]->get_child(  29 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  24 ]->get_child(  32 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  24 ]->get_child(  33 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  24 ]->get_child(  36 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  24 ]->get_child(  37 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  24 ]->get_child(  40 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  24 ]->get_child(  41 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  24 ]->get_child(  44 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  24 ]->get_child(  45 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  24 ]->get_child(  48 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  24 ]->get_child(  49 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  24 ]->get_child(  52 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  24 ]->get_child(  53 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  24 ]->get_child(  56 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  24 ]->get_child(  57 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  24 ]->get_child(  60 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  24 ]->get_child(  61 );
            }
        }

        if ( mBasis[  25 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  25 ]->get_child(   8 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  25 ]->get_child(   9 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  25 ]->get_child(  10 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  25 ]->get_child(  11 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  25 ]->get_child(  12 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  25 ]->get_child(  13 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  25 ]->get_child(  14 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  25 ]->get_child(  15 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  25 ]->get_child(  24 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  25 ]->get_child(  25 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  25 ]->get_child(  26 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  25 ]->get_child(  27 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  25 ]->get_child(  28 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  25 ]->get_child(  29 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  25 ]->get_child(  30 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  25 ]->get_child(  31 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  25 ]->get_child(  40 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  25 ]->get_child(  41 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  25 ]->get_child(  42 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  25 ]->get_child(  43 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  25 ]->get_child(  44 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  25 ]->get_child(  45 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  25 ]->get_child(  46 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  25 ]->get_child(  47 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  25 ]->get_child(  56 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  25 ]->get_child(  57 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  25 ]->get_child(  58 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  25 ]->get_child(  59 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  25 ]->get_child(  60 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  25 ]->get_child(  61 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  25 ]->get_child(  62 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  25 ]->get_child(  63 );
            }
        }

        if ( mBasis[  26 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  26 ]->get_child(   0 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  26 ]->get_child(   1 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  26 ]->get_child(   2 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  26 ]->get_child(   3 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  26 ]->get_child(   4 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  26 ]->get_child(   5 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  26 ]->get_child(   6 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  26 ]->get_child(   7 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  26 ]->get_child(  16 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  26 ]->get_child(  17 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  26 ]->get_child(  18 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  26 ]->get_child(  19 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  26 ]->get_child(  20 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  26 ]->get_child(  21 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  26 ]->get_child(  22 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  26 ]->get_child(  23 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  26 ]->get_child(  32 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  26 ]->get_child(  33 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  26 ]->get_child(  34 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  26 ]->get_child(  35 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  26 ]->get_child(  36 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  26 ]->get_child(  37 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  26 ]->get_child(  38 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  26 ]->get_child(  39 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  26 ]->get_child(  48 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  26 ]->get_child(  49 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  26 ]->get_child(  50 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  26 ]->get_child(  51 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  26 ]->get_child(  52 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  26 ]->get_child(  53 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  26 ]->get_child(  54 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  26 ]->get_child(  55 );
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
        tChildren[ 0 ]->insert_basis(  1, tBasis[   2 ] );
        tChildren[ 0 ]->insert_basis(  2, tBasis[  10 ] );
        tChildren[ 0 ]->insert_basis(  3, tBasis[   8 ] );
        tChildren[ 0 ]->insert_basis(  4, tBasis[  32 ] );
        tChildren[ 0 ]->insert_basis(  5, tBasis[  34 ] );
        tChildren[ 0 ]->insert_basis(  6, tBasis[  42 ] );
        tChildren[ 0 ]->insert_basis(  7, tBasis[  40 ] );
        tChildren[ 0 ]->insert_basis(  8, tBasis[   1 ] );
        tChildren[ 0 ]->insert_basis(  9, tBasis[   6 ] );
        tChildren[ 0 ]->insert_basis( 10, tBasis[   9 ] );
        tChildren[ 0 ]->insert_basis( 11, tBasis[   4 ] );
        tChildren[ 0 ]->insert_basis( 12, tBasis[  16 ] );
        tChildren[ 0 ]->insert_basis( 13, tBasis[  18 ] );
        tChildren[ 0 ]->insert_basis( 14, tBasis[  26 ] );
        tChildren[ 0 ]->insert_basis( 15, tBasis[  24 ] );
        tChildren[ 0 ]->insert_basis( 16, tBasis[  33 ] );
        tChildren[ 0 ]->insert_basis( 17, tBasis[  38 ] );
        tChildren[ 0 ]->insert_basis( 18, tBasis[  41 ] );
        tChildren[ 0 ]->insert_basis( 19, tBasis[  36 ] );
        tChildren[ 0 ]->insert_basis( 20, tBasis[  21 ] );
        tChildren[ 0 ]->insert_basis( 21, tBasis[   5 ] );
        tChildren[ 0 ]->insert_basis( 22, tBasis[  37 ] );
        tChildren[ 0 ]->insert_basis( 23, tBasis[  20 ] );
        tChildren[ 0 ]->insert_basis( 24, tBasis[  22 ] );
        tChildren[ 0 ]->insert_basis( 25, tBasis[  17 ] );
        tChildren[ 0 ]->insert_basis( 26, tBasis[  25 ] );

        // assign basis to child 2
        tChildren[ 1 ]->insert_basis(  0, tBasis[   1 ] );
        tChildren[ 1 ]->insert_basis(  1, tBasis[   3 ] );
        tChildren[ 1 ]->insert_basis(  2, tBasis[  11 ] );
        tChildren[ 1 ]->insert_basis(  3, tBasis[   9 ] );
        tChildren[ 1 ]->insert_basis(  4, tBasis[  33 ] );
        tChildren[ 1 ]->insert_basis(  5, tBasis[  35 ] );
        tChildren[ 1 ]->insert_basis(  6, tBasis[  43 ] );
        tChildren[ 1 ]->insert_basis(  7, tBasis[  41 ] );
        tChildren[ 1 ]->insert_basis(  8, tBasis[   2 ] );
        tChildren[ 1 ]->insert_basis(  9, tBasis[   7 ] );
        tChildren[ 1 ]->insert_basis( 10, tBasis[  10 ] );
        tChildren[ 1 ]->insert_basis( 11, tBasis[   5 ] );
        tChildren[ 1 ]->insert_basis( 12, tBasis[  17 ] );
        tChildren[ 1 ]->insert_basis( 13, tBasis[  19 ] );
        tChildren[ 1 ]->insert_basis( 14, tBasis[  27 ] );
        tChildren[ 1 ]->insert_basis( 15, tBasis[  25 ] );
        tChildren[ 1 ]->insert_basis( 16, tBasis[  34 ] );
        tChildren[ 1 ]->insert_basis( 17, tBasis[  39 ] );
        tChildren[ 1 ]->insert_basis( 18, tBasis[  42 ] );
        tChildren[ 1 ]->insert_basis( 19, tBasis[  37 ] );
        tChildren[ 1 ]->insert_basis( 20, tBasis[  22 ] );
        tChildren[ 1 ]->insert_basis( 21, tBasis[   6 ] );
        tChildren[ 1 ]->insert_basis( 22, tBasis[  38 ] );
        tChildren[ 1 ]->insert_basis( 23, tBasis[  21 ] );
        tChildren[ 1 ]->insert_basis( 24, tBasis[  23 ] );
        tChildren[ 1 ]->insert_basis( 25, tBasis[  18 ] );
        tChildren[ 1 ]->insert_basis( 26, tBasis[  26 ] );

        // assign basis to child 3
        tChildren[ 2 ]->insert_basis(  0, tBasis[   4 ] );
        tChildren[ 2 ]->insert_basis(  1, tBasis[   6 ] );
        tChildren[ 2 ]->insert_basis(  2, tBasis[  14 ] );
        tChildren[ 2 ]->insert_basis(  3, tBasis[  12 ] );
        tChildren[ 2 ]->insert_basis(  4, tBasis[  36 ] );
        tChildren[ 2 ]->insert_basis(  5, tBasis[  38 ] );
        tChildren[ 2 ]->insert_basis(  6, tBasis[  46 ] );
        tChildren[ 2 ]->insert_basis(  7, tBasis[  44 ] );
        tChildren[ 2 ]->insert_basis(  8, tBasis[   5 ] );
        tChildren[ 2 ]->insert_basis(  9, tBasis[  10 ] );
        tChildren[ 2 ]->insert_basis( 10, tBasis[  13 ] );
        tChildren[ 2 ]->insert_basis( 11, tBasis[   8 ] );
        tChildren[ 2 ]->insert_basis( 12, tBasis[  20 ] );
        tChildren[ 2 ]->insert_basis( 13, tBasis[  22 ] );
        tChildren[ 2 ]->insert_basis( 14, tBasis[  30 ] );
        tChildren[ 2 ]->insert_basis( 15, tBasis[  28 ] );
        tChildren[ 2 ]->insert_basis( 16, tBasis[  37 ] );
        tChildren[ 2 ]->insert_basis( 17, tBasis[  42 ] );
        tChildren[ 2 ]->insert_basis( 18, tBasis[  45 ] );
        tChildren[ 2 ]->insert_basis( 19, tBasis[  40 ] );
        tChildren[ 2 ]->insert_basis( 20, tBasis[  25 ] );
        tChildren[ 2 ]->insert_basis( 21, tBasis[   9 ] );
        tChildren[ 2 ]->insert_basis( 22, tBasis[  41 ] );
        tChildren[ 2 ]->insert_basis( 23, tBasis[  24 ] );
        tChildren[ 2 ]->insert_basis( 24, tBasis[  26 ] );
        tChildren[ 2 ]->insert_basis( 25, tBasis[  21 ] );
        tChildren[ 2 ]->insert_basis( 26, tBasis[  29 ] );

        // assign basis to child 4
        tChildren[ 3 ]->insert_basis(  0, tBasis[   5 ] );
        tChildren[ 3 ]->insert_basis(  1, tBasis[   7 ] );
        tChildren[ 3 ]->insert_basis(  2, tBasis[  15 ] );
        tChildren[ 3 ]->insert_basis(  3, tBasis[  13 ] );
        tChildren[ 3 ]->insert_basis(  4, tBasis[  37 ] );
        tChildren[ 3 ]->insert_basis(  5, tBasis[  39 ] );
        tChildren[ 3 ]->insert_basis(  6, tBasis[  47 ] );
        tChildren[ 3 ]->insert_basis(  7, tBasis[  45 ] );
        tChildren[ 3 ]->insert_basis(  8, tBasis[   6 ] );
        tChildren[ 3 ]->insert_basis(  9, tBasis[  11 ] );
        tChildren[ 3 ]->insert_basis( 10, tBasis[  14 ] );
        tChildren[ 3 ]->insert_basis( 11, tBasis[   9 ] );
        tChildren[ 3 ]->insert_basis( 12, tBasis[  21 ] );
        tChildren[ 3 ]->insert_basis( 13, tBasis[  23 ] );
        tChildren[ 3 ]->insert_basis( 14, tBasis[  31 ] );
        tChildren[ 3 ]->insert_basis( 15, tBasis[  29 ] );
        tChildren[ 3 ]->insert_basis( 16, tBasis[  38 ] );
        tChildren[ 3 ]->insert_basis( 17, tBasis[  43 ] );
        tChildren[ 3 ]->insert_basis( 18, tBasis[  46 ] );
        tChildren[ 3 ]->insert_basis( 19, tBasis[  41 ] );
        tChildren[ 3 ]->insert_basis( 20, tBasis[  26 ] );
        tChildren[ 3 ]->insert_basis( 21, tBasis[  10 ] );
        tChildren[ 3 ]->insert_basis( 22, tBasis[  42 ] );
        tChildren[ 3 ]->insert_basis( 23, tBasis[  25 ] );
        tChildren[ 3 ]->insert_basis( 24, tBasis[  27 ] );
        tChildren[ 3 ]->insert_basis( 25, tBasis[  22 ] );
        tChildren[ 3 ]->insert_basis( 26, tBasis[  30 ] );

        // assign basis to child 5
        tChildren[ 4 ]->insert_basis(  0, tBasis[  16 ] );
        tChildren[ 4 ]->insert_basis(  1, tBasis[  18 ] );
        tChildren[ 4 ]->insert_basis(  2, tBasis[  26 ] );
        tChildren[ 4 ]->insert_basis(  3, tBasis[  24 ] );
        tChildren[ 4 ]->insert_basis(  4, tBasis[  48 ] );
        tChildren[ 4 ]->insert_basis(  5, tBasis[  50 ] );
        tChildren[ 4 ]->insert_basis(  6, tBasis[  58 ] );
        tChildren[ 4 ]->insert_basis(  7, tBasis[  56 ] );
        tChildren[ 4 ]->insert_basis(  8, tBasis[  17 ] );
        tChildren[ 4 ]->insert_basis(  9, tBasis[  22 ] );
        tChildren[ 4 ]->insert_basis( 10, tBasis[  25 ] );
        tChildren[ 4 ]->insert_basis( 11, tBasis[  20 ] );
        tChildren[ 4 ]->insert_basis( 12, tBasis[  32 ] );
        tChildren[ 4 ]->insert_basis( 13, tBasis[  34 ] );
        tChildren[ 4 ]->insert_basis( 14, tBasis[  42 ] );
        tChildren[ 4 ]->insert_basis( 15, tBasis[  40 ] );
        tChildren[ 4 ]->insert_basis( 16, tBasis[  49 ] );
        tChildren[ 4 ]->insert_basis( 17, tBasis[  54 ] );
        tChildren[ 4 ]->insert_basis( 18, tBasis[  57 ] );
        tChildren[ 4 ]->insert_basis( 19, tBasis[  52 ] );
        tChildren[ 4 ]->insert_basis( 20, tBasis[  37 ] );
        tChildren[ 4 ]->insert_basis( 21, tBasis[  21 ] );
        tChildren[ 4 ]->insert_basis( 22, tBasis[  53 ] );
        tChildren[ 4 ]->insert_basis( 23, tBasis[  36 ] );
        tChildren[ 4 ]->insert_basis( 24, tBasis[  38 ] );
        tChildren[ 4 ]->insert_basis( 25, tBasis[  33 ] );
        tChildren[ 4 ]->insert_basis( 26, tBasis[  41 ] );

        // assign basis to child 6
        tChildren[ 5 ]->insert_basis(  0, tBasis[  17 ] );
        tChildren[ 5 ]->insert_basis(  1, tBasis[  19 ] );
        tChildren[ 5 ]->insert_basis(  2, tBasis[  27 ] );
        tChildren[ 5 ]->insert_basis(  3, tBasis[  25 ] );
        tChildren[ 5 ]->insert_basis(  4, tBasis[  49 ] );
        tChildren[ 5 ]->insert_basis(  5, tBasis[  51 ] );
        tChildren[ 5 ]->insert_basis(  6, tBasis[  59 ] );
        tChildren[ 5 ]->insert_basis(  7, tBasis[  57 ] );
        tChildren[ 5 ]->insert_basis(  8, tBasis[  18 ] );
        tChildren[ 5 ]->insert_basis(  9, tBasis[  23 ] );
        tChildren[ 5 ]->insert_basis( 10, tBasis[  26 ] );
        tChildren[ 5 ]->insert_basis( 11, tBasis[  21 ] );
        tChildren[ 5 ]->insert_basis( 12, tBasis[  33 ] );
        tChildren[ 5 ]->insert_basis( 13, tBasis[  35 ] );
        tChildren[ 5 ]->insert_basis( 14, tBasis[  43 ] );
        tChildren[ 5 ]->insert_basis( 15, tBasis[  41 ] );
        tChildren[ 5 ]->insert_basis( 16, tBasis[  50 ] );
        tChildren[ 5 ]->insert_basis( 17, tBasis[  55 ] );
        tChildren[ 5 ]->insert_basis( 18, tBasis[  58 ] );
        tChildren[ 5 ]->insert_basis( 19, tBasis[  53 ] );
        tChildren[ 5 ]->insert_basis( 20, tBasis[  38 ] );
        tChildren[ 5 ]->insert_basis( 21, tBasis[  22 ] );
        tChildren[ 5 ]->insert_basis( 22, tBasis[  54 ] );
        tChildren[ 5 ]->insert_basis( 23, tBasis[  37 ] );
        tChildren[ 5 ]->insert_basis( 24, tBasis[  39 ] );
        tChildren[ 5 ]->insert_basis( 25, tBasis[  34 ] );
        tChildren[ 5 ]->insert_basis( 26, tBasis[  42 ] );

        // assign basis to child 7
        tChildren[ 6 ]->insert_basis(  0, tBasis[  20 ] );
        tChildren[ 6 ]->insert_basis(  1, tBasis[  22 ] );
        tChildren[ 6 ]->insert_basis(  2, tBasis[  30 ] );
        tChildren[ 6 ]->insert_basis(  3, tBasis[  28 ] );
        tChildren[ 6 ]->insert_basis(  4, tBasis[  52 ] );
        tChildren[ 6 ]->insert_basis(  5, tBasis[  54 ] );
        tChildren[ 6 ]->insert_basis(  6, tBasis[  62 ] );
        tChildren[ 6 ]->insert_basis(  7, tBasis[  60 ] );
        tChildren[ 6 ]->insert_basis(  8, tBasis[  21 ] );
        tChildren[ 6 ]->insert_basis(  9, tBasis[  26 ] );
        tChildren[ 6 ]->insert_basis( 10, tBasis[  29 ] );
        tChildren[ 6 ]->insert_basis( 11, tBasis[  24 ] );
        tChildren[ 6 ]->insert_basis( 12, tBasis[  36 ] );
        tChildren[ 6 ]->insert_basis( 13, tBasis[  38 ] );
        tChildren[ 6 ]->insert_basis( 14, tBasis[  46 ] );
        tChildren[ 6 ]->insert_basis( 15, tBasis[  44 ] );
        tChildren[ 6 ]->insert_basis( 16, tBasis[  53 ] );
        tChildren[ 6 ]->insert_basis( 17, tBasis[  58 ] );
        tChildren[ 6 ]->insert_basis( 18, tBasis[  61 ] );
        tChildren[ 6 ]->insert_basis( 19, tBasis[  56 ] );
        tChildren[ 6 ]->insert_basis( 20, tBasis[  41 ] );
        tChildren[ 6 ]->insert_basis( 21, tBasis[  25 ] );
        tChildren[ 6 ]->insert_basis( 22, tBasis[  57 ] );
        tChildren[ 6 ]->insert_basis( 23, tBasis[  40 ] );
        tChildren[ 6 ]->insert_basis( 24, tBasis[  42 ] );
        tChildren[ 6 ]->insert_basis( 25, tBasis[  37 ] );
        tChildren[ 6 ]->insert_basis( 26, tBasis[  45 ] );

        // assign basis to child 8
        tChildren[ 7 ]->insert_basis(  0, tBasis[  21 ] );
        tChildren[ 7 ]->insert_basis(  1, tBasis[  23 ] );
        tChildren[ 7 ]->insert_basis(  2, tBasis[  31 ] );
        tChildren[ 7 ]->insert_basis(  3, tBasis[  29 ] );
        tChildren[ 7 ]->insert_basis(  4, tBasis[  53 ] );
        tChildren[ 7 ]->insert_basis(  5, tBasis[  55 ] );
        tChildren[ 7 ]->insert_basis(  6, tBasis[  63 ] );
        tChildren[ 7 ]->insert_basis(  7, tBasis[  61 ] );
        tChildren[ 7 ]->insert_basis(  8, tBasis[  22 ] );
        tChildren[ 7 ]->insert_basis(  9, tBasis[  27 ] );
        tChildren[ 7 ]->insert_basis( 10, tBasis[  30 ] );
        tChildren[ 7 ]->insert_basis( 11, tBasis[  25 ] );
        tChildren[ 7 ]->insert_basis( 12, tBasis[  37 ] );
        tChildren[ 7 ]->insert_basis( 13, tBasis[  39 ] );
        tChildren[ 7 ]->insert_basis( 14, tBasis[  47 ] );
        tChildren[ 7 ]->insert_basis( 15, tBasis[  45 ] );
        tChildren[ 7 ]->insert_basis( 16, tBasis[  54 ] );
        tChildren[ 7 ]->insert_basis( 17, tBasis[  59 ] );
        tChildren[ 7 ]->insert_basis( 18, tBasis[  62 ] );
        tChildren[ 7 ]->insert_basis( 19, tBasis[  57 ] );
        tChildren[ 7 ]->insert_basis( 20, tBasis[  42 ] );
        tChildren[ 7 ]->insert_basis( 21, tBasis[  26 ] );
        tChildren[ 7 ]->insert_basis( 22, tBasis[  58 ] );
        tChildren[ 7 ]->insert_basis( 23, tBasis[  41 ] );
        tChildren[ 7 ]->insert_basis( 24, tBasis[  43 ] );
        tChildren[ 7 ]->insert_basis( 25, tBasis[  38 ] );
        tChildren[ 7 ]->insert_basis( 26, tBasis[  46 ] );

        // set basis flag of element
        mChildrenBasisFlag = true;

        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX27_HPP_ */

