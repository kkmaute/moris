/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element_Hex64.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX64_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX64_HPP_

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
    BSpline_Element< 3, 3, 3 >::get_geometry_type() const
    {
        return mtk::Geometry_Type::HEX;
    }

// ----------------------------------------------------------------------------

    /**
     * returns the ijk position of a given basis function
     *
     * @param[in]  aBasisNumber   element local number of basis function
     * @param[out] aIJK           proc local ijk position of this basis function
     *
     * @return void
     *
     */
    template<>
    void
    BSpline_Element< 3, 3, 3 >::get_ijk_of_basis(
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
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  2 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  3 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case(  4 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case(  5 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case(  6 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case(  7 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 3 ;
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
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 10 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 11 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
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
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 14 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 15 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 16 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 17 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 18 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 19 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 20 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 21 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 22 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 23 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 24 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 25 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 26 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 27 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 28 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 29 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 30 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 31 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 32 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 33 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 34 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 35 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 0 ;
                break;
            }
            case( 36 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 37 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 38 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 39 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 40 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 41 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 42 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 43 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 44 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 45 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 46 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 47 ) :
            {
                aIJK[ 0 ] = 3 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 48 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 49 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 50 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 51 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 3 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 52 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 53 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 54 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 55 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 3 ;
                break;
            }
            case( 56 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 57 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 58 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 59 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 1 ;
                break;
            }
            case( 60 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 61 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 62 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
                break;
            }
            case( 63 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                aIJK[ 2 ] = 2 ;
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
    BSpline_Element< 3, 3, 3 >::link_basis_with_neighbors(
          moris::Cell< Element* > & aAllElementsOnProc )
    {
         // initialize frame of basis around basis from this element
         Basis_Function* tBasis[ 152 ] = { nullptr };

         // get pointer to neighbor  0
         Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

         // test if neighbor  0 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  37 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  38 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  39 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  40 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  57 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  60 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  77 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  80 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  97 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  98 ] = tNeighbor->get_basis_function(  24 );
             tBasis[  99 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 100 ] = tNeighbor->get_basis_function(   5 );
         }

         // get pointer to neighbor  1
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

         // test if neighbor  1 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  43 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  45 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  47 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  49 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  63 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  69 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  83 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  89 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 103 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 109 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  2
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

         // test if neighbor  2 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  51 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  52 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  53 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  54 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  71 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  74 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  91 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  94 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 111 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 114 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  3
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

         // test if neighbor  3 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  42 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  44 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  46 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  48 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  62 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  68 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  82 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  88 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 102 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 108 ] = tNeighbor->get_basis_function(   7 );
         }

         // get pointer to neighbor  4
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

         // test if neighbor  4 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   7 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   8 ] = tNeighbor->get_basis_function(   8 );
             tBasis[   9 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  10 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  13 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  16 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  19 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  22 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  25 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  26 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  27 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  28 ] = tNeighbor->get_basis_function(   2 );
         }

         // get pointer to neighbor  5
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

         // test if neighbor  5 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[ 123 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 126 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 141 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 144 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  6
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

         // test if neighbor  6 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   1 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   2 ] = tNeighbor->get_basis_function(   8 );
             tBasis[   3 ] = tNeighbor->get_basis_function(   9 );
             tBasis[   4 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   7 ] = tNeighbor->get_basis_function(  10 );
             tBasis[   8 ] = tNeighbor->get_basis_function(  32 );
             tBasis[   9 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  10 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  13 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  16 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  19 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  22 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  37 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  38 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  39 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  40 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  57 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  60 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  77 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  24 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  25 );
             tBasis[  80 ] = tNeighbor->get_basis_function(   5 );
         }

         // get pointer to neighbor  7
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

         // test if neighbor  7 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   8 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   9 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  10 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  11 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  16 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  17 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  22 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  23 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  26 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  27 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  28 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  29 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  43 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  45 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  47 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  49 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  63 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  69 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  83 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  28 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  29 );
             tBasis[  89 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  8
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 8 );

         // test if neighbor  8 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  13 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  14 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  15 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  16 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  19 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  22 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  25 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  26 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  27 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  28 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  31 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  32 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  33 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  34 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  51 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  52 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  53 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  54 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  71 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  74 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  91 ] = tNeighbor->get_basis_function(   7 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  31 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  30 );
             tBasis[  94 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  9
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 9 );

         // test if neighbor  9 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   6 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   7 ] = tNeighbor->get_basis_function(   8 );
             tBasis[   8 ] = tNeighbor->get_basis_function(   9 );
             tBasis[   9 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  12 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  13 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  18 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  19 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  24 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  25 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  26 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  27 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  42 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  44 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  46 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  48 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  62 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  68 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  82 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  26 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  27 );
             tBasis[  88 ] = tNeighbor->get_basis_function(   7 );
         }

         // get pointer to neighbor  10
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 10 );

         // test if neighbor  10 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  36 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  37 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  38 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  39 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  42 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  44 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  46 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  56 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  57 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  62 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  76 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  77 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  82 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  96 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  97 ] = tNeighbor->get_basis_function(  24 );
             tBasis[  98 ] = tNeighbor->get_basis_function(  25 );
             tBasis[  99 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 102 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(   7 );
         }

         // get pointer to neighbor  11
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 11 );

         // test if neighbor  11 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  38 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  39 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  40 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  41 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  43 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  45 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  47 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  60 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  61 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  63 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  80 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  81 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  83 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  98 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  99 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 100 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 101 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 103 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  12
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 12 );

         // test if neighbor  12 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  45 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  47 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  49 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  52 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  53 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  54 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  55 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  69 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  74 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  75 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  89 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  94 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  95 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 109 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 114 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 115 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  13
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 13 );

         // test if neighbor  13 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  44 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  46 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  48 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  50 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  51 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  52 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  53 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  68 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  70 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  71 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  88 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  90 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  91 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 108 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 110 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 111 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  14
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 14 );

         // test if neighbor  14 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  57 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  58 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  59 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  60 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  77 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  80 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  97 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  98 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  99 ] = tNeighbor->get_basis_function(  38 );
             tBasis[ 100 ] = tNeighbor->get_basis_function(  17 );
             tBasis[ 117 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 118 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 119 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 120 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 123 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 126 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  15
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 15 );

         // test if neighbor  15 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  63 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  69 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  83 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  89 ] = tNeighbor->get_basis_function(  20 );
             tBasis[ 103 ] = tNeighbor->get_basis_function(  17 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(  47 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(  46 );
             tBasis[ 109 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 126 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 127 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 133 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 139 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 144 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 145 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  16
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 16 );

         // test if neighbor  16 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  71 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  74 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  91 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  94 ] = tNeighbor->get_basis_function(  20 );
             tBasis[ 111 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(  50 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(  51 );
             tBasis[ 114 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 141 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 144 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 147 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 148 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 149 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 150 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  17
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 17 );

         // test if neighbor  17 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  62 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  68 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  82 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  88 ] = tNeighbor->get_basis_function(  22 );
             tBasis[ 102 ] = tNeighbor->get_basis_function(  13 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(  41 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(  42 );
             tBasis[ 108 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 122 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 123 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 128 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 134 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 140 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 141 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  18
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 18 );

         // test if neighbor  18 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   0 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   1 ] = tNeighbor->get_basis_function(   8 );
             tBasis[   2 ] = tNeighbor->get_basis_function(   9 );
             tBasis[   3 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   6 ] = tNeighbor->get_basis_function(  10 );
             tBasis[   7 ] = tNeighbor->get_basis_function(  32 );
             tBasis[   8 ] = tNeighbor->get_basis_function(  35 );
             tBasis[   9 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  12 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  13 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  18 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  19 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  21 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  36 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  37 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  38 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  39 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  42 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  44 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  46 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  56 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  57 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  62 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  76 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  77 ] = tNeighbor->get_basis_function(  24 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  25 );
             tBasis[  79 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  82 ] = tNeighbor->get_basis_function(  26 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  27 );
             tBasis[  86 ] = tNeighbor->get_basis_function(   7 );
         }

         // get pointer to neighbor  19
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 19 );

         // test if neighbor  19 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   2 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   3 ] = tNeighbor->get_basis_function(   8 );
             tBasis[   4 ] = tNeighbor->get_basis_function(   9 );
             tBasis[   5 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   8 ] = tNeighbor->get_basis_function(  10 );
             tBasis[   9 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  10 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  11 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  14 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  15 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  16 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  17 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  20 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  22 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  23 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  38 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  39 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  40 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  41 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  43 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  45 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  47 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  58 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  59 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  60 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  61 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  63 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  78 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  24 );
             tBasis[  80 ] = tNeighbor->get_basis_function(  25 );
             tBasis[  81 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  83 ] = tNeighbor->get_basis_function(  28 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  29 );
             tBasis[  87 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  20
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 20 );

         // test if neighbor  20 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  14 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  15 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  16 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  17 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  22 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  23 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  26 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  27 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  28 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  29 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  32 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  33 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  34 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  35 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  45 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  47 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  49 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  52 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  53 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  54 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  55 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  17 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  47 );
             tBasis[  69 ] = tNeighbor->get_basis_function(  46 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  74 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  75 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  85 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  28 );
             tBasis[  89 ] = tNeighbor->get_basis_function(  29 );
             tBasis[  92 ] = tNeighbor->get_basis_function(   7 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  31 );
             tBasis[  94 ] = tNeighbor->get_basis_function(  30 );
             tBasis[  95 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  21
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 21 );

         // test if neighbor  21 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  12 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  13 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  14 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  15 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  18 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  19 ] = tNeighbor->get_basis_function(  32 );
             tBasis[  20 ] = tNeighbor->get_basis_function(  35 );
             tBasis[  21 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  24 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  25 ] = tNeighbor->get_basis_function(  33 );
             tBasis[  26 ] = tNeighbor->get_basis_function(  34 );
             tBasis[  27 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  30 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  31 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  32 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  33 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  44 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  46 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  48 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  50 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  51 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  52 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  53 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  41 );
             tBasis[  68 ] = tNeighbor->get_basis_function(  42 );
             tBasis[  70 ] = tNeighbor->get_basis_function(  23 );
             tBasis[  71 ] = tNeighbor->get_basis_function(  50 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  51 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  21 );
             tBasis[  84 ] = tNeighbor->get_basis_function(   4 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  26 );
             tBasis[  88 ] = tNeighbor->get_basis_function(  27 );
             tBasis[  90 ] = tNeighbor->get_basis_function(   7 );
             tBasis[  91 ] = tNeighbor->get_basis_function(  31 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  30 );
             tBasis[  93 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  22
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 22 );

         // test if neighbor  22 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  56 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  57 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  58 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  59 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  62 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  64 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  66 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  76 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  77 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  82 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  96 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  97 ] = tNeighbor->get_basis_function(  39 );
             tBasis[  98 ] = tNeighbor->get_basis_function(  38 );
             tBasis[  99 ] = tNeighbor->get_basis_function(  17 );
             tBasis[ 102 ] = tNeighbor->get_basis_function(  41 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(  42 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 116 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 117 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 118 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 119 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 122 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 123 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 128 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 134 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  23
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 23 );

         // test if neighbor  23 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  58 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  59 ] = tNeighbor->get_basis_function(   8 );
             tBasis[  60 ] = tNeighbor->get_basis_function(   9 );
             tBasis[  61 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  63 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  65 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  67 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  78 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  79 ] = tNeighbor->get_basis_function(  36 );
             tBasis[  80 ] = tNeighbor->get_basis_function(  37 );
             tBasis[  81 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  83 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  20 );
             tBasis[  98 ] = tNeighbor->get_basis_function(  13 );
             tBasis[  99 ] = tNeighbor->get_basis_function(  39 );
             tBasis[ 100 ] = tNeighbor->get_basis_function(  38 );
             tBasis[ 101 ] = tNeighbor->get_basis_function(  17 );
             tBasis[ 103 ] = tNeighbor->get_basis_function(  47 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(  46 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 118 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 119 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 120 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 121 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 124 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 125 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 126 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 127 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 133 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 139 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  24
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 24 );

         // test if neighbor  24 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  65 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  67 ] = tNeighbor->get_basis_function(  14 );
             tBasis[  69 ] = tNeighbor->get_basis_function(  15 );
             tBasis[  72 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  73 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  74 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  75 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  85 ] = tNeighbor->get_basis_function(  16 );
             tBasis[  87 ] = tNeighbor->get_basis_function(  44 );
             tBasis[  89 ] = tNeighbor->get_basis_function(  45 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  94 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  95 ] = tNeighbor->get_basis_function(  20 );
             tBasis[ 105 ] = tNeighbor->get_basis_function(  17 );
             tBasis[ 107 ] = tNeighbor->get_basis_function(  47 );
             tBasis[ 109 ] = tNeighbor->get_basis_function(  46 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(  50 );
             tBasis[ 114 ] = tNeighbor->get_basis_function(  51 );
             tBasis[ 115 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 132 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 133 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 138 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 139 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 144 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 145 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 148 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 149 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 150 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 151 ] = tNeighbor->get_basis_function(   6 );
         }

         // get pointer to neighbor  25
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 25 );

         // test if neighbor  25 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  64 ] = tNeighbor->get_basis_function(   0 );
             tBasis[  66 ] = tNeighbor->get_basis_function(  10 );
             tBasis[  68 ] = tNeighbor->get_basis_function(  11 );
             tBasis[  70 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  71 ] = tNeighbor->get_basis_function(  19 );
             tBasis[  72 ] = tNeighbor->get_basis_function(  18 );
             tBasis[  73 ] = tNeighbor->get_basis_function(   2 );
             tBasis[  84 ] = tNeighbor->get_basis_function(  12 );
             tBasis[  86 ] = tNeighbor->get_basis_function(  40 );
             tBasis[  88 ] = tNeighbor->get_basis_function(  43 );
             tBasis[  90 ] = tNeighbor->get_basis_function(  22 );
             tBasis[  91 ] = tNeighbor->get_basis_function(  49 );
             tBasis[  92 ] = tNeighbor->get_basis_function(  48 );
             tBasis[  93 ] = tNeighbor->get_basis_function(  20 );
             tBasis[ 104 ] = tNeighbor->get_basis_function(  13 );
             tBasis[ 106 ] = tNeighbor->get_basis_function(  41 );
             tBasis[ 108 ] = tNeighbor->get_basis_function(  42 );
             tBasis[ 110 ] = tNeighbor->get_basis_function(  23 );
             tBasis[ 111 ] = tNeighbor->get_basis_function(  50 );
             tBasis[ 112 ] = tNeighbor->get_basis_function(  51 );
             tBasis[ 113 ] = tNeighbor->get_basis_function(  21 );
             tBasis[ 128 ] = tNeighbor->get_basis_function(   4 );
             tBasis[ 129 ] = tNeighbor->get_basis_function(  24 );
             tBasis[ 130 ] = tNeighbor->get_basis_function(  25 );
             tBasis[ 131 ] = tNeighbor->get_basis_function(   5 );
             tBasis[ 134 ] = tNeighbor->get_basis_function(  26 );
             tBasis[ 135 ] = tNeighbor->get_basis_function(  52 );
             tBasis[ 136 ] = tNeighbor->get_basis_function(  53 );
             tBasis[ 137 ] = tNeighbor->get_basis_function(  28 );
             tBasis[ 140 ] = tNeighbor->get_basis_function(  27 );
             tBasis[ 141 ] = tNeighbor->get_basis_function(  55 );
             tBasis[ 142 ] = tNeighbor->get_basis_function(  54 );
             tBasis[ 143 ] = tNeighbor->get_basis_function(  29 );
             tBasis[ 146 ] = tNeighbor->get_basis_function(   7 );
             tBasis[ 147 ] = tNeighbor->get_basis_function(  31 );
             tBasis[ 148 ] = tNeighbor->get_basis_function(  30 );
             tBasis[ 149 ] = tNeighbor->get_basis_function(   6 );
         }

         // test if basis 0 exists
         if ( mBasis[   0 ] != nullptr )
         {
             // test if basis 0 has been processed
             if ( ! mBasis[   0 ]->is_flagged() )
             {
                 // link neighbors of basis 0
                 mBasis[   0 ]->insert_neighbor(  0, tBasis[  37 ] );
                 mBasis[   0 ]->insert_neighbor(  1, mBasis[   8 ] );
                 mBasis[   0 ]->insert_neighbor(  2, mBasis[  10 ] );
                 mBasis[   0 ]->insert_neighbor(  3, tBasis[  42 ] );
                 mBasis[   0 ]->insert_neighbor(  4, tBasis[   7 ] );
                 mBasis[   0 ]->insert_neighbor(  5, mBasis[  12 ] );
                 mBasis[   0 ]->insert_neighbor(  6, tBasis[   1 ] );
                 mBasis[   0 ]->insert_neighbor(  7, tBasis[   8 ] );
                 mBasis[   0 ]->insert_neighbor(  8, tBasis[  13 ] );
                 mBasis[   0 ]->insert_neighbor(  9, tBasis[   6 ] );
                 mBasis[   0 ]->insert_neighbor( 10, tBasis[  36 ] );
                 mBasis[   0 ]->insert_neighbor( 11, tBasis[  38 ] );
                 mBasis[   0 ]->insert_neighbor( 12, mBasis[  32 ] );
                 mBasis[   0 ]->insert_neighbor( 13, tBasis[  44 ] );
                 mBasis[   0 ]->insert_neighbor( 14, tBasis[  57 ] );
                 mBasis[   0 ]->insert_neighbor( 15, mBasis[  36 ] );
                 mBasis[   0 ]->insert_neighbor( 16, mBasis[  40 ] );
                 mBasis[   0 ]->insert_neighbor( 17, tBasis[  62 ] );
                 mBasis[   0 ]->insert_neighbor( 18, tBasis[   0 ] );
                 mBasis[   0 ]->insert_neighbor( 19, tBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor( 20, tBasis[  14 ] );
                 mBasis[   0 ]->insert_neighbor( 21, tBasis[  12 ] );
                 mBasis[   0 ]->insert_neighbor( 22, tBasis[  56 ] );
                 mBasis[   0 ]->insert_neighbor( 23, tBasis[  58 ] );
                 mBasis[   0 ]->insert_neighbor( 24, mBasis[  56 ] );
                 mBasis[   0 ]->insert_neighbor( 25, tBasis[  64 ] );

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
                 mBasis[   1 ]->insert_neighbor(  0, tBasis[  40 ] );
                 mBasis[   1 ]->insert_neighbor(  1, tBasis[  43 ] );
                 mBasis[   1 ]->insert_neighbor(  2, mBasis[  14 ] );
                 mBasis[   1 ]->insert_neighbor(  3, mBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor(  4, tBasis[  10 ] );
                 mBasis[   1 ]->insert_neighbor(  5, mBasis[  16 ] );
                 mBasis[   1 ]->insert_neighbor(  6, tBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor(  7, tBasis[  11 ] );
                 mBasis[   1 ]->insert_neighbor(  8, tBasis[  16 ] );
                 mBasis[   1 ]->insert_neighbor(  9, tBasis[   9 ] );
                 mBasis[   1 ]->insert_neighbor( 10, tBasis[  39 ] );
                 mBasis[   1 ]->insert_neighbor( 11, tBasis[  41 ] );
                 mBasis[   1 ]->insert_neighbor( 12, tBasis[  45 ] );
                 mBasis[   1 ]->insert_neighbor( 13, mBasis[  35 ] );
                 mBasis[   1 ]->insert_neighbor( 14, tBasis[  60 ] );
                 mBasis[   1 ]->insert_neighbor( 15, tBasis[  63 ] );
                 mBasis[   1 ]->insert_neighbor( 16, mBasis[  44 ] );
                 mBasis[   1 ]->insert_neighbor( 17, mBasis[  37 ] );
                 mBasis[   1 ]->insert_neighbor( 18, tBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor( 19, tBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor( 20, tBasis[  17 ] );
                 mBasis[   1 ]->insert_neighbor( 21, tBasis[  15 ] );
                 mBasis[   1 ]->insert_neighbor( 22, tBasis[  59 ] );
                 mBasis[   1 ]->insert_neighbor( 23, tBasis[  61 ] );
                 mBasis[   1 ]->insert_neighbor( 24, tBasis[  65 ] );
                 mBasis[   1 ]->insert_neighbor( 25, mBasis[  57 ] );

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
                 mBasis[   2 ]->insert_neighbor(  0, mBasis[  15 ] );
                 mBasis[   2 ]->insert_neighbor(  1, tBasis[  49 ] );
                 mBasis[   2 ]->insert_neighbor(  2, tBasis[  54 ] );
                 mBasis[   2 ]->insert_neighbor(  3, mBasis[  18 ] );
                 mBasis[   2 ]->insert_neighbor(  4, tBasis[  28 ] );
                 mBasis[   2 ]->insert_neighbor(  5, mBasis[  20 ] );
                 mBasis[   2 ]->insert_neighbor(  6, tBasis[  22 ] );
                 mBasis[   2 ]->insert_neighbor(  7, tBasis[  29 ] );
                 mBasis[   2 ]->insert_neighbor(  8, tBasis[  34 ] );
                 mBasis[   2 ]->insert_neighbor(  9, tBasis[  27 ] );
                 mBasis[   2 ]->insert_neighbor( 10, mBasis[  34 ] );
                 mBasis[   2 ]->insert_neighbor( 11, tBasis[  47 ] );
                 mBasis[   2 ]->insert_neighbor( 12, tBasis[  55 ] );
                 mBasis[   2 ]->insert_neighbor( 13, tBasis[  53 ] );
                 mBasis[   2 ]->insert_neighbor( 14, mBasis[  45 ] );
                 mBasis[   2 ]->insert_neighbor( 15, tBasis[  69 ] );
                 mBasis[   2 ]->insert_neighbor( 16, tBasis[  74 ] );
                 mBasis[   2 ]->insert_neighbor( 17, mBasis[  48 ] );
                 mBasis[   2 ]->insert_neighbor( 18, tBasis[  21 ] );
                 mBasis[   2 ]->insert_neighbor( 19, tBasis[  23 ] );
                 mBasis[   2 ]->insert_neighbor( 20, tBasis[  35 ] );
                 mBasis[   2 ]->insert_neighbor( 21, tBasis[  33 ] );
                 mBasis[   2 ]->insert_neighbor( 22, mBasis[  58 ] );
                 mBasis[   2 ]->insert_neighbor( 23, tBasis[  67 ] );
                 mBasis[   2 ]->insert_neighbor( 24, tBasis[  75 ] );
                 mBasis[   2 ]->insert_neighbor( 25, tBasis[  73 ] );

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
                 mBasis[   3 ]->insert_neighbor(  1, mBasis[  19 ] );
                 mBasis[   3 ]->insert_neighbor(  2, tBasis[  51 ] );
                 mBasis[   3 ]->insert_neighbor(  3, tBasis[  48 ] );
                 mBasis[   3 ]->insert_neighbor(  4, tBasis[  25 ] );
                 mBasis[   3 ]->insert_neighbor(  5, mBasis[  22 ] );
                 mBasis[   3 ]->insert_neighbor(  6, tBasis[  19 ] );
                 mBasis[   3 ]->insert_neighbor(  7, tBasis[  26 ] );
                 mBasis[   3 ]->insert_neighbor(  8, tBasis[  31 ] );
                 mBasis[   3 ]->insert_neighbor(  9, tBasis[  24 ] );
                 mBasis[   3 ]->insert_neighbor( 10, tBasis[  46 ] );
                 mBasis[   3 ]->insert_neighbor( 11, mBasis[  33 ] );
                 mBasis[   3 ]->insert_neighbor( 12, tBasis[  52 ] );
                 mBasis[   3 ]->insert_neighbor( 13, tBasis[  50 ] );
                 mBasis[   3 ]->insert_neighbor( 14, mBasis[  43 ] );
                 mBasis[   3 ]->insert_neighbor( 15, mBasis[  49 ] );
                 mBasis[   3 ]->insert_neighbor( 16, tBasis[  71 ] );
                 mBasis[   3 ]->insert_neighbor( 17, tBasis[  68 ] );
                 mBasis[   3 ]->insert_neighbor( 18, tBasis[  18 ] );
                 mBasis[   3 ]->insert_neighbor( 19, tBasis[  20 ] );
                 mBasis[   3 ]->insert_neighbor( 20, tBasis[  32 ] );
                 mBasis[   3 ]->insert_neighbor( 21, tBasis[  30 ] );
                 mBasis[   3 ]->insert_neighbor( 22, tBasis[  66 ] );
                 mBasis[   3 ]->insert_neighbor( 23, mBasis[  59 ] );
                 mBasis[   3 ]->insert_neighbor( 24, tBasis[  72 ] );
                 mBasis[   3 ]->insert_neighbor( 25, tBasis[  70 ] );

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
                 mBasis[   4 ]->insert_neighbor(  0, tBasis[  97 ] );
                 mBasis[   4 ]->insert_neighbor(  1, mBasis[  24 ] );
                 mBasis[   4 ]->insert_neighbor(  2, mBasis[  26 ] );
                 mBasis[   4 ]->insert_neighbor(  3, tBasis[ 102 ] );
                 mBasis[   4 ]->insert_neighbor(  4, mBasis[  13 ] );
                 mBasis[   4 ]->insert_neighbor(  5, tBasis[ 123 ] );
                 mBasis[   4 ]->insert_neighbor(  6, tBasis[  77 ] );
                 mBasis[   4 ]->insert_neighbor(  7, mBasis[  39 ] );
                 mBasis[   4 ]->insert_neighbor(  8, mBasis[  41 ] );
                 mBasis[   4 ]->insert_neighbor(  9, tBasis[  82 ] );
                 mBasis[   4 ]->insert_neighbor( 10, tBasis[  96 ] );
                 mBasis[   4 ]->insert_neighbor( 11, tBasis[  98 ] );
                 mBasis[   4 ]->insert_neighbor( 12, mBasis[  52 ] );
                 mBasis[   4 ]->insert_neighbor( 13, tBasis[ 104 ] );
                 mBasis[   4 ]->insert_neighbor( 14, tBasis[ 117 ] );
                 mBasis[   4 ]->insert_neighbor( 15, tBasis[ 124 ] );
                 mBasis[   4 ]->insert_neighbor( 16, tBasis[ 129 ] );
                 mBasis[   4 ]->insert_neighbor( 17, tBasis[ 122 ] );
                 mBasis[   4 ]->insert_neighbor( 18, tBasis[  76 ] );
                 mBasis[   4 ]->insert_neighbor( 19, tBasis[  78 ] );
                 mBasis[   4 ]->insert_neighbor( 20, mBasis[  60 ] );
                 mBasis[   4 ]->insert_neighbor( 21, tBasis[  84 ] );
                 mBasis[   4 ]->insert_neighbor( 22, tBasis[ 116 ] );
                 mBasis[   4 ]->insert_neighbor( 23, tBasis[ 118 ] );
                 mBasis[   4 ]->insert_neighbor( 24, tBasis[ 130 ] );
                 mBasis[   4 ]->insert_neighbor( 25, tBasis[ 128 ] );

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
                 mBasis[   5 ]->insert_neighbor(  0, tBasis[ 100 ] );
                 mBasis[   5 ]->insert_neighbor(  1, tBasis[ 103 ] );
                 mBasis[   5 ]->insert_neighbor(  2, mBasis[  28 ] );
                 mBasis[   5 ]->insert_neighbor(  3, mBasis[  25 ] );
                 mBasis[   5 ]->insert_neighbor(  4, mBasis[  17 ] );
                 mBasis[   5 ]->insert_neighbor(  5, tBasis[ 126 ] );
                 mBasis[   5 ]->insert_neighbor(  6, tBasis[  80 ] );
                 mBasis[   5 ]->insert_neighbor(  7, tBasis[  83 ] );
                 mBasis[   5 ]->insert_neighbor(  8, mBasis[  47 ] );
                 mBasis[   5 ]->insert_neighbor(  9, mBasis[  38 ] );
                 mBasis[   5 ]->insert_neighbor( 10, tBasis[  99 ] );
                 mBasis[   5 ]->insert_neighbor( 11, tBasis[ 101 ] );
                 mBasis[   5 ]->insert_neighbor( 12, tBasis[ 105 ] );
                 mBasis[   5 ]->insert_neighbor( 13, mBasis[  53 ] );
                 mBasis[   5 ]->insert_neighbor( 14, tBasis[ 120 ] );
                 mBasis[   5 ]->insert_neighbor( 15, tBasis[ 127 ] );
                 mBasis[   5 ]->insert_neighbor( 16, tBasis[ 132 ] );
                 mBasis[   5 ]->insert_neighbor( 17, tBasis[ 125 ] );
                 mBasis[   5 ]->insert_neighbor( 18, tBasis[  79 ] );
                 mBasis[   5 ]->insert_neighbor( 19, tBasis[  81 ] );
                 mBasis[   5 ]->insert_neighbor( 20, tBasis[  85 ] );
                 mBasis[   5 ]->insert_neighbor( 21, mBasis[  61 ] );
                 mBasis[   5 ]->insert_neighbor( 22, tBasis[ 119 ] );
                 mBasis[   5 ]->insert_neighbor( 23, tBasis[ 121 ] );
                 mBasis[   5 ]->insert_neighbor( 24, tBasis[ 133 ] );
                 mBasis[   5 ]->insert_neighbor( 25, tBasis[ 131 ] );

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
                 mBasis[   6 ]->insert_neighbor(  0, mBasis[  29 ] );
                 mBasis[   6 ]->insert_neighbor(  1, tBasis[ 109 ] );
                 mBasis[   6 ]->insert_neighbor(  2, tBasis[ 114 ] );
                 mBasis[   6 ]->insert_neighbor(  3, mBasis[  30 ] );
                 mBasis[   6 ]->insert_neighbor(  4, mBasis[  21 ] );
                 mBasis[   6 ]->insert_neighbor(  5, tBasis[ 144 ] );
                 mBasis[   6 ]->insert_neighbor(  6, mBasis[  46 ] );
                 mBasis[   6 ]->insert_neighbor(  7, tBasis[  89 ] );
                 mBasis[   6 ]->insert_neighbor(  8, tBasis[  94 ] );
                 mBasis[   6 ]->insert_neighbor(  9, mBasis[  51 ] );
                 mBasis[   6 ]->insert_neighbor( 10, mBasis[  54 ] );
                 mBasis[   6 ]->insert_neighbor( 11, tBasis[ 107 ] );
                 mBasis[   6 ]->insert_neighbor( 12, tBasis[ 115 ] );
                 mBasis[   6 ]->insert_neighbor( 13, tBasis[ 113 ] );
                 mBasis[   6 ]->insert_neighbor( 14, tBasis[ 138 ] );
                 mBasis[   6 ]->insert_neighbor( 15, tBasis[ 145 ] );
                 mBasis[   6 ]->insert_neighbor( 16, tBasis[ 150 ] );
                 mBasis[   6 ]->insert_neighbor( 17, tBasis[ 143 ] );
                 mBasis[   6 ]->insert_neighbor( 18, mBasis[  62 ] );
                 mBasis[   6 ]->insert_neighbor( 19, tBasis[  87 ] );
                 mBasis[   6 ]->insert_neighbor( 20, tBasis[  95 ] );
                 mBasis[   6 ]->insert_neighbor( 21, tBasis[  93 ] );
                 mBasis[   6 ]->insert_neighbor( 22, tBasis[ 137 ] );
                 mBasis[   6 ]->insert_neighbor( 23, tBasis[ 139 ] );
                 mBasis[   6 ]->insert_neighbor( 24, tBasis[ 151 ] );
                 mBasis[   6 ]->insert_neighbor( 25, tBasis[ 149 ] );

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
                 mBasis[   7 ]->insert_neighbor(  0, mBasis[  27 ] );
                 mBasis[   7 ]->insert_neighbor(  1, mBasis[  31 ] );
                 mBasis[   7 ]->insert_neighbor(  2, tBasis[ 111 ] );
                 mBasis[   7 ]->insert_neighbor(  3, tBasis[ 108 ] );
                 mBasis[   7 ]->insert_neighbor(  4, mBasis[  23 ] );
                 mBasis[   7 ]->insert_neighbor(  5, tBasis[ 141 ] );
                 mBasis[   7 ]->insert_neighbor(  6, mBasis[  42 ] );
                 mBasis[   7 ]->insert_neighbor(  7, mBasis[  50 ] );
                 mBasis[   7 ]->insert_neighbor(  8, tBasis[  91 ] );
                 mBasis[   7 ]->insert_neighbor(  9, tBasis[  88 ] );
                 mBasis[   7 ]->insert_neighbor( 10, tBasis[ 106 ] );
                 mBasis[   7 ]->insert_neighbor( 11, mBasis[  55 ] );
                 mBasis[   7 ]->insert_neighbor( 12, tBasis[ 112 ] );
                 mBasis[   7 ]->insert_neighbor( 13, tBasis[ 110 ] );
                 mBasis[   7 ]->insert_neighbor( 14, tBasis[ 135 ] );
                 mBasis[   7 ]->insert_neighbor( 15, tBasis[ 142 ] );
                 mBasis[   7 ]->insert_neighbor( 16, tBasis[ 147 ] );
                 mBasis[   7 ]->insert_neighbor( 17, tBasis[ 140 ] );
                 mBasis[   7 ]->insert_neighbor( 18, tBasis[  86 ] );
                 mBasis[   7 ]->insert_neighbor( 19, mBasis[  63 ] );
                 mBasis[   7 ]->insert_neighbor( 20, tBasis[  92 ] );
                 mBasis[   7 ]->insert_neighbor( 21, tBasis[  90 ] );
                 mBasis[   7 ]->insert_neighbor( 22, tBasis[ 134 ] );
                 mBasis[   7 ]->insert_neighbor( 23, tBasis[ 136 ] );
                 mBasis[   7 ]->insert_neighbor( 24, tBasis[ 148 ] );
                 mBasis[   7 ]->insert_neighbor( 25, tBasis[ 146 ] );

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
                 mBasis[   8 ]->insert_neighbor(  0, tBasis[  38 ] );
                 mBasis[   8 ]->insert_neighbor(  1, mBasis[   9 ] );
                 mBasis[   8 ]->insert_neighbor(  2, mBasis[  32 ] );
                 mBasis[   8 ]->insert_neighbor(  3, mBasis[   0 ] );
                 mBasis[   8 ]->insert_neighbor(  4, tBasis[   8 ] );
                 mBasis[   8 ]->insert_neighbor(  5, mBasis[  36 ] );
                 mBasis[   8 ]->insert_neighbor(  6, tBasis[   2 ] );
                 mBasis[   8 ]->insert_neighbor(  7, tBasis[   9 ] );
                 mBasis[   8 ]->insert_neighbor(  8, tBasis[  14 ] );
                 mBasis[   8 ]->insert_neighbor(  9, tBasis[   7 ] );
                 mBasis[   8 ]->insert_neighbor( 10, tBasis[  37 ] );
                 mBasis[   8 ]->insert_neighbor( 11, tBasis[  39 ] );
                 mBasis[   8 ]->insert_neighbor( 12, mBasis[  35 ] );
                 mBasis[   8 ]->insert_neighbor( 13, mBasis[  10 ] );
                 mBasis[   8 ]->insert_neighbor( 14, tBasis[  58 ] );
                 mBasis[   8 ]->insert_neighbor( 15, mBasis[  37 ] );
                 mBasis[   8 ]->insert_neighbor( 16, mBasis[  56 ] );
                 mBasis[   8 ]->insert_neighbor( 17, mBasis[  12 ] );
                 mBasis[   8 ]->insert_neighbor( 18, tBasis[   1 ] );
                 mBasis[   8 ]->insert_neighbor( 19, tBasis[   3 ] );
                 mBasis[   8 ]->insert_neighbor( 20, tBasis[  15 ] );
                 mBasis[   8 ]->insert_neighbor( 21, tBasis[  13 ] );
                 mBasis[   8 ]->insert_neighbor( 22, tBasis[  57 ] );
                 mBasis[   8 ]->insert_neighbor( 23, tBasis[  59 ] );
                 mBasis[   8 ]->insert_neighbor( 24, mBasis[  57 ] );
                 mBasis[   8 ]->insert_neighbor( 25, mBasis[  40 ] );

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
                 mBasis[   9 ]->insert_neighbor(  0, tBasis[  39 ] );
                 mBasis[   9 ]->insert_neighbor(  1, mBasis[   1 ] );
                 mBasis[   9 ]->insert_neighbor(  2, mBasis[  35 ] );
                 mBasis[   9 ]->insert_neighbor(  3, mBasis[   8 ] );
                 mBasis[   9 ]->insert_neighbor(  4, tBasis[   9 ] );
                 mBasis[   9 ]->insert_neighbor(  5, mBasis[  37 ] );
                 mBasis[   9 ]->insert_neighbor(  6, tBasis[   3 ] );
                 mBasis[   9 ]->insert_neighbor(  7, tBasis[  10 ] );
                 mBasis[   9 ]->insert_neighbor(  8, tBasis[  15 ] );
                 mBasis[   9 ]->insert_neighbor(  9, tBasis[   8 ] );
                 mBasis[   9 ]->insert_neighbor( 10, tBasis[  38 ] );
                 mBasis[   9 ]->insert_neighbor( 11, tBasis[  40 ] );
                 mBasis[   9 ]->insert_neighbor( 12, mBasis[  14 ] );
                 mBasis[   9 ]->insert_neighbor( 13, mBasis[  32 ] );
                 mBasis[   9 ]->insert_neighbor( 14, tBasis[  59 ] );
                 mBasis[   9 ]->insert_neighbor( 15, mBasis[  16 ] );
                 mBasis[   9 ]->insert_neighbor( 16, mBasis[  57 ] );
                 mBasis[   9 ]->insert_neighbor( 17, mBasis[  36 ] );
                 mBasis[   9 ]->insert_neighbor( 18, tBasis[   2 ] );
                 mBasis[   9 ]->insert_neighbor( 19, tBasis[   4 ] );
                 mBasis[   9 ]->insert_neighbor( 20, tBasis[  16 ] );
                 mBasis[   9 ]->insert_neighbor( 21, tBasis[  14 ] );
                 mBasis[   9 ]->insert_neighbor( 22, tBasis[  58 ] );
                 mBasis[   9 ]->insert_neighbor( 23, tBasis[  60 ] );
                 mBasis[   9 ]->insert_neighbor( 24, mBasis[  44 ] );
                 mBasis[   9 ]->insert_neighbor( 25, mBasis[  56 ] );

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
                 mBasis[  10 ]->insert_neighbor(  0, mBasis[   0 ] );
                 mBasis[  10 ]->insert_neighbor(  1, mBasis[  32 ] );
                 mBasis[  10 ]->insert_neighbor(  2, mBasis[  11 ] );
                 mBasis[  10 ]->insert_neighbor(  3, tBasis[  44 ] );
                 mBasis[  10 ]->insert_neighbor(  4, tBasis[  13 ] );
                 mBasis[  10 ]->insert_neighbor(  5, mBasis[  40 ] );
                 mBasis[  10 ]->insert_neighbor(  6, tBasis[   7 ] );
                 mBasis[  10 ]->insert_neighbor(  7, tBasis[  14 ] );
                 mBasis[  10 ]->insert_neighbor(  8, tBasis[  19 ] );
                 mBasis[  10 ]->insert_neighbor(  9, tBasis[  12 ] );
                 mBasis[  10 ]->insert_neighbor( 10, tBasis[  42 ] );
                 mBasis[  10 ]->insert_neighbor( 11, mBasis[   8 ] );
                 mBasis[  10 ]->insert_neighbor( 12, mBasis[  33 ] );
                 mBasis[  10 ]->insert_neighbor( 13, tBasis[  46 ] );
                 mBasis[  10 ]->insert_neighbor( 14, mBasis[  12 ] );
                 mBasis[  10 ]->insert_neighbor( 15, mBasis[  56 ] );
                 mBasis[  10 ]->insert_neighbor( 16, mBasis[  43 ] );
                 mBasis[  10 ]->insert_neighbor( 17, tBasis[  64 ] );
                 mBasis[  10 ]->insert_neighbor( 18, tBasis[   6 ] );
                 mBasis[  10 ]->insert_neighbor( 19, tBasis[   8 ] );
                 mBasis[  10 ]->insert_neighbor( 20, tBasis[  20 ] );
                 mBasis[  10 ]->insert_neighbor( 21, tBasis[  18 ] );
                 mBasis[  10 ]->insert_neighbor( 22, tBasis[  62 ] );
                 mBasis[  10 ]->insert_neighbor( 23, mBasis[  36 ] );
                 mBasis[  10 ]->insert_neighbor( 24, mBasis[  59 ] );
                 mBasis[  10 ]->insert_neighbor( 25, tBasis[  66 ] );

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
                 mBasis[  11 ]->insert_neighbor(  0, mBasis[  10 ] );
                 mBasis[  11 ]->insert_neighbor(  1, mBasis[  33 ] );
                 mBasis[  11 ]->insert_neighbor(  2, mBasis[   3 ] );
                 mBasis[  11 ]->insert_neighbor(  3, tBasis[  46 ] );
                 mBasis[  11 ]->insert_neighbor(  4, tBasis[  19 ] );
                 mBasis[  11 ]->insert_neighbor(  5, mBasis[  43 ] );
                 mBasis[  11 ]->insert_neighbor(  6, tBasis[  13 ] );
                 mBasis[  11 ]->insert_neighbor(  7, tBasis[  20 ] );
                 mBasis[  11 ]->insert_neighbor(  8, tBasis[  25 ] );
                 mBasis[  11 ]->insert_neighbor(  9, tBasis[  18 ] );
                 mBasis[  11 ]->insert_neighbor( 10, tBasis[  44 ] );
                 mBasis[  11 ]->insert_neighbor( 11, mBasis[  32 ] );
                 mBasis[  11 ]->insert_neighbor( 12, mBasis[  19 ] );
                 mBasis[  11 ]->insert_neighbor( 13, tBasis[  48 ] );
                 mBasis[  11 ]->insert_neighbor( 14, mBasis[  40 ] );
                 mBasis[  11 ]->insert_neighbor( 15, mBasis[  59 ] );
                 mBasis[  11 ]->insert_neighbor( 16, mBasis[  22 ] );
                 mBasis[  11 ]->insert_neighbor( 17, tBasis[  66 ] );
                 mBasis[  11 ]->insert_neighbor( 18, tBasis[  12 ] );
                 mBasis[  11 ]->insert_neighbor( 19, tBasis[  14 ] );
                 mBasis[  11 ]->insert_neighbor( 20, tBasis[  26 ] );
                 mBasis[  11 ]->insert_neighbor( 21, tBasis[  24 ] );
                 mBasis[  11 ]->insert_neighbor( 22, tBasis[  64 ] );
                 mBasis[  11 ]->insert_neighbor( 23, mBasis[  56 ] );
                 mBasis[  11 ]->insert_neighbor( 24, mBasis[  49 ] );
                 mBasis[  11 ]->insert_neighbor( 25, tBasis[  68 ] );

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
                 mBasis[  12 ]->insert_neighbor(  0, tBasis[  57 ] );
                 mBasis[  12 ]->insert_neighbor(  1, mBasis[  36 ] );
                 mBasis[  12 ]->insert_neighbor(  2, mBasis[  40 ] );
                 mBasis[  12 ]->insert_neighbor(  3, tBasis[  62 ] );
                 mBasis[  12 ]->insert_neighbor(  4, mBasis[   0 ] );
                 mBasis[  12 ]->insert_neighbor(  5, mBasis[  13 ] );
                 mBasis[  12 ]->insert_neighbor(  6, tBasis[  37 ] );
                 mBasis[  12 ]->insert_neighbor(  7, mBasis[   8 ] );
                 mBasis[  12 ]->insert_neighbor(  8, mBasis[  10 ] );
                 mBasis[  12 ]->insert_neighbor(  9, tBasis[  42 ] );
                 mBasis[  12 ]->insert_neighbor( 10, tBasis[  56 ] );
                 mBasis[  12 ]->insert_neighbor( 11, tBasis[  58 ] );
                 mBasis[  12 ]->insert_neighbor( 12, mBasis[  56 ] );
                 mBasis[  12 ]->insert_neighbor( 13, tBasis[  64 ] );
                 mBasis[  12 ]->insert_neighbor( 14, tBasis[  77 ] );
                 mBasis[  12 ]->insert_neighbor( 15, mBasis[  39 ] );
                 mBasis[  12 ]->insert_neighbor( 16, mBasis[  41 ] );
                 mBasis[  12 ]->insert_neighbor( 17, tBasis[  82 ] );
                 mBasis[  12 ]->insert_neighbor( 18, tBasis[  36 ] );
                 mBasis[  12 ]->insert_neighbor( 19, tBasis[  38 ] );
                 mBasis[  12 ]->insert_neighbor( 20, mBasis[  32 ] );
                 mBasis[  12 ]->insert_neighbor( 21, tBasis[  44 ] );
                 mBasis[  12 ]->insert_neighbor( 22, tBasis[  76 ] );
                 mBasis[  12 ]->insert_neighbor( 23, tBasis[  78 ] );
                 mBasis[  12 ]->insert_neighbor( 24, mBasis[  60 ] );
                 mBasis[  12 ]->insert_neighbor( 25, tBasis[  84 ] );

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
                 mBasis[  13 ]->insert_neighbor(  0, tBasis[  77 ] );
                 mBasis[  13 ]->insert_neighbor(  1, mBasis[  39 ] );
                 mBasis[  13 ]->insert_neighbor(  2, mBasis[  41 ] );
                 mBasis[  13 ]->insert_neighbor(  3, tBasis[  82 ] );
                 mBasis[  13 ]->insert_neighbor(  4, mBasis[  12 ] );
                 mBasis[  13 ]->insert_neighbor(  5, mBasis[   4 ] );
                 mBasis[  13 ]->insert_neighbor(  6, tBasis[  57 ] );
                 mBasis[  13 ]->insert_neighbor(  7, mBasis[  36 ] );
                 mBasis[  13 ]->insert_neighbor(  8, mBasis[  40 ] );
                 mBasis[  13 ]->insert_neighbor(  9, tBasis[  62 ] );
                 mBasis[  13 ]->insert_neighbor( 10, tBasis[  76 ] );
                 mBasis[  13 ]->insert_neighbor( 11, tBasis[  78 ] );
                 mBasis[  13 ]->insert_neighbor( 12, mBasis[  60 ] );
                 mBasis[  13 ]->insert_neighbor( 13, tBasis[  84 ] );
                 mBasis[  13 ]->insert_neighbor( 14, tBasis[  97 ] );
                 mBasis[  13 ]->insert_neighbor( 15, mBasis[  24 ] );
                 mBasis[  13 ]->insert_neighbor( 16, mBasis[  26 ] );
                 mBasis[  13 ]->insert_neighbor( 17, tBasis[ 102 ] );
                 mBasis[  13 ]->insert_neighbor( 18, tBasis[  56 ] );
                 mBasis[  13 ]->insert_neighbor( 19, tBasis[  58 ] );
                 mBasis[  13 ]->insert_neighbor( 20, mBasis[  56 ] );
                 mBasis[  13 ]->insert_neighbor( 21, tBasis[  64 ] );
                 mBasis[  13 ]->insert_neighbor( 22, tBasis[  96 ] );
                 mBasis[  13 ]->insert_neighbor( 23, tBasis[  98 ] );
                 mBasis[  13 ]->insert_neighbor( 24, mBasis[  52 ] );
                 mBasis[  13 ]->insert_neighbor( 25, tBasis[ 104 ] );

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
                 mBasis[  14 ]->insert_neighbor(  0, mBasis[   1 ] );
                 mBasis[  14 ]->insert_neighbor(  1, tBasis[  45 ] );
                 mBasis[  14 ]->insert_neighbor(  2, mBasis[  15 ] );
                 mBasis[  14 ]->insert_neighbor(  3, mBasis[  35 ] );
                 mBasis[  14 ]->insert_neighbor(  4, tBasis[  16 ] );
                 mBasis[  14 ]->insert_neighbor(  5, mBasis[  44 ] );
                 mBasis[  14 ]->insert_neighbor(  6, tBasis[  10 ] );
                 mBasis[  14 ]->insert_neighbor(  7, tBasis[  17 ] );
                 mBasis[  14 ]->insert_neighbor(  8, tBasis[  22 ] );
                 mBasis[  14 ]->insert_neighbor(  9, tBasis[  15 ] );
                 mBasis[  14 ]->insert_neighbor( 10, mBasis[   9 ] );
                 mBasis[  14 ]->insert_neighbor( 11, tBasis[  43 ] );
                 mBasis[  14 ]->insert_neighbor( 12, tBasis[  47 ] );
                 mBasis[  14 ]->insert_neighbor( 13, mBasis[  34 ] );
                 mBasis[  14 ]->insert_neighbor( 14, mBasis[  16 ] );
                 mBasis[  14 ]->insert_neighbor( 15, tBasis[  65 ] );
                 mBasis[  14 ]->insert_neighbor( 16, mBasis[  45 ] );
                 mBasis[  14 ]->insert_neighbor( 17, mBasis[  57 ] );
                 mBasis[  14 ]->insert_neighbor( 18, tBasis[   9 ] );
                 mBasis[  14 ]->insert_neighbor( 19, tBasis[  11 ] );
                 mBasis[  14 ]->insert_neighbor( 20, tBasis[  23 ] );
                 mBasis[  14 ]->insert_neighbor( 21, tBasis[  21 ] );
                 mBasis[  14 ]->insert_neighbor( 22, mBasis[  37 ] );
                 mBasis[  14 ]->insert_neighbor( 23, tBasis[  63 ] );
                 mBasis[  14 ]->insert_neighbor( 24, tBasis[  67 ] );
                 mBasis[  14 ]->insert_neighbor( 25, mBasis[  58 ] );

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
                 mBasis[  15 ]->insert_neighbor(  0, mBasis[  14 ] );
                 mBasis[  15 ]->insert_neighbor(  1, tBasis[  47 ] );
                 mBasis[  15 ]->insert_neighbor(  2, mBasis[   2 ] );
                 mBasis[  15 ]->insert_neighbor(  3, mBasis[  34 ] );
                 mBasis[  15 ]->insert_neighbor(  4, tBasis[  22 ] );
                 mBasis[  15 ]->insert_neighbor(  5, mBasis[  45 ] );
                 mBasis[  15 ]->insert_neighbor(  6, tBasis[  16 ] );
                 mBasis[  15 ]->insert_neighbor(  7, tBasis[  23 ] );
                 mBasis[  15 ]->insert_neighbor(  8, tBasis[  28 ] );
                 mBasis[  15 ]->insert_neighbor(  9, tBasis[  21 ] );
                 mBasis[  15 ]->insert_neighbor( 10, mBasis[  35 ] );
                 mBasis[  15 ]->insert_neighbor( 11, tBasis[  45 ] );
                 mBasis[  15 ]->insert_neighbor( 12, tBasis[  49 ] );
                 mBasis[  15 ]->insert_neighbor( 13, mBasis[  18 ] );
                 mBasis[  15 ]->insert_neighbor( 14, mBasis[  44 ] );
                 mBasis[  15 ]->insert_neighbor( 15, tBasis[  67 ] );
                 mBasis[  15 ]->insert_neighbor( 16, mBasis[  20 ] );
                 mBasis[  15 ]->insert_neighbor( 17, mBasis[  58 ] );
                 mBasis[  15 ]->insert_neighbor( 18, tBasis[  15 ] );
                 mBasis[  15 ]->insert_neighbor( 19, tBasis[  17 ] );
                 mBasis[  15 ]->insert_neighbor( 20, tBasis[  29 ] );
                 mBasis[  15 ]->insert_neighbor( 21, tBasis[  27 ] );
                 mBasis[  15 ]->insert_neighbor( 22, mBasis[  57 ] );
                 mBasis[  15 ]->insert_neighbor( 23, tBasis[  65 ] );
                 mBasis[  15 ]->insert_neighbor( 24, tBasis[  69 ] );
                 mBasis[  15 ]->insert_neighbor( 25, mBasis[  48 ] );

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
                 mBasis[  16 ]->insert_neighbor(  0, tBasis[  60 ] );
                 mBasis[  16 ]->insert_neighbor(  1, tBasis[  63 ] );
                 mBasis[  16 ]->insert_neighbor(  2, mBasis[  44 ] );
                 mBasis[  16 ]->insert_neighbor(  3, mBasis[  37 ] );
                 mBasis[  16 ]->insert_neighbor(  4, mBasis[   1 ] );
                 mBasis[  16 ]->insert_neighbor(  5, mBasis[  17 ] );
                 mBasis[  16 ]->insert_neighbor(  6, tBasis[  40 ] );
                 mBasis[  16 ]->insert_neighbor(  7, tBasis[  43 ] );
                 mBasis[  16 ]->insert_neighbor(  8, mBasis[  14 ] );
                 mBasis[  16 ]->insert_neighbor(  9, mBasis[   9 ] );
                 mBasis[  16 ]->insert_neighbor( 10, tBasis[  59 ] );
                 mBasis[  16 ]->insert_neighbor( 11, tBasis[  61 ] );
                 mBasis[  16 ]->insert_neighbor( 12, tBasis[  65 ] );
                 mBasis[  16 ]->insert_neighbor( 13, mBasis[  57 ] );
                 mBasis[  16 ]->insert_neighbor( 14, tBasis[  80 ] );
                 mBasis[  16 ]->insert_neighbor( 15, tBasis[  83 ] );
                 mBasis[  16 ]->insert_neighbor( 16, mBasis[  47 ] );
                 mBasis[  16 ]->insert_neighbor( 17, mBasis[  38 ] );
                 mBasis[  16 ]->insert_neighbor( 18, tBasis[  39 ] );
                 mBasis[  16 ]->insert_neighbor( 19, tBasis[  41 ] );
                 mBasis[  16 ]->insert_neighbor( 20, tBasis[  45 ] );
                 mBasis[  16 ]->insert_neighbor( 21, mBasis[  35 ] );
                 mBasis[  16 ]->insert_neighbor( 22, tBasis[  79 ] );
                 mBasis[  16 ]->insert_neighbor( 23, tBasis[  81 ] );
                 mBasis[  16 ]->insert_neighbor( 24, tBasis[  85 ] );
                 mBasis[  16 ]->insert_neighbor( 25, mBasis[  61 ] );

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
                 mBasis[  17 ]->insert_neighbor(  0, tBasis[  80 ] );
                 mBasis[  17 ]->insert_neighbor(  1, tBasis[  83 ] );
                 mBasis[  17 ]->insert_neighbor(  2, mBasis[  47 ] );
                 mBasis[  17 ]->insert_neighbor(  3, mBasis[  38 ] );
                 mBasis[  17 ]->insert_neighbor(  4, mBasis[  16 ] );
                 mBasis[  17 ]->insert_neighbor(  5, mBasis[   5 ] );
                 mBasis[  17 ]->insert_neighbor(  6, tBasis[  60 ] );
                 mBasis[  17 ]->insert_neighbor(  7, tBasis[  63 ] );
                 mBasis[  17 ]->insert_neighbor(  8, mBasis[  44 ] );
                 mBasis[  17 ]->insert_neighbor(  9, mBasis[  37 ] );
                 mBasis[  17 ]->insert_neighbor( 10, tBasis[  79 ] );
                 mBasis[  17 ]->insert_neighbor( 11, tBasis[  81 ] );
                 mBasis[  17 ]->insert_neighbor( 12, tBasis[  85 ] );
                 mBasis[  17 ]->insert_neighbor( 13, mBasis[  61 ] );
                 mBasis[  17 ]->insert_neighbor( 14, tBasis[ 100 ] );
                 mBasis[  17 ]->insert_neighbor( 15, tBasis[ 103 ] );
                 mBasis[  17 ]->insert_neighbor( 16, mBasis[  28 ] );
                 mBasis[  17 ]->insert_neighbor( 17, mBasis[  25 ] );
                 mBasis[  17 ]->insert_neighbor( 18, tBasis[  59 ] );
                 mBasis[  17 ]->insert_neighbor( 19, tBasis[  61 ] );
                 mBasis[  17 ]->insert_neighbor( 20, tBasis[  65 ] );
                 mBasis[  17 ]->insert_neighbor( 21, mBasis[  57 ] );
                 mBasis[  17 ]->insert_neighbor( 22, tBasis[  99 ] );
                 mBasis[  17 ]->insert_neighbor( 23, tBasis[ 101 ] );
                 mBasis[  17 ]->insert_neighbor( 24, tBasis[ 105 ] );
                 mBasis[  17 ]->insert_neighbor( 25, mBasis[  53 ] );

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
                 mBasis[  18 ]->insert_neighbor(  0, mBasis[  34 ] );
                 mBasis[  18 ]->insert_neighbor(  1, mBasis[   2 ] );
                 mBasis[  18 ]->insert_neighbor(  2, tBasis[  53 ] );
                 mBasis[  18 ]->insert_neighbor(  3, mBasis[  19 ] );
                 mBasis[  18 ]->insert_neighbor(  4, tBasis[  27 ] );
                 mBasis[  18 ]->insert_neighbor(  5, mBasis[  48 ] );
                 mBasis[  18 ]->insert_neighbor(  6, tBasis[  21 ] );
                 mBasis[  18 ]->insert_neighbor(  7, tBasis[  28 ] );
                 mBasis[  18 ]->insert_neighbor(  8, tBasis[  33 ] );
                 mBasis[  18 ]->insert_neighbor(  9, tBasis[  26 ] );
                 mBasis[  18 ]->insert_neighbor( 10, mBasis[  33 ] );
                 mBasis[  18 ]->insert_neighbor( 11, mBasis[  15 ] );
                 mBasis[  18 ]->insert_neighbor( 12, tBasis[  54 ] );
                 mBasis[  18 ]->insert_neighbor( 13, tBasis[  52 ] );
                 mBasis[  18 ]->insert_neighbor( 14, mBasis[  58 ] );
                 mBasis[  18 ]->insert_neighbor( 15, mBasis[  20 ] );
                 mBasis[  18 ]->insert_neighbor( 16, tBasis[  73 ] );
                 mBasis[  18 ]->insert_neighbor( 17, mBasis[  49 ] );
                 mBasis[  18 ]->insert_neighbor( 18, tBasis[  20 ] );
                 mBasis[  18 ]->insert_neighbor( 19, tBasis[  22 ] );
                 mBasis[  18 ]->insert_neighbor( 20, tBasis[  34 ] );
                 mBasis[  18 ]->insert_neighbor( 21, tBasis[  32 ] );
                 mBasis[  18 ]->insert_neighbor( 22, mBasis[  59 ] );
                 mBasis[  18 ]->insert_neighbor( 23, mBasis[  45 ] );
                 mBasis[  18 ]->insert_neighbor( 24, tBasis[  74 ] );
                 mBasis[  18 ]->insert_neighbor( 25, tBasis[  72 ] );

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
                 mBasis[  19 ]->insert_neighbor(  0, mBasis[  33 ] );
                 mBasis[  19 ]->insert_neighbor(  1, mBasis[  18 ] );
                 mBasis[  19 ]->insert_neighbor(  2, tBasis[  52 ] );
                 mBasis[  19 ]->insert_neighbor(  3, mBasis[   3 ] );
                 mBasis[  19 ]->insert_neighbor(  4, tBasis[  26 ] );
                 mBasis[  19 ]->insert_neighbor(  5, mBasis[  49 ] );
                 mBasis[  19 ]->insert_neighbor(  6, tBasis[  20 ] );
                 mBasis[  19 ]->insert_neighbor(  7, tBasis[  27 ] );
                 mBasis[  19 ]->insert_neighbor(  8, tBasis[  32 ] );
                 mBasis[  19 ]->insert_neighbor(  9, tBasis[  25 ] );
                 mBasis[  19 ]->insert_neighbor( 10, mBasis[  11 ] );
                 mBasis[  19 ]->insert_neighbor( 11, mBasis[  34 ] );
                 mBasis[  19 ]->insert_neighbor( 12, tBasis[  53 ] );
                 mBasis[  19 ]->insert_neighbor( 13, tBasis[  51 ] );
                 mBasis[  19 ]->insert_neighbor( 14, mBasis[  59 ] );
                 mBasis[  19 ]->insert_neighbor( 15, mBasis[  48 ] );
                 mBasis[  19 ]->insert_neighbor( 16, tBasis[  72 ] );
                 mBasis[  19 ]->insert_neighbor( 17, mBasis[  22 ] );
                 mBasis[  19 ]->insert_neighbor( 18, tBasis[  19 ] );
                 mBasis[  19 ]->insert_neighbor( 19, tBasis[  21 ] );
                 mBasis[  19 ]->insert_neighbor( 20, tBasis[  33 ] );
                 mBasis[  19 ]->insert_neighbor( 21, tBasis[  31 ] );
                 mBasis[  19 ]->insert_neighbor( 22, mBasis[  43 ] );
                 mBasis[  19 ]->insert_neighbor( 23, mBasis[  58 ] );
                 mBasis[  19 ]->insert_neighbor( 24, tBasis[  73 ] );
                 mBasis[  19 ]->insert_neighbor( 25, tBasis[  71 ] );

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
                 mBasis[  20 ]->insert_neighbor(  0, mBasis[  45 ] );
                 mBasis[  20 ]->insert_neighbor(  1, tBasis[  69 ] );
                 mBasis[  20 ]->insert_neighbor(  2, tBasis[  74 ] );
                 mBasis[  20 ]->insert_neighbor(  3, mBasis[  48 ] );
                 mBasis[  20 ]->insert_neighbor(  4, mBasis[   2 ] );
                 mBasis[  20 ]->insert_neighbor(  5, mBasis[  21 ] );
                 mBasis[  20 ]->insert_neighbor(  6, mBasis[  15 ] );
                 mBasis[  20 ]->insert_neighbor(  7, tBasis[  49 ] );
                 mBasis[  20 ]->insert_neighbor(  8, tBasis[  54 ] );
                 mBasis[  20 ]->insert_neighbor(  9, mBasis[  18 ] );
                 mBasis[  20 ]->insert_neighbor( 10, mBasis[  58 ] );
                 mBasis[  20 ]->insert_neighbor( 11, tBasis[  67 ] );
                 mBasis[  20 ]->insert_neighbor( 12, tBasis[  75 ] );
                 mBasis[  20 ]->insert_neighbor( 13, tBasis[  73 ] );
                 mBasis[  20 ]->insert_neighbor( 14, mBasis[  46 ] );
                 mBasis[  20 ]->insert_neighbor( 15, tBasis[  89 ] );
                 mBasis[  20 ]->insert_neighbor( 16, tBasis[  94 ] );
                 mBasis[  20 ]->insert_neighbor( 17, mBasis[  51 ] );
                 mBasis[  20 ]->insert_neighbor( 18, mBasis[  34 ] );
                 mBasis[  20 ]->insert_neighbor( 19, tBasis[  47 ] );
                 mBasis[  20 ]->insert_neighbor( 20, tBasis[  55 ] );
                 mBasis[  20 ]->insert_neighbor( 21, tBasis[  53 ] );
                 mBasis[  20 ]->insert_neighbor( 22, mBasis[  62 ] );
                 mBasis[  20 ]->insert_neighbor( 23, tBasis[  87 ] );
                 mBasis[  20 ]->insert_neighbor( 24, tBasis[  95 ] );
                 mBasis[  20 ]->insert_neighbor( 25, tBasis[  93 ] );

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
                 mBasis[  21 ]->insert_neighbor(  0, mBasis[  46 ] );
                 mBasis[  21 ]->insert_neighbor(  1, tBasis[  89 ] );
                 mBasis[  21 ]->insert_neighbor(  2, tBasis[  94 ] );
                 mBasis[  21 ]->insert_neighbor(  3, mBasis[  51 ] );
                 mBasis[  21 ]->insert_neighbor(  4, mBasis[  20 ] );
                 mBasis[  21 ]->insert_neighbor(  5, mBasis[   6 ] );
                 mBasis[  21 ]->insert_neighbor(  6, mBasis[  45 ] );
                 mBasis[  21 ]->insert_neighbor(  7, tBasis[  69 ] );
                 mBasis[  21 ]->insert_neighbor(  8, tBasis[  74 ] );
                 mBasis[  21 ]->insert_neighbor(  9, mBasis[  48 ] );
                 mBasis[  21 ]->insert_neighbor( 10, mBasis[  62 ] );
                 mBasis[  21 ]->insert_neighbor( 11, tBasis[  87 ] );
                 mBasis[  21 ]->insert_neighbor( 12, tBasis[  95 ] );
                 mBasis[  21 ]->insert_neighbor( 13, tBasis[  93 ] );
                 mBasis[  21 ]->insert_neighbor( 14, mBasis[  29 ] );
                 mBasis[  21 ]->insert_neighbor( 15, tBasis[ 109 ] );
                 mBasis[  21 ]->insert_neighbor( 16, tBasis[ 114 ] );
                 mBasis[  21 ]->insert_neighbor( 17, mBasis[  30 ] );
                 mBasis[  21 ]->insert_neighbor( 18, mBasis[  58 ] );
                 mBasis[  21 ]->insert_neighbor( 19, tBasis[  67 ] );
                 mBasis[  21 ]->insert_neighbor( 20, tBasis[  75 ] );
                 mBasis[  21 ]->insert_neighbor( 21, tBasis[  73 ] );
                 mBasis[  21 ]->insert_neighbor( 22, mBasis[  54 ] );
                 mBasis[  21 ]->insert_neighbor( 23, tBasis[ 107 ] );
                 mBasis[  21 ]->insert_neighbor( 24, tBasis[ 115 ] );
                 mBasis[  21 ]->insert_neighbor( 25, tBasis[ 113 ] );

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
                 mBasis[  22 ]->insert_neighbor(  0, mBasis[  43 ] );
                 mBasis[  22 ]->insert_neighbor(  1, mBasis[  49 ] );
                 mBasis[  22 ]->insert_neighbor(  2, tBasis[  71 ] );
                 mBasis[  22 ]->insert_neighbor(  3, tBasis[  68 ] );
                 mBasis[  22 ]->insert_neighbor(  4, mBasis[   3 ] );
                 mBasis[  22 ]->insert_neighbor(  5, mBasis[  23 ] );
                 mBasis[  22 ]->insert_neighbor(  6, mBasis[  11 ] );
                 mBasis[  22 ]->insert_neighbor(  7, mBasis[  19 ] );
                 mBasis[  22 ]->insert_neighbor(  8, tBasis[  51 ] );
                 mBasis[  22 ]->insert_neighbor(  9, tBasis[  48 ] );
                 mBasis[  22 ]->insert_neighbor( 10, tBasis[  66 ] );
                 mBasis[  22 ]->insert_neighbor( 11, mBasis[  59 ] );
                 mBasis[  22 ]->insert_neighbor( 12, tBasis[  72 ] );
                 mBasis[  22 ]->insert_neighbor( 13, tBasis[  70 ] );
                 mBasis[  22 ]->insert_neighbor( 14, mBasis[  42 ] );
                 mBasis[  22 ]->insert_neighbor( 15, mBasis[  50 ] );
                 mBasis[  22 ]->insert_neighbor( 16, tBasis[  91 ] );
                 mBasis[  22 ]->insert_neighbor( 17, tBasis[  88 ] );
                 mBasis[  22 ]->insert_neighbor( 18, tBasis[  46 ] );
                 mBasis[  22 ]->insert_neighbor( 19, mBasis[  33 ] );
                 mBasis[  22 ]->insert_neighbor( 20, tBasis[  52 ] );
                 mBasis[  22 ]->insert_neighbor( 21, tBasis[  50 ] );
                 mBasis[  22 ]->insert_neighbor( 22, tBasis[  86 ] );
                 mBasis[  22 ]->insert_neighbor( 23, mBasis[  63 ] );
                 mBasis[  22 ]->insert_neighbor( 24, tBasis[  92 ] );
                 mBasis[  22 ]->insert_neighbor( 25, tBasis[  90 ] );

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
                 mBasis[  23 ]->insert_neighbor(  0, mBasis[  42 ] );
                 mBasis[  23 ]->insert_neighbor(  1, mBasis[  50 ] );
                 mBasis[  23 ]->insert_neighbor(  2, tBasis[  91 ] );
                 mBasis[  23 ]->insert_neighbor(  3, tBasis[  88 ] );
                 mBasis[  23 ]->insert_neighbor(  4, mBasis[  22 ] );
                 mBasis[  23 ]->insert_neighbor(  5, mBasis[   7 ] );
                 mBasis[  23 ]->insert_neighbor(  6, mBasis[  43 ] );
                 mBasis[  23 ]->insert_neighbor(  7, mBasis[  49 ] );
                 mBasis[  23 ]->insert_neighbor(  8, tBasis[  71 ] );
                 mBasis[  23 ]->insert_neighbor(  9, tBasis[  68 ] );
                 mBasis[  23 ]->insert_neighbor( 10, tBasis[  86 ] );
                 mBasis[  23 ]->insert_neighbor( 11, mBasis[  63 ] );
                 mBasis[  23 ]->insert_neighbor( 12, tBasis[  92 ] );
                 mBasis[  23 ]->insert_neighbor( 13, tBasis[  90 ] );
                 mBasis[  23 ]->insert_neighbor( 14, mBasis[  27 ] );
                 mBasis[  23 ]->insert_neighbor( 15, mBasis[  31 ] );
                 mBasis[  23 ]->insert_neighbor( 16, tBasis[ 111 ] );
                 mBasis[  23 ]->insert_neighbor( 17, tBasis[ 108 ] );
                 mBasis[  23 ]->insert_neighbor( 18, tBasis[  66 ] );
                 mBasis[  23 ]->insert_neighbor( 19, mBasis[  59 ] );
                 mBasis[  23 ]->insert_neighbor( 20, tBasis[  72 ] );
                 mBasis[  23 ]->insert_neighbor( 21, tBasis[  70 ] );
                 mBasis[  23 ]->insert_neighbor( 22, tBasis[ 106 ] );
                 mBasis[  23 ]->insert_neighbor( 23, mBasis[  55 ] );
                 mBasis[  23 ]->insert_neighbor( 24, tBasis[ 112 ] );
                 mBasis[  23 ]->insert_neighbor( 25, tBasis[ 110 ] );

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
                 mBasis[  24 ]->insert_neighbor(  0, tBasis[  98 ] );
                 mBasis[  24 ]->insert_neighbor(  1, mBasis[  25 ] );
                 mBasis[  24 ]->insert_neighbor(  2, mBasis[  52 ] );
                 mBasis[  24 ]->insert_neighbor(  3, mBasis[   4 ] );
                 mBasis[  24 ]->insert_neighbor(  4, mBasis[  39 ] );
                 mBasis[  24 ]->insert_neighbor(  5, tBasis[ 124 ] );
                 mBasis[  24 ]->insert_neighbor(  6, tBasis[  78 ] );
                 mBasis[  24 ]->insert_neighbor(  7, mBasis[  38 ] );
                 mBasis[  24 ]->insert_neighbor(  8, mBasis[  60 ] );
                 mBasis[  24 ]->insert_neighbor(  9, mBasis[  13 ] );
                 mBasis[  24 ]->insert_neighbor( 10, tBasis[  97 ] );
                 mBasis[  24 ]->insert_neighbor( 11, tBasis[  99 ] );
                 mBasis[  24 ]->insert_neighbor( 12, mBasis[  53 ] );
                 mBasis[  24 ]->insert_neighbor( 13, mBasis[  26 ] );
                 mBasis[  24 ]->insert_neighbor( 14, tBasis[ 118 ] );
                 mBasis[  24 ]->insert_neighbor( 15, tBasis[ 125 ] );
                 mBasis[  24 ]->insert_neighbor( 16, tBasis[ 130 ] );
                 mBasis[  24 ]->insert_neighbor( 17, tBasis[ 123 ] );
                 mBasis[  24 ]->insert_neighbor( 18, tBasis[  77 ] );
                 mBasis[  24 ]->insert_neighbor( 19, tBasis[  79 ] );
                 mBasis[  24 ]->insert_neighbor( 20, mBasis[  61 ] );
                 mBasis[  24 ]->insert_neighbor( 21, mBasis[  41 ] );
                 mBasis[  24 ]->insert_neighbor( 22, tBasis[ 117 ] );
                 mBasis[  24 ]->insert_neighbor( 23, tBasis[ 119 ] );
                 mBasis[  24 ]->insert_neighbor( 24, tBasis[ 131 ] );
                 mBasis[  24 ]->insert_neighbor( 25, tBasis[ 129 ] );

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
                 mBasis[  25 ]->insert_neighbor(  0, tBasis[  99 ] );
                 mBasis[  25 ]->insert_neighbor(  1, mBasis[   5 ] );
                 mBasis[  25 ]->insert_neighbor(  2, mBasis[  53 ] );
                 mBasis[  25 ]->insert_neighbor(  3, mBasis[  24 ] );
                 mBasis[  25 ]->insert_neighbor(  4, mBasis[  38 ] );
                 mBasis[  25 ]->insert_neighbor(  5, tBasis[ 125 ] );
                 mBasis[  25 ]->insert_neighbor(  6, tBasis[  79 ] );
                 mBasis[  25 ]->insert_neighbor(  7, mBasis[  17 ] );
                 mBasis[  25 ]->insert_neighbor(  8, mBasis[  61 ] );
                 mBasis[  25 ]->insert_neighbor(  9, mBasis[  39 ] );
                 mBasis[  25 ]->insert_neighbor( 10, tBasis[  98 ] );
                 mBasis[  25 ]->insert_neighbor( 11, tBasis[ 100 ] );
                 mBasis[  25 ]->insert_neighbor( 12, mBasis[  28 ] );
                 mBasis[  25 ]->insert_neighbor( 13, mBasis[  52 ] );
                 mBasis[  25 ]->insert_neighbor( 14, tBasis[ 119 ] );
                 mBasis[  25 ]->insert_neighbor( 15, tBasis[ 126 ] );
                 mBasis[  25 ]->insert_neighbor( 16, tBasis[ 131 ] );
                 mBasis[  25 ]->insert_neighbor( 17, tBasis[ 124 ] );
                 mBasis[  25 ]->insert_neighbor( 18, tBasis[  78 ] );
                 mBasis[  25 ]->insert_neighbor( 19, tBasis[  80 ] );
                 mBasis[  25 ]->insert_neighbor( 20, mBasis[  47 ] );
                 mBasis[  25 ]->insert_neighbor( 21, mBasis[  60 ] );
                 mBasis[  25 ]->insert_neighbor( 22, tBasis[ 118 ] );
                 mBasis[  25 ]->insert_neighbor( 23, tBasis[ 120 ] );
                 mBasis[  25 ]->insert_neighbor( 24, tBasis[ 132 ] );
                 mBasis[  25 ]->insert_neighbor( 25, tBasis[ 130 ] );

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
                 mBasis[  26 ]->insert_neighbor(  0, mBasis[   4 ] );
                 mBasis[  26 ]->insert_neighbor(  1, mBasis[  52 ] );
                 mBasis[  26 ]->insert_neighbor(  2, mBasis[  27 ] );
                 mBasis[  26 ]->insert_neighbor(  3, tBasis[ 104 ] );
                 mBasis[  26 ]->insert_neighbor(  4, mBasis[  41 ] );
                 mBasis[  26 ]->insert_neighbor(  5, tBasis[ 129 ] );
                 mBasis[  26 ]->insert_neighbor(  6, mBasis[  13 ] );
                 mBasis[  26 ]->insert_neighbor(  7, mBasis[  60 ] );
                 mBasis[  26 ]->insert_neighbor(  8, mBasis[  42 ] );
                 mBasis[  26 ]->insert_neighbor(  9, tBasis[  84 ] );
                 mBasis[  26 ]->insert_neighbor( 10, tBasis[ 102 ] );
                 mBasis[  26 ]->insert_neighbor( 11, mBasis[  24 ] );
                 mBasis[  26 ]->insert_neighbor( 12, mBasis[  55 ] );
                 mBasis[  26 ]->insert_neighbor( 13, tBasis[ 106 ] );
                 mBasis[  26 ]->insert_neighbor( 14, tBasis[ 123 ] );
                 mBasis[  26 ]->insert_neighbor( 15, tBasis[ 130 ] );
                 mBasis[  26 ]->insert_neighbor( 16, tBasis[ 135 ] );
                 mBasis[  26 ]->insert_neighbor( 17, tBasis[ 128 ] );
                 mBasis[  26 ]->insert_neighbor( 18, tBasis[  82 ] );
                 mBasis[  26 ]->insert_neighbor( 19, mBasis[  39 ] );
                 mBasis[  26 ]->insert_neighbor( 20, mBasis[  63 ] );
                 mBasis[  26 ]->insert_neighbor( 21, tBasis[  86 ] );
                 mBasis[  26 ]->insert_neighbor( 22, tBasis[ 122 ] );
                 mBasis[  26 ]->insert_neighbor( 23, tBasis[ 124 ] );
                 mBasis[  26 ]->insert_neighbor( 24, tBasis[ 136 ] );
                 mBasis[  26 ]->insert_neighbor( 25, tBasis[ 134 ] );

                 // flag this basis
                 mBasis[  26 ]->flag();
             }

         }

         // test if basis 27 exists
         if ( mBasis[  27 ] != nullptr )
         {
             // test if basis 27 has been processed
             if ( ! mBasis[  27 ]->is_flagged() )
             {
                 // link neighbors of basis 27
                 mBasis[  27 ]->insert_neighbor(  0, mBasis[  26 ] );
                 mBasis[  27 ]->insert_neighbor(  1, mBasis[  55 ] );
                 mBasis[  27 ]->insert_neighbor(  2, mBasis[   7 ] );
                 mBasis[  27 ]->insert_neighbor(  3, tBasis[ 106 ] );
                 mBasis[  27 ]->insert_neighbor(  4, mBasis[  42 ] );
                 mBasis[  27 ]->insert_neighbor(  5, tBasis[ 135 ] );
                 mBasis[  27 ]->insert_neighbor(  6, mBasis[  41 ] );
                 mBasis[  27 ]->insert_neighbor(  7, mBasis[  63 ] );
                 mBasis[  27 ]->insert_neighbor(  8, mBasis[  23 ] );
                 mBasis[  27 ]->insert_neighbor(  9, tBasis[  86 ] );
                 mBasis[  27 ]->insert_neighbor( 10, tBasis[ 104 ] );
                 mBasis[  27 ]->insert_neighbor( 11, mBasis[  52 ] );
                 mBasis[  27 ]->insert_neighbor( 12, mBasis[  31 ] );
                 mBasis[  27 ]->insert_neighbor( 13, tBasis[ 108 ] );
                 mBasis[  27 ]->insert_neighbor( 14, tBasis[ 129 ] );
                 mBasis[  27 ]->insert_neighbor( 15, tBasis[ 136 ] );
                 mBasis[  27 ]->insert_neighbor( 16, tBasis[ 141 ] );
                 mBasis[  27 ]->insert_neighbor( 17, tBasis[ 134 ] );
                 mBasis[  27 ]->insert_neighbor( 18, tBasis[  84 ] );
                 mBasis[  27 ]->insert_neighbor( 19, mBasis[  60 ] );
                 mBasis[  27 ]->insert_neighbor( 20, mBasis[  50 ] );
                 mBasis[  27 ]->insert_neighbor( 21, tBasis[  88 ] );
                 mBasis[  27 ]->insert_neighbor( 22, tBasis[ 128 ] );
                 mBasis[  27 ]->insert_neighbor( 23, tBasis[ 130 ] );
                 mBasis[  27 ]->insert_neighbor( 24, tBasis[ 142 ] );
                 mBasis[  27 ]->insert_neighbor( 25, tBasis[ 140 ] );

                 // flag this basis
                 mBasis[  27 ]->flag();
             }

         }

         // test if basis 28 exists
         if ( mBasis[  28 ] != nullptr )
         {
             // test if basis 28 has been processed
             if ( ! mBasis[  28 ]->is_flagged() )
             {
                 // link neighbors of basis 28
                 mBasis[  28 ]->insert_neighbor(  0, mBasis[   5 ] );
                 mBasis[  28 ]->insert_neighbor(  1, tBasis[ 105 ] );
                 mBasis[  28 ]->insert_neighbor(  2, mBasis[  29 ] );
                 mBasis[  28 ]->insert_neighbor(  3, mBasis[  53 ] );
                 mBasis[  28 ]->insert_neighbor(  4, mBasis[  47 ] );
                 mBasis[  28 ]->insert_neighbor(  5, tBasis[ 132 ] );
                 mBasis[  28 ]->insert_neighbor(  6, mBasis[  17 ] );
                 mBasis[  28 ]->insert_neighbor(  7, tBasis[  85 ] );
                 mBasis[  28 ]->insert_neighbor(  8, mBasis[  46 ] );
                 mBasis[  28 ]->insert_neighbor(  9, mBasis[  61 ] );
                 mBasis[  28 ]->insert_neighbor( 10, mBasis[  25 ] );
                 mBasis[  28 ]->insert_neighbor( 11, tBasis[ 103 ] );
                 mBasis[  28 ]->insert_neighbor( 12, tBasis[ 107 ] );
                 mBasis[  28 ]->insert_neighbor( 13, mBasis[  54 ] );
                 mBasis[  28 ]->insert_neighbor( 14, tBasis[ 126 ] );
                 mBasis[  28 ]->insert_neighbor( 15, tBasis[ 133 ] );
                 mBasis[  28 ]->insert_neighbor( 16, tBasis[ 138 ] );
                 mBasis[  28 ]->insert_neighbor( 17, tBasis[ 131 ] );
                 mBasis[  28 ]->insert_neighbor( 18, mBasis[  38 ] );
                 mBasis[  28 ]->insert_neighbor( 19, tBasis[  83 ] );
                 mBasis[  28 ]->insert_neighbor( 20, tBasis[  87 ] );
                 mBasis[  28 ]->insert_neighbor( 21, mBasis[  62 ] );
                 mBasis[  28 ]->insert_neighbor( 22, tBasis[ 125 ] );
                 mBasis[  28 ]->insert_neighbor( 23, tBasis[ 127 ] );
                 mBasis[  28 ]->insert_neighbor( 24, tBasis[ 139 ] );
                 mBasis[  28 ]->insert_neighbor( 25, tBasis[ 137 ] );

                 // flag this basis
                 mBasis[  28 ]->flag();
             }

         }

         // test if basis 29 exists
         if ( mBasis[  29 ] != nullptr )
         {
             // test if basis 29 has been processed
             if ( ! mBasis[  29 ]->is_flagged() )
             {
                 // link neighbors of basis 29
                 mBasis[  29 ]->insert_neighbor(  0, mBasis[  28 ] );
                 mBasis[  29 ]->insert_neighbor(  1, tBasis[ 107 ] );
                 mBasis[  29 ]->insert_neighbor(  2, mBasis[   6 ] );
                 mBasis[  29 ]->insert_neighbor(  3, mBasis[  54 ] );
                 mBasis[  29 ]->insert_neighbor(  4, mBasis[  46 ] );
                 mBasis[  29 ]->insert_neighbor(  5, tBasis[ 138 ] );
                 mBasis[  29 ]->insert_neighbor(  6, mBasis[  47 ] );
                 mBasis[  29 ]->insert_neighbor(  7, tBasis[  87 ] );
                 mBasis[  29 ]->insert_neighbor(  8, mBasis[  21 ] );
                 mBasis[  29 ]->insert_neighbor(  9, mBasis[  62 ] );
                 mBasis[  29 ]->insert_neighbor( 10, mBasis[  53 ] );
                 mBasis[  29 ]->insert_neighbor( 11, tBasis[ 105 ] );
                 mBasis[  29 ]->insert_neighbor( 12, tBasis[ 109 ] );
                 mBasis[  29 ]->insert_neighbor( 13, mBasis[  30 ] );
                 mBasis[  29 ]->insert_neighbor( 14, tBasis[ 132 ] );
                 mBasis[  29 ]->insert_neighbor( 15, tBasis[ 139 ] );
                 mBasis[  29 ]->insert_neighbor( 16, tBasis[ 144 ] );
                 mBasis[  29 ]->insert_neighbor( 17, tBasis[ 137 ] );
                 mBasis[  29 ]->insert_neighbor( 18, mBasis[  61 ] );
                 mBasis[  29 ]->insert_neighbor( 19, tBasis[  85 ] );
                 mBasis[  29 ]->insert_neighbor( 20, tBasis[  89 ] );
                 mBasis[  29 ]->insert_neighbor( 21, mBasis[  51 ] );
                 mBasis[  29 ]->insert_neighbor( 22, tBasis[ 131 ] );
                 mBasis[  29 ]->insert_neighbor( 23, tBasis[ 133 ] );
                 mBasis[  29 ]->insert_neighbor( 24, tBasis[ 145 ] );
                 mBasis[  29 ]->insert_neighbor( 25, tBasis[ 143 ] );

                 // flag this basis
                 mBasis[  29 ]->flag();
             }

         }

         // test if basis 30 exists
         if ( mBasis[  30 ] != nullptr )
         {
             // test if basis 30 has been processed
             if ( ! mBasis[  30 ]->is_flagged() )
             {
                 // link neighbors of basis 30
                 mBasis[  30 ]->insert_neighbor(  0, mBasis[  54 ] );
                 mBasis[  30 ]->insert_neighbor(  1, mBasis[   6 ] );
                 mBasis[  30 ]->insert_neighbor(  2, tBasis[ 113 ] );
                 mBasis[  30 ]->insert_neighbor(  3, mBasis[  31 ] );
                 mBasis[  30 ]->insert_neighbor(  4, mBasis[  51 ] );
                 mBasis[  30 ]->insert_neighbor(  5, tBasis[ 143 ] );
                 mBasis[  30 ]->insert_neighbor(  6, mBasis[  62 ] );
                 mBasis[  30 ]->insert_neighbor(  7, mBasis[  21 ] );
                 mBasis[  30 ]->insert_neighbor(  8, tBasis[  93 ] );
                 mBasis[  30 ]->insert_neighbor(  9, mBasis[  50 ] );
                 mBasis[  30 ]->insert_neighbor( 10, mBasis[  55 ] );
                 mBasis[  30 ]->insert_neighbor( 11, mBasis[  29 ] );
                 mBasis[  30 ]->insert_neighbor( 12, tBasis[ 114 ] );
                 mBasis[  30 ]->insert_neighbor( 13, tBasis[ 112 ] );
                 mBasis[  30 ]->insert_neighbor( 14, tBasis[ 137 ] );
                 mBasis[  30 ]->insert_neighbor( 15, tBasis[ 144 ] );
                 mBasis[  30 ]->insert_neighbor( 16, tBasis[ 149 ] );
                 mBasis[  30 ]->insert_neighbor( 17, tBasis[ 142 ] );
                 mBasis[  30 ]->insert_neighbor( 18, mBasis[  63 ] );
                 mBasis[  30 ]->insert_neighbor( 19, mBasis[  46 ] );
                 mBasis[  30 ]->insert_neighbor( 20, tBasis[  94 ] );
                 mBasis[  30 ]->insert_neighbor( 21, tBasis[  92 ] );
                 mBasis[  30 ]->insert_neighbor( 22, tBasis[ 136 ] );
                 mBasis[  30 ]->insert_neighbor( 23, tBasis[ 138 ] );
                 mBasis[  30 ]->insert_neighbor( 24, tBasis[ 150 ] );
                 mBasis[  30 ]->insert_neighbor( 25, tBasis[ 148 ] );

                 // flag this basis
                 mBasis[  30 ]->flag();
             }

         }

         // test if basis 31 exists
         if ( mBasis[  31 ] != nullptr )
         {
             // test if basis 31 has been processed
             if ( ! mBasis[  31 ]->is_flagged() )
             {
                 // link neighbors of basis 31
                 mBasis[  31 ]->insert_neighbor(  0, mBasis[  55 ] );
                 mBasis[  31 ]->insert_neighbor(  1, mBasis[  30 ] );
                 mBasis[  31 ]->insert_neighbor(  2, tBasis[ 112 ] );
                 mBasis[  31 ]->insert_neighbor(  3, mBasis[   7 ] );
                 mBasis[  31 ]->insert_neighbor(  4, mBasis[  50 ] );
                 mBasis[  31 ]->insert_neighbor(  5, tBasis[ 142 ] );
                 mBasis[  31 ]->insert_neighbor(  6, mBasis[  63 ] );
                 mBasis[  31 ]->insert_neighbor(  7, mBasis[  51 ] );
                 mBasis[  31 ]->insert_neighbor(  8, tBasis[  92 ] );
                 mBasis[  31 ]->insert_neighbor(  9, mBasis[  23 ] );
                 mBasis[  31 ]->insert_neighbor( 10, mBasis[  27 ] );
                 mBasis[  31 ]->insert_neighbor( 11, mBasis[  54 ] );
                 mBasis[  31 ]->insert_neighbor( 12, tBasis[ 113 ] );
                 mBasis[  31 ]->insert_neighbor( 13, tBasis[ 111 ] );
                 mBasis[  31 ]->insert_neighbor( 14, tBasis[ 136 ] );
                 mBasis[  31 ]->insert_neighbor( 15, tBasis[ 143 ] );
                 mBasis[  31 ]->insert_neighbor( 16, tBasis[ 148 ] );
                 mBasis[  31 ]->insert_neighbor( 17, tBasis[ 141 ] );
                 mBasis[  31 ]->insert_neighbor( 18, mBasis[  42 ] );
                 mBasis[  31 ]->insert_neighbor( 19, mBasis[  62 ] );
                 mBasis[  31 ]->insert_neighbor( 20, tBasis[  93 ] );
                 mBasis[  31 ]->insert_neighbor( 21, tBasis[  91 ] );
                 mBasis[  31 ]->insert_neighbor( 22, tBasis[ 135 ] );
                 mBasis[  31 ]->insert_neighbor( 23, tBasis[ 137 ] );
                 mBasis[  31 ]->insert_neighbor( 24, tBasis[ 149 ] );
                 mBasis[  31 ]->insert_neighbor( 25, tBasis[ 147 ] );

                 // flag this basis
                 mBasis[  31 ]->flag();
             }

         }

         // test if basis 32 exists
         if ( mBasis[  32 ] != nullptr )
         {
             // test if basis 32 has been processed
             if ( ! mBasis[  32 ]->is_flagged() )
             {
                 // link neighbors of basis 32
                 mBasis[  32 ]->insert_neighbor(  0, mBasis[   8 ] );
                 mBasis[  32 ]->insert_neighbor(  1, mBasis[  35 ] );
                 mBasis[  32 ]->insert_neighbor(  2, mBasis[  33 ] );
                 mBasis[  32 ]->insert_neighbor(  3, mBasis[  10 ] );
                 mBasis[  32 ]->insert_neighbor(  4, tBasis[  14 ] );
                 mBasis[  32 ]->insert_neighbor(  5, mBasis[  56 ] );
                 mBasis[  32 ]->insert_neighbor(  6, tBasis[   8 ] );
                 mBasis[  32 ]->insert_neighbor(  7, tBasis[  15 ] );
                 mBasis[  32 ]->insert_neighbor(  8, tBasis[  20 ] );
                 mBasis[  32 ]->insert_neighbor(  9, tBasis[  13 ] );
                 mBasis[  32 ]->insert_neighbor( 10, mBasis[   0 ] );
                 mBasis[  32 ]->insert_neighbor( 11, mBasis[   9 ] );
                 mBasis[  32 ]->insert_neighbor( 12, mBasis[  34 ] );
                 mBasis[  32 ]->insert_neighbor( 13, mBasis[  11 ] );
                 mBasis[  32 ]->insert_neighbor( 14, mBasis[  36 ] );
                 mBasis[  32 ]->insert_neighbor( 15, mBasis[  57 ] );
                 mBasis[  32 ]->insert_neighbor( 16, mBasis[  59 ] );
                 mBasis[  32 ]->insert_neighbor( 17, mBasis[  40 ] );
                 mBasis[  32 ]->insert_neighbor( 18, tBasis[   7 ] );
                 mBasis[  32 ]->insert_neighbor( 19, tBasis[   9 ] );
                 mBasis[  32 ]->insert_neighbor( 20, tBasis[  21 ] );
                 mBasis[  32 ]->insert_neighbor( 21, tBasis[  19 ] );
                 mBasis[  32 ]->insert_neighbor( 22, mBasis[  12 ] );
                 mBasis[  32 ]->insert_neighbor( 23, mBasis[  37 ] );
                 mBasis[  32 ]->insert_neighbor( 24, mBasis[  58 ] );
                 mBasis[  32 ]->insert_neighbor( 25, mBasis[  43 ] );

                 // flag this basis
                 mBasis[  32 ]->flag();
             }

         }

         // test if basis 33 exists
         if ( mBasis[  33 ] != nullptr )
         {
             // test if basis 33 has been processed
             if ( ! mBasis[  33 ]->is_flagged() )
             {
                 // link neighbors of basis 33
                 mBasis[  33 ]->insert_neighbor(  0, mBasis[  32 ] );
                 mBasis[  33 ]->insert_neighbor(  1, mBasis[  34 ] );
                 mBasis[  33 ]->insert_neighbor(  2, mBasis[  19 ] );
                 mBasis[  33 ]->insert_neighbor(  3, mBasis[  11 ] );
                 mBasis[  33 ]->insert_neighbor(  4, tBasis[  20 ] );
                 mBasis[  33 ]->insert_neighbor(  5, mBasis[  59 ] );
                 mBasis[  33 ]->insert_neighbor(  6, tBasis[  14 ] );
                 mBasis[  33 ]->insert_neighbor(  7, tBasis[  21 ] );
                 mBasis[  33 ]->insert_neighbor(  8, tBasis[  26 ] );
                 mBasis[  33 ]->insert_neighbor(  9, tBasis[  19 ] );
                 mBasis[  33 ]->insert_neighbor( 10, mBasis[  10 ] );
                 mBasis[  33 ]->insert_neighbor( 11, mBasis[  35 ] );
                 mBasis[  33 ]->insert_neighbor( 12, mBasis[  18 ] );
                 mBasis[  33 ]->insert_neighbor( 13, mBasis[   3 ] );
                 mBasis[  33 ]->insert_neighbor( 14, mBasis[  56 ] );
                 mBasis[  33 ]->insert_neighbor( 15, mBasis[  58 ] );
                 mBasis[  33 ]->insert_neighbor( 16, mBasis[  49 ] );
                 mBasis[  33 ]->insert_neighbor( 17, mBasis[  43 ] );
                 mBasis[  33 ]->insert_neighbor( 18, tBasis[  13 ] );
                 mBasis[  33 ]->insert_neighbor( 19, tBasis[  15 ] );
                 mBasis[  33 ]->insert_neighbor( 20, tBasis[  27 ] );
                 mBasis[  33 ]->insert_neighbor( 21, tBasis[  25 ] );
                 mBasis[  33 ]->insert_neighbor( 22, mBasis[  40 ] );
                 mBasis[  33 ]->insert_neighbor( 23, mBasis[  57 ] );
                 mBasis[  33 ]->insert_neighbor( 24, mBasis[  48 ] );
                 mBasis[  33 ]->insert_neighbor( 25, mBasis[  22 ] );

                 // flag this basis
                 mBasis[  33 ]->flag();
             }

         }

         // test if basis 34 exists
         if ( mBasis[  34 ] != nullptr )
         {
             // test if basis 34 has been processed
             if ( ! mBasis[  34 ]->is_flagged() )
             {
                 // link neighbors of basis 34
                 mBasis[  34 ]->insert_neighbor(  0, mBasis[  35 ] );
                 mBasis[  34 ]->insert_neighbor(  1, mBasis[  15 ] );
                 mBasis[  34 ]->insert_neighbor(  2, mBasis[  18 ] );
                 mBasis[  34 ]->insert_neighbor(  3, mBasis[  33 ] );
                 mBasis[  34 ]->insert_neighbor(  4, tBasis[  21 ] );
                 mBasis[  34 ]->insert_neighbor(  5, mBasis[  58 ] );
                 mBasis[  34 ]->insert_neighbor(  6, tBasis[  15 ] );
                 mBasis[  34 ]->insert_neighbor(  7, tBasis[  22 ] );
                 mBasis[  34 ]->insert_neighbor(  8, tBasis[  27 ] );
                 mBasis[  34 ]->insert_neighbor(  9, tBasis[  20 ] );
                 mBasis[  34 ]->insert_neighbor( 10, mBasis[  32 ] );
                 mBasis[  34 ]->insert_neighbor( 11, mBasis[  14 ] );
                 mBasis[  34 ]->insert_neighbor( 12, mBasis[   2 ] );
                 mBasis[  34 ]->insert_neighbor( 13, mBasis[  19 ] );
                 mBasis[  34 ]->insert_neighbor( 14, mBasis[  57 ] );
                 mBasis[  34 ]->insert_neighbor( 15, mBasis[  45 ] );
                 mBasis[  34 ]->insert_neighbor( 16, mBasis[  48 ] );
                 mBasis[  34 ]->insert_neighbor( 17, mBasis[  59 ] );
                 mBasis[  34 ]->insert_neighbor( 18, tBasis[  14 ] );
                 mBasis[  34 ]->insert_neighbor( 19, tBasis[  16 ] );
                 mBasis[  34 ]->insert_neighbor( 20, tBasis[  28 ] );
                 mBasis[  34 ]->insert_neighbor( 21, tBasis[  26 ] );
                 mBasis[  34 ]->insert_neighbor( 22, mBasis[  56 ] );
                 mBasis[  34 ]->insert_neighbor( 23, mBasis[  44 ] );
                 mBasis[  34 ]->insert_neighbor( 24, mBasis[  20 ] );
                 mBasis[  34 ]->insert_neighbor( 25, mBasis[  49 ] );

                 // flag this basis
                 mBasis[  34 ]->flag();
             }

         }

         // test if basis 35 exists
         if ( mBasis[  35 ] != nullptr )
         {
             // test if basis 35 has been processed
             if ( ! mBasis[  35 ]->is_flagged() )
             {
                 // link neighbors of basis 35
                 mBasis[  35 ]->insert_neighbor(  0, mBasis[   9 ] );
                 mBasis[  35 ]->insert_neighbor(  1, mBasis[  14 ] );
                 mBasis[  35 ]->insert_neighbor(  2, mBasis[  34 ] );
                 mBasis[  35 ]->insert_neighbor(  3, mBasis[  32 ] );
                 mBasis[  35 ]->insert_neighbor(  4, tBasis[  15 ] );
                 mBasis[  35 ]->insert_neighbor(  5, mBasis[  57 ] );
                 mBasis[  35 ]->insert_neighbor(  6, tBasis[   9 ] );
                 mBasis[  35 ]->insert_neighbor(  7, tBasis[  16 ] );
                 mBasis[  35 ]->insert_neighbor(  8, tBasis[  21 ] );
                 mBasis[  35 ]->insert_neighbor(  9, tBasis[  14 ] );
                 mBasis[  35 ]->insert_neighbor( 10, mBasis[   8 ] );
                 mBasis[  35 ]->insert_neighbor( 11, mBasis[   1 ] );
                 mBasis[  35 ]->insert_neighbor( 12, mBasis[  15 ] );
                 mBasis[  35 ]->insert_neighbor( 13, mBasis[  33 ] );
                 mBasis[  35 ]->insert_neighbor( 14, mBasis[  37 ] );
                 mBasis[  35 ]->insert_neighbor( 15, mBasis[  44 ] );
                 mBasis[  35 ]->insert_neighbor( 16, mBasis[  58 ] );
                 mBasis[  35 ]->insert_neighbor( 17, mBasis[  56 ] );
                 mBasis[  35 ]->insert_neighbor( 18, tBasis[   8 ] );
                 mBasis[  35 ]->insert_neighbor( 19, tBasis[  10 ] );
                 mBasis[  35 ]->insert_neighbor( 20, tBasis[  22 ] );
                 mBasis[  35 ]->insert_neighbor( 21, tBasis[  20 ] );
                 mBasis[  35 ]->insert_neighbor( 22, mBasis[  36 ] );
                 mBasis[  35 ]->insert_neighbor( 23, mBasis[  16 ] );
                 mBasis[  35 ]->insert_neighbor( 24, mBasis[  45 ] );
                 mBasis[  35 ]->insert_neighbor( 25, mBasis[  59 ] );

                 // flag this basis
                 mBasis[  35 ]->flag();
             }

         }

         // test if basis 36 exists
         if ( mBasis[  36 ] != nullptr )
         {
             // test if basis 36 has been processed
             if ( ! mBasis[  36 ]->is_flagged() )
             {
                 // link neighbors of basis 36
                 mBasis[  36 ]->insert_neighbor(  0, tBasis[  58 ] );
                 mBasis[  36 ]->insert_neighbor(  1, mBasis[  37 ] );
                 mBasis[  36 ]->insert_neighbor(  2, mBasis[  56 ] );
                 mBasis[  36 ]->insert_neighbor(  3, mBasis[  12 ] );
                 mBasis[  36 ]->insert_neighbor(  4, mBasis[   8 ] );
                 mBasis[  36 ]->insert_neighbor(  5, mBasis[  39 ] );
                 mBasis[  36 ]->insert_neighbor(  6, tBasis[  38 ] );
                 mBasis[  36 ]->insert_neighbor(  7, mBasis[   9 ] );
                 mBasis[  36 ]->insert_neighbor(  8, mBasis[  32 ] );
                 mBasis[  36 ]->insert_neighbor(  9, mBasis[   0 ] );
                 mBasis[  36 ]->insert_neighbor( 10, tBasis[  57 ] );
                 mBasis[  36 ]->insert_neighbor( 11, tBasis[  59 ] );
                 mBasis[  36 ]->insert_neighbor( 12, mBasis[  57 ] );
                 mBasis[  36 ]->insert_neighbor( 13, mBasis[  40 ] );
                 mBasis[  36 ]->insert_neighbor( 14, tBasis[  78 ] );
                 mBasis[  36 ]->insert_neighbor( 15, mBasis[  38 ] );
                 mBasis[  36 ]->insert_neighbor( 16, mBasis[  60 ] );
                 mBasis[  36 ]->insert_neighbor( 17, mBasis[  13 ] );
                 mBasis[  36 ]->insert_neighbor( 18, tBasis[  37 ] );
                 mBasis[  36 ]->insert_neighbor( 19, tBasis[  39 ] );
                 mBasis[  36 ]->insert_neighbor( 20, mBasis[  35 ] );
                 mBasis[  36 ]->insert_neighbor( 21, mBasis[  10 ] );
                 mBasis[  36 ]->insert_neighbor( 22, tBasis[  77 ] );
                 mBasis[  36 ]->insert_neighbor( 23, tBasis[  79 ] );
                 mBasis[  36 ]->insert_neighbor( 24, mBasis[  61 ] );
                 mBasis[  36 ]->insert_neighbor( 25, mBasis[  41 ] );

                 // flag this basis
                 mBasis[  36 ]->flag();
             }

         }

         // test if basis 37 exists
         if ( mBasis[  37 ] != nullptr )
         {
             // test if basis 37 has been processed
             if ( ! mBasis[  37 ]->is_flagged() )
             {
                 // link neighbors of basis 37
                 mBasis[  37 ]->insert_neighbor(  0, tBasis[  59 ] );
                 mBasis[  37 ]->insert_neighbor(  1, mBasis[  16 ] );
                 mBasis[  37 ]->insert_neighbor(  2, mBasis[  57 ] );
                 mBasis[  37 ]->insert_neighbor(  3, mBasis[  36 ] );
                 mBasis[  37 ]->insert_neighbor(  4, mBasis[   9 ] );
                 mBasis[  37 ]->insert_neighbor(  5, mBasis[  38 ] );
                 mBasis[  37 ]->insert_neighbor(  6, tBasis[  39 ] );
                 mBasis[  37 ]->insert_neighbor(  7, mBasis[   1 ] );
                 mBasis[  37 ]->insert_neighbor(  8, mBasis[  35 ] );
                 mBasis[  37 ]->insert_neighbor(  9, mBasis[   8 ] );
                 mBasis[  37 ]->insert_neighbor( 10, tBasis[  58 ] );
                 mBasis[  37 ]->insert_neighbor( 11, tBasis[  60 ] );
                 mBasis[  37 ]->insert_neighbor( 12, mBasis[  44 ] );
                 mBasis[  37 ]->insert_neighbor( 13, mBasis[  56 ] );
                 mBasis[  37 ]->insert_neighbor( 14, tBasis[  79 ] );
                 mBasis[  37 ]->insert_neighbor( 15, mBasis[  17 ] );
                 mBasis[  37 ]->insert_neighbor( 16, mBasis[  61 ] );
                 mBasis[  37 ]->insert_neighbor( 17, mBasis[  39 ] );
                 mBasis[  37 ]->insert_neighbor( 18, tBasis[  38 ] );
                 mBasis[  37 ]->insert_neighbor( 19, tBasis[  40 ] );
                 mBasis[  37 ]->insert_neighbor( 20, mBasis[  14 ] );
                 mBasis[  37 ]->insert_neighbor( 21, mBasis[  32 ] );
                 mBasis[  37 ]->insert_neighbor( 22, tBasis[  78 ] );
                 mBasis[  37 ]->insert_neighbor( 23, tBasis[  80 ] );
                 mBasis[  37 ]->insert_neighbor( 24, mBasis[  47 ] );
                 mBasis[  37 ]->insert_neighbor( 25, mBasis[  60 ] );

                 // flag this basis
                 mBasis[  37 ]->flag();
             }

         }

         // test if basis 38 exists
         if ( mBasis[  38 ] != nullptr )
         {
             // test if basis 38 has been processed
             if ( ! mBasis[  38 ]->is_flagged() )
             {
                 // link neighbors of basis 38
                 mBasis[  38 ]->insert_neighbor(  0, tBasis[  79 ] );
                 mBasis[  38 ]->insert_neighbor(  1, mBasis[  17 ] );
                 mBasis[  38 ]->insert_neighbor(  2, mBasis[  61 ] );
                 mBasis[  38 ]->insert_neighbor(  3, mBasis[  39 ] );
                 mBasis[  38 ]->insert_neighbor(  4, mBasis[  37 ] );
                 mBasis[  38 ]->insert_neighbor(  5, mBasis[  25 ] );
                 mBasis[  38 ]->insert_neighbor(  6, tBasis[  59 ] );
                 mBasis[  38 ]->insert_neighbor(  7, mBasis[  16 ] );
                 mBasis[  38 ]->insert_neighbor(  8, mBasis[  57 ] );
                 mBasis[  38 ]->insert_neighbor(  9, mBasis[  36 ] );
                 mBasis[  38 ]->insert_neighbor( 10, tBasis[  78 ] );
                 mBasis[  38 ]->insert_neighbor( 11, tBasis[  80 ] );
                 mBasis[  38 ]->insert_neighbor( 12, mBasis[  47 ] );
                 mBasis[  38 ]->insert_neighbor( 13, mBasis[  60 ] );
                 mBasis[  38 ]->insert_neighbor( 14, tBasis[  99 ] );
                 mBasis[  38 ]->insert_neighbor( 15, mBasis[   5 ] );
                 mBasis[  38 ]->insert_neighbor( 16, mBasis[  53 ] );
                 mBasis[  38 ]->insert_neighbor( 17, mBasis[  24 ] );
                 mBasis[  38 ]->insert_neighbor( 18, tBasis[  58 ] );
                 mBasis[  38 ]->insert_neighbor( 19, tBasis[  60 ] );
                 mBasis[  38 ]->insert_neighbor( 20, mBasis[  44 ] );
                 mBasis[  38 ]->insert_neighbor( 21, mBasis[  56 ] );
                 mBasis[  38 ]->insert_neighbor( 22, tBasis[  98 ] );
                 mBasis[  38 ]->insert_neighbor( 23, tBasis[ 100 ] );
                 mBasis[  38 ]->insert_neighbor( 24, mBasis[  28 ] );
                 mBasis[  38 ]->insert_neighbor( 25, mBasis[  52 ] );

                 // flag this basis
                 mBasis[  38 ]->flag();
             }

         }

         // test if basis 39 exists
         if ( mBasis[  39 ] != nullptr )
         {
             // test if basis 39 has been processed
             if ( ! mBasis[  39 ]->is_flagged() )
             {
                 // link neighbors of basis 39
                 mBasis[  39 ]->insert_neighbor(  0, tBasis[  78 ] );
                 mBasis[  39 ]->insert_neighbor(  1, mBasis[  38 ] );
                 mBasis[  39 ]->insert_neighbor(  2, mBasis[  60 ] );
                 mBasis[  39 ]->insert_neighbor(  3, mBasis[  13 ] );
                 mBasis[  39 ]->insert_neighbor(  4, mBasis[  36 ] );
                 mBasis[  39 ]->insert_neighbor(  5, mBasis[  24 ] );
                 mBasis[  39 ]->insert_neighbor(  6, tBasis[  58 ] );
                 mBasis[  39 ]->insert_neighbor(  7, mBasis[  37 ] );
                 mBasis[  39 ]->insert_neighbor(  8, mBasis[  56 ] );
                 mBasis[  39 ]->insert_neighbor(  9, mBasis[  12 ] );
                 mBasis[  39 ]->insert_neighbor( 10, tBasis[  77 ] );
                 mBasis[  39 ]->insert_neighbor( 11, tBasis[  79 ] );
                 mBasis[  39 ]->insert_neighbor( 12, mBasis[  61 ] );
                 mBasis[  39 ]->insert_neighbor( 13, mBasis[  41 ] );
                 mBasis[  39 ]->insert_neighbor( 14, tBasis[  98 ] );
                 mBasis[  39 ]->insert_neighbor( 15, mBasis[  25 ] );
                 mBasis[  39 ]->insert_neighbor( 16, mBasis[  52 ] );
                 mBasis[  39 ]->insert_neighbor( 17, mBasis[   4 ] );
                 mBasis[  39 ]->insert_neighbor( 18, tBasis[  57 ] );
                 mBasis[  39 ]->insert_neighbor( 19, tBasis[  59 ] );
                 mBasis[  39 ]->insert_neighbor( 20, mBasis[  57 ] );
                 mBasis[  39 ]->insert_neighbor( 21, mBasis[  40 ] );
                 mBasis[  39 ]->insert_neighbor( 22, tBasis[  97 ] );
                 mBasis[  39 ]->insert_neighbor( 23, tBasis[  99 ] );
                 mBasis[  39 ]->insert_neighbor( 24, mBasis[  53 ] );
                 mBasis[  39 ]->insert_neighbor( 25, mBasis[  26 ] );

                 // flag this basis
                 mBasis[  39 ]->flag();
             }

         }

         // test if basis 40 exists
         if ( mBasis[  40 ] != nullptr )
         {
             // test if basis 40 has been processed
             if ( ! mBasis[  40 ]->is_flagged() )
             {
                 // link neighbors of basis 40
                 mBasis[  40 ]->insert_neighbor(  0, mBasis[  12 ] );
                 mBasis[  40 ]->insert_neighbor(  1, mBasis[  56 ] );
                 mBasis[  40 ]->insert_neighbor(  2, mBasis[  43 ] );
                 mBasis[  40 ]->insert_neighbor(  3, tBasis[  64 ] );
                 mBasis[  40 ]->insert_neighbor(  4, mBasis[  10 ] );
                 mBasis[  40 ]->insert_neighbor(  5, mBasis[  41 ] );
                 mBasis[  40 ]->insert_neighbor(  6, mBasis[   0 ] );
                 mBasis[  40 ]->insert_neighbor(  7, mBasis[  32 ] );
                 mBasis[  40 ]->insert_neighbor(  8, mBasis[  11 ] );
                 mBasis[  40 ]->insert_neighbor(  9, tBasis[  44 ] );
                 mBasis[  40 ]->insert_neighbor( 10, tBasis[  62 ] );
                 mBasis[  40 ]->insert_neighbor( 11, mBasis[  36 ] );
                 mBasis[  40 ]->insert_neighbor( 12, mBasis[  59 ] );
                 mBasis[  40 ]->insert_neighbor( 13, tBasis[  66 ] );
                 mBasis[  40 ]->insert_neighbor( 14, mBasis[  13 ] );
                 mBasis[  40 ]->insert_neighbor( 15, mBasis[  60 ] );
                 mBasis[  40 ]->insert_neighbor( 16, mBasis[  42 ] );
                 mBasis[  40 ]->insert_neighbor( 17, tBasis[  84 ] );
                 mBasis[  40 ]->insert_neighbor( 18, tBasis[  42 ] );
                 mBasis[  40 ]->insert_neighbor( 19, mBasis[   8 ] );
                 mBasis[  40 ]->insert_neighbor( 20, mBasis[  33 ] );
                 mBasis[  40 ]->insert_neighbor( 21, tBasis[  46 ] );
                 mBasis[  40 ]->insert_neighbor( 22, tBasis[  82 ] );
                 mBasis[  40 ]->insert_neighbor( 23, mBasis[  39 ] );
                 mBasis[  40 ]->insert_neighbor( 24, mBasis[  63 ] );
                 mBasis[  40 ]->insert_neighbor( 25, tBasis[  86 ] );

                 // flag this basis
                 mBasis[  40 ]->flag();
             }

         }

         // test if basis 41 exists
         if ( mBasis[  41 ] != nullptr )
         {
             // test if basis 41 has been processed
             if ( ! mBasis[  41 ]->is_flagged() )
             {
                 // link neighbors of basis 41
                 mBasis[  41 ]->insert_neighbor(  0, mBasis[  13 ] );
                 mBasis[  41 ]->insert_neighbor(  1, mBasis[  60 ] );
                 mBasis[  41 ]->insert_neighbor(  2, mBasis[  42 ] );
                 mBasis[  41 ]->insert_neighbor(  3, tBasis[  84 ] );
                 mBasis[  41 ]->insert_neighbor(  4, mBasis[  40 ] );
                 mBasis[  41 ]->insert_neighbor(  5, mBasis[  26 ] );
                 mBasis[  41 ]->insert_neighbor(  6, mBasis[  12 ] );
                 mBasis[  41 ]->insert_neighbor(  7, mBasis[  56 ] );
                 mBasis[  41 ]->insert_neighbor(  8, mBasis[  43 ] );
                 mBasis[  41 ]->insert_neighbor(  9, tBasis[  64 ] );
                 mBasis[  41 ]->insert_neighbor( 10, tBasis[  82 ] );
                 mBasis[  41 ]->insert_neighbor( 11, mBasis[  39 ] );
                 mBasis[  41 ]->insert_neighbor( 12, mBasis[  63 ] );
                 mBasis[  41 ]->insert_neighbor( 13, tBasis[  86 ] );
                 mBasis[  41 ]->insert_neighbor( 14, mBasis[   4 ] );
                 mBasis[  41 ]->insert_neighbor( 15, mBasis[  52 ] );
                 mBasis[  41 ]->insert_neighbor( 16, mBasis[  27 ] );
                 mBasis[  41 ]->insert_neighbor( 17, tBasis[ 104 ] );
                 mBasis[  41 ]->insert_neighbor( 18, tBasis[  62 ] );
                 mBasis[  41 ]->insert_neighbor( 19, mBasis[  36 ] );
                 mBasis[  41 ]->insert_neighbor( 20, mBasis[  59 ] );
                 mBasis[  41 ]->insert_neighbor( 21, tBasis[  66 ] );
                 mBasis[  41 ]->insert_neighbor( 22, tBasis[ 102 ] );
                 mBasis[  41 ]->insert_neighbor( 23, mBasis[  24 ] );
                 mBasis[  41 ]->insert_neighbor( 24, mBasis[  55 ] );
                 mBasis[  41 ]->insert_neighbor( 25, tBasis[ 106 ] );

                 // flag this basis
                 mBasis[  41 ]->flag();
             }

         }

         // test if basis 42 exists
         if ( mBasis[  42 ] != nullptr )
         {
             // test if basis 42 has been processed
             if ( ! mBasis[  42 ]->is_flagged() )
             {
                 // link neighbors of basis 42
                 mBasis[  42 ]->insert_neighbor(  0, mBasis[  41 ] );
                 mBasis[  42 ]->insert_neighbor(  1, mBasis[  63 ] );
                 mBasis[  42 ]->insert_neighbor(  2, mBasis[  23 ] );
                 mBasis[  42 ]->insert_neighbor(  3, tBasis[  86 ] );
                 mBasis[  42 ]->insert_neighbor(  4, mBasis[  43 ] );
                 mBasis[  42 ]->insert_neighbor(  5, mBasis[  27 ] );
                 mBasis[  42 ]->insert_neighbor(  6, mBasis[  40 ] );
                 mBasis[  42 ]->insert_neighbor(  7, mBasis[  59 ] );
                 mBasis[  42 ]->insert_neighbor(  8, mBasis[  22 ] );
                 mBasis[  42 ]->insert_neighbor(  9, tBasis[  66 ] );
                 mBasis[  42 ]->insert_neighbor( 10, tBasis[  84 ] );
                 mBasis[  42 ]->insert_neighbor( 11, mBasis[  60 ] );
                 mBasis[  42 ]->insert_neighbor( 12, mBasis[  50 ] );
                 mBasis[  42 ]->insert_neighbor( 13, tBasis[  88 ] );
                 mBasis[  42 ]->insert_neighbor( 14, mBasis[  26 ] );
                 mBasis[  42 ]->insert_neighbor( 15, mBasis[  55 ] );
                 mBasis[  42 ]->insert_neighbor( 16, mBasis[   7 ] );
                 mBasis[  42 ]->insert_neighbor( 17, tBasis[ 106 ] );
                 mBasis[  42 ]->insert_neighbor( 18, tBasis[  64 ] );
                 mBasis[  42 ]->insert_neighbor( 19, mBasis[  56 ] );
                 mBasis[  42 ]->insert_neighbor( 20, mBasis[  49 ] );
                 mBasis[  42 ]->insert_neighbor( 21, tBasis[  68 ] );
                 mBasis[  42 ]->insert_neighbor( 22, tBasis[ 104 ] );
                 mBasis[  42 ]->insert_neighbor( 23, mBasis[  52 ] );
                 mBasis[  42 ]->insert_neighbor( 24, mBasis[  31 ] );
                 mBasis[  42 ]->insert_neighbor( 25, tBasis[ 108 ] );

                 // flag this basis
                 mBasis[  42 ]->flag();
             }

         }

         // test if basis 43 exists
         if ( mBasis[  43 ] != nullptr )
         {
             // test if basis 43 has been processed
             if ( ! mBasis[  43 ]->is_flagged() )
             {
                 // link neighbors of basis 43
                 mBasis[  43 ]->insert_neighbor(  0, mBasis[  40 ] );
                 mBasis[  43 ]->insert_neighbor(  1, mBasis[  59 ] );
                 mBasis[  43 ]->insert_neighbor(  2, mBasis[  22 ] );
                 mBasis[  43 ]->insert_neighbor(  3, tBasis[  66 ] );
                 mBasis[  43 ]->insert_neighbor(  4, mBasis[  11 ] );
                 mBasis[  43 ]->insert_neighbor(  5, mBasis[  42 ] );
                 mBasis[  43 ]->insert_neighbor(  6, mBasis[  10 ] );
                 mBasis[  43 ]->insert_neighbor(  7, mBasis[  33 ] );
                 mBasis[  43 ]->insert_neighbor(  8, mBasis[   3 ] );
                 mBasis[  43 ]->insert_neighbor(  9, tBasis[  46 ] );
                 mBasis[  43 ]->insert_neighbor( 10, tBasis[  64 ] );
                 mBasis[  43 ]->insert_neighbor( 11, mBasis[  56 ] );
                 mBasis[  43 ]->insert_neighbor( 12, mBasis[  49 ] );
                 mBasis[  43 ]->insert_neighbor( 13, tBasis[  68 ] );
                 mBasis[  43 ]->insert_neighbor( 14, mBasis[  41 ] );
                 mBasis[  43 ]->insert_neighbor( 15, mBasis[  63 ] );
                 mBasis[  43 ]->insert_neighbor( 16, mBasis[  23 ] );
                 mBasis[  43 ]->insert_neighbor( 17, tBasis[  86 ] );
                 mBasis[  43 ]->insert_neighbor( 18, tBasis[  44 ] );
                 mBasis[  43 ]->insert_neighbor( 19, mBasis[  32 ] );
                 mBasis[  43 ]->insert_neighbor( 20, mBasis[  19 ] );
                 mBasis[  43 ]->insert_neighbor( 21, tBasis[  48 ] );
                 mBasis[  43 ]->insert_neighbor( 22, tBasis[  84 ] );
                 mBasis[  43 ]->insert_neighbor( 23, mBasis[  60 ] );
                 mBasis[  43 ]->insert_neighbor( 24, mBasis[  50 ] );
                 mBasis[  43 ]->insert_neighbor( 25, tBasis[  88 ] );

                 // flag this basis
                 mBasis[  43 ]->flag();
             }

         }

         // test if basis 44 exists
         if ( mBasis[  44 ] != nullptr )
         {
             // test if basis 44 has been processed
             if ( ! mBasis[  44 ]->is_flagged() )
             {
                 // link neighbors of basis 44
                 mBasis[  44 ]->insert_neighbor(  0, mBasis[  16 ] );
                 mBasis[  44 ]->insert_neighbor(  1, tBasis[  65 ] );
                 mBasis[  44 ]->insert_neighbor(  2, mBasis[  45 ] );
                 mBasis[  44 ]->insert_neighbor(  3, mBasis[  57 ] );
                 mBasis[  44 ]->insert_neighbor(  4, mBasis[  14 ] );
                 mBasis[  44 ]->insert_neighbor(  5, mBasis[  47 ] );
                 mBasis[  44 ]->insert_neighbor(  6, mBasis[   1 ] );
                 mBasis[  44 ]->insert_neighbor(  7, tBasis[  45 ] );
                 mBasis[  44 ]->insert_neighbor(  8, mBasis[  15 ] );
                 mBasis[  44 ]->insert_neighbor(  9, mBasis[  35 ] );
                 mBasis[  44 ]->insert_neighbor( 10, mBasis[  37 ] );
                 mBasis[  44 ]->insert_neighbor( 11, tBasis[  63 ] );
                 mBasis[  44 ]->insert_neighbor( 12, tBasis[  67 ] );
                 mBasis[  44 ]->insert_neighbor( 13, mBasis[  58 ] );
                 mBasis[  44 ]->insert_neighbor( 14, mBasis[  17 ] );
                 mBasis[  44 ]->insert_neighbor( 15, tBasis[  85 ] );
                 mBasis[  44 ]->insert_neighbor( 16, mBasis[  46 ] );
                 mBasis[  44 ]->insert_neighbor( 17, mBasis[  61 ] );
                 mBasis[  44 ]->insert_neighbor( 18, mBasis[   9 ] );
                 mBasis[  44 ]->insert_neighbor( 19, tBasis[  43 ] );
                 mBasis[  44 ]->insert_neighbor( 20, tBasis[  47 ] );
                 mBasis[  44 ]->insert_neighbor( 21, mBasis[  34 ] );
                 mBasis[  44 ]->insert_neighbor( 22, mBasis[  38 ] );
                 mBasis[  44 ]->insert_neighbor( 23, tBasis[  83 ] );
                 mBasis[  44 ]->insert_neighbor( 24, tBasis[  87 ] );
                 mBasis[  44 ]->insert_neighbor( 25, mBasis[  62 ] );

                 // flag this basis
                 mBasis[  44 ]->flag();
             }

         }

         // test if basis 45 exists
         if ( mBasis[  45 ] != nullptr )
         {
             // test if basis 45 has been processed
             if ( ! mBasis[  45 ]->is_flagged() )
             {
                 // link neighbors of basis 45
                 mBasis[  45 ]->insert_neighbor(  0, mBasis[  44 ] );
                 mBasis[  45 ]->insert_neighbor(  1, tBasis[  67 ] );
                 mBasis[  45 ]->insert_neighbor(  2, mBasis[  20 ] );
                 mBasis[  45 ]->insert_neighbor(  3, mBasis[  58 ] );
                 mBasis[  45 ]->insert_neighbor(  4, mBasis[  15 ] );
                 mBasis[  45 ]->insert_neighbor(  5, mBasis[  46 ] );
                 mBasis[  45 ]->insert_neighbor(  6, mBasis[  14 ] );
                 mBasis[  45 ]->insert_neighbor(  7, tBasis[  47 ] );
                 mBasis[  45 ]->insert_neighbor(  8, mBasis[   2 ] );
                 mBasis[  45 ]->insert_neighbor(  9, mBasis[  34 ] );
                 mBasis[  45 ]->insert_neighbor( 10, mBasis[  57 ] );
                 mBasis[  45 ]->insert_neighbor( 11, tBasis[  65 ] );
                 mBasis[  45 ]->insert_neighbor( 12, tBasis[  69 ] );
                 mBasis[  45 ]->insert_neighbor( 13, mBasis[  48 ] );
                 mBasis[  45 ]->insert_neighbor( 14, mBasis[  47 ] );
                 mBasis[  45 ]->insert_neighbor( 15, tBasis[  87 ] );
                 mBasis[  45 ]->insert_neighbor( 16, mBasis[  21 ] );
                 mBasis[  45 ]->insert_neighbor( 17, mBasis[  62 ] );
                 mBasis[  45 ]->insert_neighbor( 18, mBasis[  35 ] );
                 mBasis[  45 ]->insert_neighbor( 19, tBasis[  45 ] );
                 mBasis[  45 ]->insert_neighbor( 20, tBasis[  49 ] );
                 mBasis[  45 ]->insert_neighbor( 21, mBasis[  18 ] );
                 mBasis[  45 ]->insert_neighbor( 22, mBasis[  61 ] );
                 mBasis[  45 ]->insert_neighbor( 23, tBasis[  85 ] );
                 mBasis[  45 ]->insert_neighbor( 24, tBasis[  89 ] );
                 mBasis[  45 ]->insert_neighbor( 25, mBasis[  51 ] );

                 // flag this basis
                 mBasis[  45 ]->flag();
             }

         }

         // test if basis 46 exists
         if ( mBasis[  46 ] != nullptr )
         {
             // test if basis 46 has been processed
             if ( ! mBasis[  46 ]->is_flagged() )
             {
                 // link neighbors of basis 46
                 mBasis[  46 ]->insert_neighbor(  0, mBasis[  47 ] );
                 mBasis[  46 ]->insert_neighbor(  1, tBasis[  87 ] );
                 mBasis[  46 ]->insert_neighbor(  2, mBasis[  21 ] );
                 mBasis[  46 ]->insert_neighbor(  3, mBasis[  62 ] );
                 mBasis[  46 ]->insert_neighbor(  4, mBasis[  45 ] );
                 mBasis[  46 ]->insert_neighbor(  5, mBasis[  29 ] );
                 mBasis[  46 ]->insert_neighbor(  6, mBasis[  44 ] );
                 mBasis[  46 ]->insert_neighbor(  7, tBasis[  67 ] );
                 mBasis[  46 ]->insert_neighbor(  8, mBasis[  20 ] );
                 mBasis[  46 ]->insert_neighbor(  9, mBasis[  58 ] );
                 mBasis[  46 ]->insert_neighbor( 10, mBasis[  61 ] );
                 mBasis[  46 ]->insert_neighbor( 11, tBasis[  85 ] );
                 mBasis[  46 ]->insert_neighbor( 12, tBasis[  89 ] );
                 mBasis[  46 ]->insert_neighbor( 13, mBasis[  51 ] );
                 mBasis[  46 ]->insert_neighbor( 14, mBasis[  28 ] );
                 mBasis[  46 ]->insert_neighbor( 15, tBasis[ 107 ] );
                 mBasis[  46 ]->insert_neighbor( 16, mBasis[   6 ] );
                 mBasis[  46 ]->insert_neighbor( 17, mBasis[  54 ] );
                 mBasis[  46 ]->insert_neighbor( 18, mBasis[  57 ] );
                 mBasis[  46 ]->insert_neighbor( 19, tBasis[  65 ] );
                 mBasis[  46 ]->insert_neighbor( 20, tBasis[  69 ] );
                 mBasis[  46 ]->insert_neighbor( 21, mBasis[  48 ] );
                 mBasis[  46 ]->insert_neighbor( 22, mBasis[  53 ] );
                 mBasis[  46 ]->insert_neighbor( 23, tBasis[ 105 ] );
                 mBasis[  46 ]->insert_neighbor( 24, tBasis[ 109 ] );
                 mBasis[  46 ]->insert_neighbor( 25, mBasis[  30 ] );

                 // flag this basis
                 mBasis[  46 ]->flag();
             }

         }

         // test if basis 47 exists
         if ( mBasis[  47 ] != nullptr )
         {
             // test if basis 47 has been processed
             if ( ! mBasis[  47 ]->is_flagged() )
             {
                 // link neighbors of basis 47
                 mBasis[  47 ]->insert_neighbor(  0, mBasis[  17 ] );
                 mBasis[  47 ]->insert_neighbor(  1, tBasis[  85 ] );
                 mBasis[  47 ]->insert_neighbor(  2, mBasis[  46 ] );
                 mBasis[  47 ]->insert_neighbor(  3, mBasis[  61 ] );
                 mBasis[  47 ]->insert_neighbor(  4, mBasis[  44 ] );
                 mBasis[  47 ]->insert_neighbor(  5, mBasis[  28 ] );
                 mBasis[  47 ]->insert_neighbor(  6, mBasis[  16 ] );
                 mBasis[  47 ]->insert_neighbor(  7, tBasis[  65 ] );
                 mBasis[  47 ]->insert_neighbor(  8, mBasis[  45 ] );
                 mBasis[  47 ]->insert_neighbor(  9, mBasis[  57 ] );
                 mBasis[  47 ]->insert_neighbor( 10, mBasis[  38 ] );
                 mBasis[  47 ]->insert_neighbor( 11, tBasis[  83 ] );
                 mBasis[  47 ]->insert_neighbor( 12, tBasis[  87 ] );
                 mBasis[  47 ]->insert_neighbor( 13, mBasis[  62 ] );
                 mBasis[  47 ]->insert_neighbor( 14, mBasis[   5 ] );
                 mBasis[  47 ]->insert_neighbor( 15, tBasis[ 105 ] );
                 mBasis[  47 ]->insert_neighbor( 16, mBasis[  29 ] );
                 mBasis[  47 ]->insert_neighbor( 17, mBasis[  53 ] );
                 mBasis[  47 ]->insert_neighbor( 18, mBasis[  37 ] );
                 mBasis[  47 ]->insert_neighbor( 19, tBasis[  63 ] );
                 mBasis[  47 ]->insert_neighbor( 20, tBasis[  67 ] );
                 mBasis[  47 ]->insert_neighbor( 21, mBasis[  58 ] );
                 mBasis[  47 ]->insert_neighbor( 22, mBasis[  25 ] );
                 mBasis[  47 ]->insert_neighbor( 23, tBasis[ 103 ] );
                 mBasis[  47 ]->insert_neighbor( 24, tBasis[ 107 ] );
                 mBasis[  47 ]->insert_neighbor( 25, mBasis[  54 ] );

                 // flag this basis
                 mBasis[  47 ]->flag();
             }

         }

         // test if basis 48 exists
         if ( mBasis[  48 ] != nullptr )
         {
             // test if basis 48 has been processed
             if ( ! mBasis[  48 ]->is_flagged() )
             {
                 // link neighbors of basis 48
                 mBasis[  48 ]->insert_neighbor(  0, mBasis[  58 ] );
                 mBasis[  48 ]->insert_neighbor(  1, mBasis[  20 ] );
                 mBasis[  48 ]->insert_neighbor(  2, tBasis[  73 ] );
                 mBasis[  48 ]->insert_neighbor(  3, mBasis[  49 ] );
                 mBasis[  48 ]->insert_neighbor(  4, mBasis[  18 ] );
                 mBasis[  48 ]->insert_neighbor(  5, mBasis[  51 ] );
                 mBasis[  48 ]->insert_neighbor(  6, mBasis[  34 ] );
                 mBasis[  48 ]->insert_neighbor(  7, mBasis[   2 ] );
                 mBasis[  48 ]->insert_neighbor(  8, tBasis[  53 ] );
                 mBasis[  48 ]->insert_neighbor(  9, mBasis[  19 ] );
                 mBasis[  48 ]->insert_neighbor( 10, mBasis[  59 ] );
                 mBasis[  48 ]->insert_neighbor( 11, mBasis[  45 ] );
                 mBasis[  48 ]->insert_neighbor( 12, tBasis[  74 ] );
                 mBasis[  48 ]->insert_neighbor( 13, tBasis[  72 ] );
                 mBasis[  48 ]->insert_neighbor( 14, mBasis[  62 ] );
                 mBasis[  48 ]->insert_neighbor( 15, mBasis[  21 ] );
                 mBasis[  48 ]->insert_neighbor( 16, tBasis[  93 ] );
                 mBasis[  48 ]->insert_neighbor( 17, mBasis[  50 ] );
                 mBasis[  48 ]->insert_neighbor( 18, mBasis[  33 ] );
                 mBasis[  48 ]->insert_neighbor( 19, mBasis[  15 ] );
                 mBasis[  48 ]->insert_neighbor( 20, tBasis[  54 ] );
                 mBasis[  48 ]->insert_neighbor( 21, tBasis[  52 ] );
                 mBasis[  48 ]->insert_neighbor( 22, mBasis[  63 ] );
                 mBasis[  48 ]->insert_neighbor( 23, mBasis[  46 ] );
                 mBasis[  48 ]->insert_neighbor( 24, tBasis[  94 ] );
                 mBasis[  48 ]->insert_neighbor( 25, tBasis[  92 ] );

                 // flag this basis
                 mBasis[  48 ]->flag();
             }

         }

         // test if basis 49 exists
         if ( mBasis[  49 ] != nullptr )
         {
             // test if basis 49 has been processed
             if ( ! mBasis[  49 ]->is_flagged() )
             {
                 // link neighbors of basis 49
                 mBasis[  49 ]->insert_neighbor(  0, mBasis[  59 ] );
                 mBasis[  49 ]->insert_neighbor(  1, mBasis[  48 ] );
                 mBasis[  49 ]->insert_neighbor(  2, tBasis[  72 ] );
                 mBasis[  49 ]->insert_neighbor(  3, mBasis[  22 ] );
                 mBasis[  49 ]->insert_neighbor(  4, mBasis[  19 ] );
                 mBasis[  49 ]->insert_neighbor(  5, mBasis[  50 ] );
                 mBasis[  49 ]->insert_neighbor(  6, mBasis[  33 ] );
                 mBasis[  49 ]->insert_neighbor(  7, mBasis[  18 ] );
                 mBasis[  49 ]->insert_neighbor(  8, tBasis[  52 ] );
                 mBasis[  49 ]->insert_neighbor(  9, mBasis[   3 ] );
                 mBasis[  49 ]->insert_neighbor( 10, mBasis[  43 ] );
                 mBasis[  49 ]->insert_neighbor( 11, mBasis[  58 ] );
                 mBasis[  49 ]->insert_neighbor( 12, tBasis[  73 ] );
                 mBasis[  49 ]->insert_neighbor( 13, tBasis[  71 ] );
                 mBasis[  49 ]->insert_neighbor( 14, mBasis[  63 ] );
                 mBasis[  49 ]->insert_neighbor( 15, mBasis[  51 ] );
                 mBasis[  49 ]->insert_neighbor( 16, tBasis[  92 ] );
                 mBasis[  49 ]->insert_neighbor( 17, mBasis[  23 ] );
                 mBasis[  49 ]->insert_neighbor( 18, mBasis[  11 ] );
                 mBasis[  49 ]->insert_neighbor( 19, mBasis[  34 ] );
                 mBasis[  49 ]->insert_neighbor( 20, tBasis[  53 ] );
                 mBasis[  49 ]->insert_neighbor( 21, tBasis[  51 ] );
                 mBasis[  49 ]->insert_neighbor( 22, mBasis[  42 ] );
                 mBasis[  49 ]->insert_neighbor( 23, mBasis[  62 ] );
                 mBasis[  49 ]->insert_neighbor( 24, tBasis[  93 ] );
                 mBasis[  49 ]->insert_neighbor( 25, tBasis[  91 ] );

                 // flag this basis
                 mBasis[  49 ]->flag();
             }

         }

         // test if basis 50 exists
         if ( mBasis[  50 ] != nullptr )
         {
             // test if basis 50 has been processed
             if ( ! mBasis[  50 ]->is_flagged() )
             {
                 // link neighbors of basis 50
                 mBasis[  50 ]->insert_neighbor(  0, mBasis[  63 ] );
                 mBasis[  50 ]->insert_neighbor(  1, mBasis[  51 ] );
                 mBasis[  50 ]->insert_neighbor(  2, tBasis[  92 ] );
                 mBasis[  50 ]->insert_neighbor(  3, mBasis[  23 ] );
                 mBasis[  50 ]->insert_neighbor(  4, mBasis[  49 ] );
                 mBasis[  50 ]->insert_neighbor(  5, mBasis[  31 ] );
                 mBasis[  50 ]->insert_neighbor(  6, mBasis[  59 ] );
                 mBasis[  50 ]->insert_neighbor(  7, mBasis[  48 ] );
                 mBasis[  50 ]->insert_neighbor(  8, tBasis[  72 ] );
                 mBasis[  50 ]->insert_neighbor(  9, mBasis[  22 ] );
                 mBasis[  50 ]->insert_neighbor( 10, mBasis[  42 ] );
                 mBasis[  50 ]->insert_neighbor( 11, mBasis[  62 ] );
                 mBasis[  50 ]->insert_neighbor( 12, tBasis[  93 ] );
                 mBasis[  50 ]->insert_neighbor( 13, tBasis[  91 ] );
                 mBasis[  50 ]->insert_neighbor( 14, mBasis[  55 ] );
                 mBasis[  50 ]->insert_neighbor( 15, mBasis[  30 ] );
                 mBasis[  50 ]->insert_neighbor( 16, tBasis[ 112 ] );
                 mBasis[  50 ]->insert_neighbor( 17, mBasis[   7 ] );
                 mBasis[  50 ]->insert_neighbor( 18, mBasis[  43 ] );
                 mBasis[  50 ]->insert_neighbor( 19, mBasis[  58 ] );
                 mBasis[  50 ]->insert_neighbor( 20, tBasis[  73 ] );
                 mBasis[  50 ]->insert_neighbor( 21, tBasis[  71 ] );
                 mBasis[  50 ]->insert_neighbor( 22, mBasis[  27 ] );
                 mBasis[  50 ]->insert_neighbor( 23, mBasis[  54 ] );
                 mBasis[  50 ]->insert_neighbor( 24, tBasis[ 113 ] );
                 mBasis[  50 ]->insert_neighbor( 25, tBasis[ 111 ] );

                 // flag this basis
                 mBasis[  50 ]->flag();
             }

         }

         // test if basis 51 exists
         if ( mBasis[  51 ] != nullptr )
         {
             // test if basis 51 has been processed
             if ( ! mBasis[  51 ]->is_flagged() )
             {
                 // link neighbors of basis 51
                 mBasis[  51 ]->insert_neighbor(  0, mBasis[  62 ] );
                 mBasis[  51 ]->insert_neighbor(  1, mBasis[  21 ] );
                 mBasis[  51 ]->insert_neighbor(  2, tBasis[  93 ] );
                 mBasis[  51 ]->insert_neighbor(  3, mBasis[  50 ] );
                 mBasis[  51 ]->insert_neighbor(  4, mBasis[  48 ] );
                 mBasis[  51 ]->insert_neighbor(  5, mBasis[  30 ] );
                 mBasis[  51 ]->insert_neighbor(  6, mBasis[  58 ] );
                 mBasis[  51 ]->insert_neighbor(  7, mBasis[  20 ] );
                 mBasis[  51 ]->insert_neighbor(  8, tBasis[  73 ] );
                 mBasis[  51 ]->insert_neighbor(  9, mBasis[  49 ] );
                 mBasis[  51 ]->insert_neighbor( 10, mBasis[  63 ] );
                 mBasis[  51 ]->insert_neighbor( 11, mBasis[  46 ] );
                 mBasis[  51 ]->insert_neighbor( 12, tBasis[  94 ] );
                 mBasis[  51 ]->insert_neighbor( 13, tBasis[  92 ] );
                 mBasis[  51 ]->insert_neighbor( 14, mBasis[  54 ] );
                 mBasis[  51 ]->insert_neighbor( 15, mBasis[   6 ] );
                 mBasis[  51 ]->insert_neighbor( 16, tBasis[ 113 ] );
                 mBasis[  51 ]->insert_neighbor( 17, mBasis[  31 ] );
                 mBasis[  51 ]->insert_neighbor( 18, mBasis[  59 ] );
                 mBasis[  51 ]->insert_neighbor( 19, mBasis[  45 ] );
                 mBasis[  51 ]->insert_neighbor( 20, tBasis[  74 ] );
                 mBasis[  51 ]->insert_neighbor( 21, tBasis[  72 ] );
                 mBasis[  51 ]->insert_neighbor( 22, mBasis[  55 ] );
                 mBasis[  51 ]->insert_neighbor( 23, mBasis[  29 ] );
                 mBasis[  51 ]->insert_neighbor( 24, tBasis[ 114 ] );
                 mBasis[  51 ]->insert_neighbor( 25, tBasis[ 112 ] );

                 // flag this basis
                 mBasis[  51 ]->flag();
             }

         }

         // test if basis 52 exists
         if ( mBasis[  52 ] != nullptr )
         {
             // test if basis 52 has been processed
             if ( ! mBasis[  52 ]->is_flagged() )
             {
                 // link neighbors of basis 52
                 mBasis[  52 ]->insert_neighbor(  0, mBasis[  24 ] );
                 mBasis[  52 ]->insert_neighbor(  1, mBasis[  53 ] );
                 mBasis[  52 ]->insert_neighbor(  2, mBasis[  55 ] );
                 mBasis[  52 ]->insert_neighbor(  3, mBasis[  26 ] );
                 mBasis[  52 ]->insert_neighbor(  4, mBasis[  60 ] );
                 mBasis[  52 ]->insert_neighbor(  5, tBasis[ 130 ] );
                 mBasis[  52 ]->insert_neighbor(  6, mBasis[  39 ] );
                 mBasis[  52 ]->insert_neighbor(  7, mBasis[  61 ] );
                 mBasis[  52 ]->insert_neighbor(  8, mBasis[  63 ] );
                 mBasis[  52 ]->insert_neighbor(  9, mBasis[  41 ] );
                 mBasis[  52 ]->insert_neighbor( 10, mBasis[   4 ] );
                 mBasis[  52 ]->insert_neighbor( 11, mBasis[  25 ] );
                 mBasis[  52 ]->insert_neighbor( 12, mBasis[  54 ] );
                 mBasis[  52 ]->insert_neighbor( 13, mBasis[  27 ] );
                 mBasis[  52 ]->insert_neighbor( 14, tBasis[ 124 ] );
                 mBasis[  52 ]->insert_neighbor( 15, tBasis[ 131 ] );
                 mBasis[  52 ]->insert_neighbor( 16, tBasis[ 136 ] );
                 mBasis[  52 ]->insert_neighbor( 17, tBasis[ 129 ] );
                 mBasis[  52 ]->insert_neighbor( 18, mBasis[  13 ] );
                 mBasis[  52 ]->insert_neighbor( 19, mBasis[  38 ] );
                 mBasis[  52 ]->insert_neighbor( 20, mBasis[  62 ] );
                 mBasis[  52 ]->insert_neighbor( 21, mBasis[  42 ] );
                 mBasis[  52 ]->insert_neighbor( 22, tBasis[ 123 ] );
                 mBasis[  52 ]->insert_neighbor( 23, tBasis[ 125 ] );
                 mBasis[  52 ]->insert_neighbor( 24, tBasis[ 137 ] );
                 mBasis[  52 ]->insert_neighbor( 25, tBasis[ 135 ] );

                 // flag this basis
                 mBasis[  52 ]->flag();
             }

         }

         // test if basis 53 exists
         if ( mBasis[  53 ] != nullptr )
         {
             // test if basis 53 has been processed
             if ( ! mBasis[  53 ]->is_flagged() )
             {
                 // link neighbors of basis 53
                 mBasis[  53 ]->insert_neighbor(  0, mBasis[  25 ] );
                 mBasis[  53 ]->insert_neighbor(  1, mBasis[  28 ] );
                 mBasis[  53 ]->insert_neighbor(  2, mBasis[  54 ] );
                 mBasis[  53 ]->insert_neighbor(  3, mBasis[  52 ] );
                 mBasis[  53 ]->insert_neighbor(  4, mBasis[  61 ] );
                 mBasis[  53 ]->insert_neighbor(  5, tBasis[ 131 ] );
                 mBasis[  53 ]->insert_neighbor(  6, mBasis[  38 ] );
                 mBasis[  53 ]->insert_neighbor(  7, mBasis[  47 ] );
                 mBasis[  53 ]->insert_neighbor(  8, mBasis[  62 ] );
                 mBasis[  53 ]->insert_neighbor(  9, mBasis[  60 ] );
                 mBasis[  53 ]->insert_neighbor( 10, mBasis[  24 ] );
                 mBasis[  53 ]->insert_neighbor( 11, mBasis[   5 ] );
                 mBasis[  53 ]->insert_neighbor( 12, mBasis[  29 ] );
                 mBasis[  53 ]->insert_neighbor( 13, mBasis[  55 ] );
                 mBasis[  53 ]->insert_neighbor( 14, tBasis[ 125 ] );
                 mBasis[  53 ]->insert_neighbor( 15, tBasis[ 132 ] );
                 mBasis[  53 ]->insert_neighbor( 16, tBasis[ 137 ] );
                 mBasis[  53 ]->insert_neighbor( 17, tBasis[ 130 ] );
                 mBasis[  53 ]->insert_neighbor( 18, mBasis[  39 ] );
                 mBasis[  53 ]->insert_neighbor( 19, mBasis[  17 ] );
                 mBasis[  53 ]->insert_neighbor( 20, mBasis[  46 ] );
                 mBasis[  53 ]->insert_neighbor( 21, mBasis[  63 ] );
                 mBasis[  53 ]->insert_neighbor( 22, tBasis[ 124 ] );
                 mBasis[  53 ]->insert_neighbor( 23, tBasis[ 126 ] );
                 mBasis[  53 ]->insert_neighbor( 24, tBasis[ 138 ] );
                 mBasis[  53 ]->insert_neighbor( 25, tBasis[ 136 ] );

                 // flag this basis
                 mBasis[  53 ]->flag();
             }

         }

         // test if basis 54 exists
         if ( mBasis[  54 ] != nullptr )
         {
             // test if basis 54 has been processed
             if ( ! mBasis[  54 ]->is_flagged() )
             {
                 // link neighbors of basis 54
                 mBasis[  54 ]->insert_neighbor(  0, mBasis[  53 ] );
                 mBasis[  54 ]->insert_neighbor(  1, mBasis[  29 ] );
                 mBasis[  54 ]->insert_neighbor(  2, mBasis[  30 ] );
                 mBasis[  54 ]->insert_neighbor(  3, mBasis[  55 ] );
                 mBasis[  54 ]->insert_neighbor(  4, mBasis[  62 ] );
                 mBasis[  54 ]->insert_neighbor(  5, tBasis[ 137 ] );
                 mBasis[  54 ]->insert_neighbor(  6, mBasis[  61 ] );
                 mBasis[  54 ]->insert_neighbor(  7, mBasis[  46 ] );
                 mBasis[  54 ]->insert_neighbor(  8, mBasis[  51 ] );
                 mBasis[  54 ]->insert_neighbor(  9, mBasis[  63 ] );
                 mBasis[  54 ]->insert_neighbor( 10, mBasis[  52 ] );
                 mBasis[  54 ]->insert_neighbor( 11, mBasis[  28 ] );
                 mBasis[  54 ]->insert_neighbor( 12, mBasis[   6 ] );
                 mBasis[  54 ]->insert_neighbor( 13, mBasis[  31 ] );
                 mBasis[  54 ]->insert_neighbor( 14, tBasis[ 131 ] );
                 mBasis[  54 ]->insert_neighbor( 15, tBasis[ 138 ] );
                 mBasis[  54 ]->insert_neighbor( 16, tBasis[ 143 ] );
                 mBasis[  54 ]->insert_neighbor( 17, tBasis[ 136 ] );
                 mBasis[  54 ]->insert_neighbor( 18, mBasis[  60 ] );
                 mBasis[  54 ]->insert_neighbor( 19, mBasis[  47 ] );
                 mBasis[  54 ]->insert_neighbor( 20, mBasis[  21 ] );
                 mBasis[  54 ]->insert_neighbor( 21, mBasis[  50 ] );
                 mBasis[  54 ]->insert_neighbor( 22, tBasis[ 130 ] );
                 mBasis[  54 ]->insert_neighbor( 23, tBasis[ 132 ] );
                 mBasis[  54 ]->insert_neighbor( 24, tBasis[ 144 ] );
                 mBasis[  54 ]->insert_neighbor( 25, tBasis[ 142 ] );

                 // flag this basis
                 mBasis[  54 ]->flag();
             }

         }

         // test if basis 55 exists
         if ( mBasis[  55 ] != nullptr )
         {
             // test if basis 55 has been processed
             if ( ! mBasis[  55 ]->is_flagged() )
             {
                 // link neighbors of basis 55
                 mBasis[  55 ]->insert_neighbor(  0, mBasis[  52 ] );
                 mBasis[  55 ]->insert_neighbor(  1, mBasis[  54 ] );
                 mBasis[  55 ]->insert_neighbor(  2, mBasis[  31 ] );
                 mBasis[  55 ]->insert_neighbor(  3, mBasis[  27 ] );
                 mBasis[  55 ]->insert_neighbor(  4, mBasis[  63 ] );
                 mBasis[  55 ]->insert_neighbor(  5, tBasis[ 136 ] );
                 mBasis[  55 ]->insert_neighbor(  6, mBasis[  60 ] );
                 mBasis[  55 ]->insert_neighbor(  7, mBasis[  62 ] );
                 mBasis[  55 ]->insert_neighbor(  8, mBasis[  50 ] );
                 mBasis[  55 ]->insert_neighbor(  9, mBasis[  42 ] );
                 mBasis[  55 ]->insert_neighbor( 10, mBasis[  26 ] );
                 mBasis[  55 ]->insert_neighbor( 11, mBasis[  53 ] );
                 mBasis[  55 ]->insert_neighbor( 12, mBasis[  30 ] );
                 mBasis[  55 ]->insert_neighbor( 13, mBasis[   7 ] );
                 mBasis[  55 ]->insert_neighbor( 14, tBasis[ 130 ] );
                 mBasis[  55 ]->insert_neighbor( 15, tBasis[ 137 ] );
                 mBasis[  55 ]->insert_neighbor( 16, tBasis[ 142 ] );
                 mBasis[  55 ]->insert_neighbor( 17, tBasis[ 135 ] );
                 mBasis[  55 ]->insert_neighbor( 18, mBasis[  41 ] );
                 mBasis[  55 ]->insert_neighbor( 19, mBasis[  61 ] );
                 mBasis[  55 ]->insert_neighbor( 20, mBasis[  51 ] );
                 mBasis[  55 ]->insert_neighbor( 21, mBasis[  23 ] );
                 mBasis[  55 ]->insert_neighbor( 22, tBasis[ 129 ] );
                 mBasis[  55 ]->insert_neighbor( 23, tBasis[ 131 ] );
                 mBasis[  55 ]->insert_neighbor( 24, tBasis[ 143 ] );
                 mBasis[  55 ]->insert_neighbor( 25, tBasis[ 141 ] );

                 // flag this basis
                 mBasis[  55 ]->flag();
             }

         }

         // test if basis 56 exists
         if ( mBasis[  56 ] != nullptr )
         {
             // test if basis 56 has been processed
             if ( ! mBasis[  56 ]->is_flagged() )
             {
                 // link neighbors of basis 56
                 mBasis[  56 ]->insert_neighbor(  0, mBasis[  36 ] );
                 mBasis[  56 ]->insert_neighbor(  1, mBasis[  57 ] );
                 mBasis[  56 ]->insert_neighbor(  2, mBasis[  59 ] );
                 mBasis[  56 ]->insert_neighbor(  3, mBasis[  40 ] );
                 mBasis[  56 ]->insert_neighbor(  4, mBasis[  32 ] );
                 mBasis[  56 ]->insert_neighbor(  5, mBasis[  60 ] );
                 mBasis[  56 ]->insert_neighbor(  6, mBasis[   8 ] );
                 mBasis[  56 ]->insert_neighbor(  7, mBasis[  35 ] );
                 mBasis[  56 ]->insert_neighbor(  8, mBasis[  33 ] );
                 mBasis[  56 ]->insert_neighbor(  9, mBasis[  10 ] );
                 mBasis[  56 ]->insert_neighbor( 10, mBasis[  12 ] );
                 mBasis[  56 ]->insert_neighbor( 11, mBasis[  37 ] );
                 mBasis[  56 ]->insert_neighbor( 12, mBasis[  58 ] );
                 mBasis[  56 ]->insert_neighbor( 13, mBasis[  43 ] );
                 mBasis[  56 ]->insert_neighbor( 14, mBasis[  39 ] );
                 mBasis[  56 ]->insert_neighbor( 15, mBasis[  61 ] );
                 mBasis[  56 ]->insert_neighbor( 16, mBasis[  63 ] );
                 mBasis[  56 ]->insert_neighbor( 17, mBasis[  41 ] );
                 mBasis[  56 ]->insert_neighbor( 18, mBasis[   0 ] );
                 mBasis[  56 ]->insert_neighbor( 19, mBasis[   9 ] );
                 mBasis[  56 ]->insert_neighbor( 20, mBasis[  34 ] );
                 mBasis[  56 ]->insert_neighbor( 21, mBasis[  11 ] );
                 mBasis[  56 ]->insert_neighbor( 22, mBasis[  13 ] );
                 mBasis[  56 ]->insert_neighbor( 23, mBasis[  38 ] );
                 mBasis[  56 ]->insert_neighbor( 24, mBasis[  62 ] );
                 mBasis[  56 ]->insert_neighbor( 25, mBasis[  42 ] );

                 // flag this basis
                 mBasis[  56 ]->flag();
             }

         }

         // test if basis 57 exists
         if ( mBasis[  57 ] != nullptr )
         {
             // test if basis 57 has been processed
             if ( ! mBasis[  57 ]->is_flagged() )
             {
                 // link neighbors of basis 57
                 mBasis[  57 ]->insert_neighbor(  0, mBasis[  37 ] );
                 mBasis[  57 ]->insert_neighbor(  1, mBasis[  44 ] );
                 mBasis[  57 ]->insert_neighbor(  2, mBasis[  58 ] );
                 mBasis[  57 ]->insert_neighbor(  3, mBasis[  56 ] );
                 mBasis[  57 ]->insert_neighbor(  4, mBasis[  35 ] );
                 mBasis[  57 ]->insert_neighbor(  5, mBasis[  61 ] );
                 mBasis[  57 ]->insert_neighbor(  6, mBasis[   9 ] );
                 mBasis[  57 ]->insert_neighbor(  7, mBasis[  14 ] );
                 mBasis[  57 ]->insert_neighbor(  8, mBasis[  34 ] );
                 mBasis[  57 ]->insert_neighbor(  9, mBasis[  32 ] );
                 mBasis[  57 ]->insert_neighbor( 10, mBasis[  36 ] );
                 mBasis[  57 ]->insert_neighbor( 11, mBasis[  16 ] );
                 mBasis[  57 ]->insert_neighbor( 12, mBasis[  45 ] );
                 mBasis[  57 ]->insert_neighbor( 13, mBasis[  59 ] );
                 mBasis[  57 ]->insert_neighbor( 14, mBasis[  38 ] );
                 mBasis[  57 ]->insert_neighbor( 15, mBasis[  47 ] );
                 mBasis[  57 ]->insert_neighbor( 16, mBasis[  62 ] );
                 mBasis[  57 ]->insert_neighbor( 17, mBasis[  60 ] );
                 mBasis[  57 ]->insert_neighbor( 18, mBasis[   8 ] );
                 mBasis[  57 ]->insert_neighbor( 19, mBasis[   1 ] );
                 mBasis[  57 ]->insert_neighbor( 20, mBasis[  15 ] );
                 mBasis[  57 ]->insert_neighbor( 21, mBasis[  33 ] );
                 mBasis[  57 ]->insert_neighbor( 22, mBasis[  39 ] );
                 mBasis[  57 ]->insert_neighbor( 23, mBasis[  17 ] );
                 mBasis[  57 ]->insert_neighbor( 24, mBasis[  46 ] );
                 mBasis[  57 ]->insert_neighbor( 25, mBasis[  63 ] );

                 // flag this basis
                 mBasis[  57 ]->flag();
             }

         }

         // test if basis 58 exists
         if ( mBasis[  58 ] != nullptr )
         {
             // test if basis 58 has been processed
             if ( ! mBasis[  58 ]->is_flagged() )
             {
                 // link neighbors of basis 58
                 mBasis[  58 ]->insert_neighbor(  0, mBasis[  57 ] );
                 mBasis[  58 ]->insert_neighbor(  1, mBasis[  45 ] );
                 mBasis[  58 ]->insert_neighbor(  2, mBasis[  48 ] );
                 mBasis[  58 ]->insert_neighbor(  3, mBasis[  59 ] );
                 mBasis[  58 ]->insert_neighbor(  4, mBasis[  34 ] );
                 mBasis[  58 ]->insert_neighbor(  5, mBasis[  62 ] );
                 mBasis[  58 ]->insert_neighbor(  6, mBasis[  35 ] );
                 mBasis[  58 ]->insert_neighbor(  7, mBasis[  15 ] );
                 mBasis[  58 ]->insert_neighbor(  8, mBasis[  18 ] );
                 mBasis[  58 ]->insert_neighbor(  9, mBasis[  33 ] );
                 mBasis[  58 ]->insert_neighbor( 10, mBasis[  56 ] );
                 mBasis[  58 ]->insert_neighbor( 11, mBasis[  44 ] );
                 mBasis[  58 ]->insert_neighbor( 12, mBasis[  20 ] );
                 mBasis[  58 ]->insert_neighbor( 13, mBasis[  49 ] );
                 mBasis[  58 ]->insert_neighbor( 14, mBasis[  61 ] );
                 mBasis[  58 ]->insert_neighbor( 15, mBasis[  46 ] );
                 mBasis[  58 ]->insert_neighbor( 16, mBasis[  51 ] );
                 mBasis[  58 ]->insert_neighbor( 17, mBasis[  63 ] );
                 mBasis[  58 ]->insert_neighbor( 18, mBasis[  32 ] );
                 mBasis[  58 ]->insert_neighbor( 19, mBasis[  14 ] );
                 mBasis[  58 ]->insert_neighbor( 20, mBasis[   2 ] );
                 mBasis[  58 ]->insert_neighbor( 21, mBasis[  19 ] );
                 mBasis[  58 ]->insert_neighbor( 22, mBasis[  60 ] );
                 mBasis[  58 ]->insert_neighbor( 23, mBasis[  47 ] );
                 mBasis[  58 ]->insert_neighbor( 24, mBasis[  21 ] );
                 mBasis[  58 ]->insert_neighbor( 25, mBasis[  50 ] );

                 // flag this basis
                 mBasis[  58 ]->flag();
             }

         }

         // test if basis 59 exists
         if ( mBasis[  59 ] != nullptr )
         {
             // test if basis 59 has been processed
             if ( ! mBasis[  59 ]->is_flagged() )
             {
                 // link neighbors of basis 59
                 mBasis[  59 ]->insert_neighbor(  0, mBasis[  56 ] );
                 mBasis[  59 ]->insert_neighbor(  1, mBasis[  58 ] );
                 mBasis[  59 ]->insert_neighbor(  2, mBasis[  49 ] );
                 mBasis[  59 ]->insert_neighbor(  3, mBasis[  43 ] );
                 mBasis[  59 ]->insert_neighbor(  4, mBasis[  33 ] );
                 mBasis[  59 ]->insert_neighbor(  5, mBasis[  63 ] );
                 mBasis[  59 ]->insert_neighbor(  6, mBasis[  32 ] );
                 mBasis[  59 ]->insert_neighbor(  7, mBasis[  34 ] );
                 mBasis[  59 ]->insert_neighbor(  8, mBasis[  19 ] );
                 mBasis[  59 ]->insert_neighbor(  9, mBasis[  11 ] );
                 mBasis[  59 ]->insert_neighbor( 10, mBasis[  40 ] );
                 mBasis[  59 ]->insert_neighbor( 11, mBasis[  57 ] );
                 mBasis[  59 ]->insert_neighbor( 12, mBasis[  48 ] );
                 mBasis[  59 ]->insert_neighbor( 13, mBasis[  22 ] );
                 mBasis[  59 ]->insert_neighbor( 14, mBasis[  60 ] );
                 mBasis[  59 ]->insert_neighbor( 15, mBasis[  62 ] );
                 mBasis[  59 ]->insert_neighbor( 16, mBasis[  50 ] );
                 mBasis[  59 ]->insert_neighbor( 17, mBasis[  42 ] );
                 mBasis[  59 ]->insert_neighbor( 18, mBasis[  10 ] );
                 mBasis[  59 ]->insert_neighbor( 19, mBasis[  35 ] );
                 mBasis[  59 ]->insert_neighbor( 20, mBasis[  18 ] );
                 mBasis[  59 ]->insert_neighbor( 21, mBasis[   3 ] );
                 mBasis[  59 ]->insert_neighbor( 22, mBasis[  41 ] );
                 mBasis[  59 ]->insert_neighbor( 23, mBasis[  61 ] );
                 mBasis[  59 ]->insert_neighbor( 24, mBasis[  51 ] );
                 mBasis[  59 ]->insert_neighbor( 25, mBasis[  23 ] );

                 // flag this basis
                 mBasis[  59 ]->flag();
             }

         }

         // test if basis 60 exists
         if ( mBasis[  60 ] != nullptr )
         {
             // test if basis 60 has been processed
             if ( ! mBasis[  60 ]->is_flagged() )
             {
                 // link neighbors of basis 60
                 mBasis[  60 ]->insert_neighbor(  0, mBasis[  39 ] );
                 mBasis[  60 ]->insert_neighbor(  1, mBasis[  61 ] );
                 mBasis[  60 ]->insert_neighbor(  2, mBasis[  63 ] );
                 mBasis[  60 ]->insert_neighbor(  3, mBasis[  41 ] );
                 mBasis[  60 ]->insert_neighbor(  4, mBasis[  56 ] );
                 mBasis[  60 ]->insert_neighbor(  5, mBasis[  52 ] );
                 mBasis[  60 ]->insert_neighbor(  6, mBasis[  36 ] );
                 mBasis[  60 ]->insert_neighbor(  7, mBasis[  57 ] );
                 mBasis[  60 ]->insert_neighbor(  8, mBasis[  59 ] );
                 mBasis[  60 ]->insert_neighbor(  9, mBasis[  40 ] );
                 mBasis[  60 ]->insert_neighbor( 10, mBasis[  13 ] );
                 mBasis[  60 ]->insert_neighbor( 11, mBasis[  38 ] );
                 mBasis[  60 ]->insert_neighbor( 12, mBasis[  62 ] );
                 mBasis[  60 ]->insert_neighbor( 13, mBasis[  42 ] );
                 mBasis[  60 ]->insert_neighbor( 14, mBasis[  24 ] );
                 mBasis[  60 ]->insert_neighbor( 15, mBasis[  53 ] );
                 mBasis[  60 ]->insert_neighbor( 16, mBasis[  55 ] );
                 mBasis[  60 ]->insert_neighbor( 17, mBasis[  26 ] );
                 mBasis[  60 ]->insert_neighbor( 18, mBasis[  12 ] );
                 mBasis[  60 ]->insert_neighbor( 19, mBasis[  37 ] );
                 mBasis[  60 ]->insert_neighbor( 20, mBasis[  58 ] );
                 mBasis[  60 ]->insert_neighbor( 21, mBasis[  43 ] );
                 mBasis[  60 ]->insert_neighbor( 22, mBasis[   4 ] );
                 mBasis[  60 ]->insert_neighbor( 23, mBasis[  25 ] );
                 mBasis[  60 ]->insert_neighbor( 24, mBasis[  54 ] );
                 mBasis[  60 ]->insert_neighbor( 25, mBasis[  27 ] );

                 // flag this basis
                 mBasis[  60 ]->flag();
             }

         }

         // test if basis 61 exists
         if ( mBasis[  61 ] != nullptr )
         {
             // test if basis 61 has been processed
             if ( ! mBasis[  61 ]->is_flagged() )
             {
                 // link neighbors of basis 61
                 mBasis[  61 ]->insert_neighbor(  0, mBasis[  38 ] );
                 mBasis[  61 ]->insert_neighbor(  1, mBasis[  47 ] );
                 mBasis[  61 ]->insert_neighbor(  2, mBasis[  62 ] );
                 mBasis[  61 ]->insert_neighbor(  3, mBasis[  60 ] );
                 mBasis[  61 ]->insert_neighbor(  4, mBasis[  57 ] );
                 mBasis[  61 ]->insert_neighbor(  5, mBasis[  53 ] );
                 mBasis[  61 ]->insert_neighbor(  6, mBasis[  37 ] );
                 mBasis[  61 ]->insert_neighbor(  7, mBasis[  44 ] );
                 mBasis[  61 ]->insert_neighbor(  8, mBasis[  58 ] );
                 mBasis[  61 ]->insert_neighbor(  9, mBasis[  56 ] );
                 mBasis[  61 ]->insert_neighbor( 10, mBasis[  39 ] );
                 mBasis[  61 ]->insert_neighbor( 11, mBasis[  17 ] );
                 mBasis[  61 ]->insert_neighbor( 12, mBasis[  46 ] );
                 mBasis[  61 ]->insert_neighbor( 13, mBasis[  63 ] );
                 mBasis[  61 ]->insert_neighbor( 14, mBasis[  25 ] );
                 mBasis[  61 ]->insert_neighbor( 15, mBasis[  28 ] );
                 mBasis[  61 ]->insert_neighbor( 16, mBasis[  54 ] );
                 mBasis[  61 ]->insert_neighbor( 17, mBasis[  52 ] );
                 mBasis[  61 ]->insert_neighbor( 18, mBasis[  36 ] );
                 mBasis[  61 ]->insert_neighbor( 19, mBasis[  16 ] );
                 mBasis[  61 ]->insert_neighbor( 20, mBasis[  45 ] );
                 mBasis[  61 ]->insert_neighbor( 21, mBasis[  59 ] );
                 mBasis[  61 ]->insert_neighbor( 22, mBasis[  24 ] );
                 mBasis[  61 ]->insert_neighbor( 23, mBasis[   5 ] );
                 mBasis[  61 ]->insert_neighbor( 24, mBasis[  29 ] );
                 mBasis[  61 ]->insert_neighbor( 25, mBasis[  55 ] );

                 // flag this basis
                 mBasis[  61 ]->flag();
             }

         }

         // test if basis 62 exists
         if ( mBasis[  62 ] != nullptr )
         {
             // test if basis 62 has been processed
             if ( ! mBasis[  62 ]->is_flagged() )
             {
                 // link neighbors of basis 62
                 mBasis[  62 ]->insert_neighbor(  0, mBasis[  61 ] );
                 mBasis[  62 ]->insert_neighbor(  1, mBasis[  46 ] );
                 mBasis[  62 ]->insert_neighbor(  2, mBasis[  51 ] );
                 mBasis[  62 ]->insert_neighbor(  3, mBasis[  63 ] );
                 mBasis[  62 ]->insert_neighbor(  4, mBasis[  58 ] );
                 mBasis[  62 ]->insert_neighbor(  5, mBasis[  54 ] );
                 mBasis[  62 ]->insert_neighbor(  6, mBasis[  57 ] );
                 mBasis[  62 ]->insert_neighbor(  7, mBasis[  45 ] );
                 mBasis[  62 ]->insert_neighbor(  8, mBasis[  48 ] );
                 mBasis[  62 ]->insert_neighbor(  9, mBasis[  59 ] );
                 mBasis[  62 ]->insert_neighbor( 10, mBasis[  60 ] );
                 mBasis[  62 ]->insert_neighbor( 11, mBasis[  47 ] );
                 mBasis[  62 ]->insert_neighbor( 12, mBasis[  21 ] );
                 mBasis[  62 ]->insert_neighbor( 13, mBasis[  50 ] );
                 mBasis[  62 ]->insert_neighbor( 14, mBasis[  53 ] );
                 mBasis[  62 ]->insert_neighbor( 15, mBasis[  29 ] );
                 mBasis[  62 ]->insert_neighbor( 16, mBasis[  30 ] );
                 mBasis[  62 ]->insert_neighbor( 17, mBasis[  55 ] );
                 mBasis[  62 ]->insert_neighbor( 18, mBasis[  56 ] );
                 mBasis[  62 ]->insert_neighbor( 19, mBasis[  44 ] );
                 mBasis[  62 ]->insert_neighbor( 20, mBasis[  20 ] );
                 mBasis[  62 ]->insert_neighbor( 21, mBasis[  49 ] );
                 mBasis[  62 ]->insert_neighbor( 22, mBasis[  52 ] );
                 mBasis[  62 ]->insert_neighbor( 23, mBasis[  28 ] );
                 mBasis[  62 ]->insert_neighbor( 24, mBasis[   6 ] );
                 mBasis[  62 ]->insert_neighbor( 25, mBasis[  31 ] );

                 // flag this basis
                 mBasis[  62 ]->flag();
             }

         }

         // test if basis 63 exists
         if ( mBasis[  63 ] != nullptr )
         {
             // test if basis 63 has been processed
             if ( ! mBasis[  63 ]->is_flagged() )
             {
                 // link neighbors of basis 63
                 mBasis[  63 ]->insert_neighbor(  0, mBasis[  60 ] );
                 mBasis[  63 ]->insert_neighbor(  1, mBasis[  62 ] );
                 mBasis[  63 ]->insert_neighbor(  2, mBasis[  50 ] );
                 mBasis[  63 ]->insert_neighbor(  3, mBasis[  42 ] );
                 mBasis[  63 ]->insert_neighbor(  4, mBasis[  59 ] );
                 mBasis[  63 ]->insert_neighbor(  5, mBasis[  55 ] );
                 mBasis[  63 ]->insert_neighbor(  6, mBasis[  56 ] );
                 mBasis[  63 ]->insert_neighbor(  7, mBasis[  58 ] );
                 mBasis[  63 ]->insert_neighbor(  8, mBasis[  49 ] );
                 mBasis[  63 ]->insert_neighbor(  9, mBasis[  43 ] );
                 mBasis[  63 ]->insert_neighbor( 10, mBasis[  41 ] );
                 mBasis[  63 ]->insert_neighbor( 11, mBasis[  61 ] );
                 mBasis[  63 ]->insert_neighbor( 12, mBasis[  51 ] );
                 mBasis[  63 ]->insert_neighbor( 13, mBasis[  23 ] );
                 mBasis[  63 ]->insert_neighbor( 14, mBasis[  52 ] );
                 mBasis[  63 ]->insert_neighbor( 15, mBasis[  54 ] );
                 mBasis[  63 ]->insert_neighbor( 16, mBasis[  31 ] );
                 mBasis[  63 ]->insert_neighbor( 17, mBasis[  27 ] );
                 mBasis[  63 ]->insert_neighbor( 18, mBasis[  40 ] );
                 mBasis[  63 ]->insert_neighbor( 19, mBasis[  57 ] );
                 mBasis[  63 ]->insert_neighbor( 20, mBasis[  48 ] );
                 mBasis[  63 ]->insert_neighbor( 21, mBasis[  22 ] );
                 mBasis[  63 ]->insert_neighbor( 22, mBasis[  26 ] );
                 mBasis[  63 ]->insert_neighbor( 23, mBasis[  53 ] );
                 mBasis[  63 ]->insert_neighbor( 24, mBasis[  30 ] );
                 mBasis[  63 ]->insert_neighbor( 25, mBasis[   7 ] );

                 // flag this basis
                 mBasis[  63 ]->flag();
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
    luint BSpline_Element< 3, 3, 3 >::refine_basis( uint aBasisNumber )
    {
        // Start basis counter
        luint tBasisCounter = 0;
         
        // get pointer to basis
        Basis_Function* tBasis = mBasis[ aBasisNumber ];

        // test if basis exists
        if ( tBasis != nullptr )
        {
            // test if basis has been refined already
            if ( ! tBasis->has_children() )
            {
                // initialize neighbor container
                Basis_Function* tNeighbors[ 124 ] = { nullptr };

                // populate neighbor container
                get_basis_neighbors_3d( tBasis, 2, tNeighbors );

                // initialize temporary child container
                Basis_Function* tChildren[ 125 ] = { nullptr };

                // test if neighbor 0 exists
                if( tNeighbors[ 0 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 0 ];

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
                       tChildren[  25 ] = tNeighbor->get_child(  35 );
                       tChildren[  26 ] = tNeighbor->get_child(  36 );
                       tChildren[  27 ] = tNeighbor->get_child(  37 );
                       tChildren[  28 ] = tNeighbor->get_child(  38 );
                       tChildren[  29 ] = tNeighbor->get_child(  39 );
                       tChildren[  30 ] = tNeighbor->get_child(  40 );
                       tChildren[  31 ] = tNeighbor->get_child(  41 );
                       tChildren[  32 ] = tNeighbor->get_child(  42 );
                       tChildren[  33 ] = tNeighbor->get_child(  43 );
                       tChildren[  34 ] = tNeighbor->get_child(  44 );
                       tChildren[  35 ] = tNeighbor->get_child(  45 );
                       tChildren[  36 ] = tNeighbor->get_child(  46 );
                       tChildren[  37 ] = tNeighbor->get_child(  47 );
                       tChildren[  38 ] = tNeighbor->get_child(  48 );
                       tChildren[  39 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  60 );
                       tChildren[  51 ] = tNeighbor->get_child(  61 );
                       tChildren[  52 ] = tNeighbor->get_child(  62 );
                       tChildren[  53 ] = tNeighbor->get_child(  63 );
                       tChildren[  54 ] = tNeighbor->get_child(  64 );
                       tChildren[  55 ] = tNeighbor->get_child(  65 );
                       tChildren[  56 ] = tNeighbor->get_child(  66 );
                       tChildren[  57 ] = tNeighbor->get_child(  67 );
                       tChildren[  58 ] = tNeighbor->get_child(  68 );
                       tChildren[  59 ] = tNeighbor->get_child(  69 );
                       tChildren[  60 ] = tNeighbor->get_child(  70 );
                       tChildren[  61 ] = tNeighbor->get_child(  71 );
                       tChildren[  62 ] = tNeighbor->get_child(  72 );
                       tChildren[  63 ] = tNeighbor->get_child(  73 );
                       tChildren[  64 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  85 );
                       tChildren[  76 ] = tNeighbor->get_child(  86 );
                       tChildren[  77 ] = tNeighbor->get_child(  87 );
                       tChildren[  78 ] = tNeighbor->get_child(  88 );
                       tChildren[  79 ] = tNeighbor->get_child(  89 );
                       tChildren[  80 ] = tNeighbor->get_child(  90 );
                       tChildren[  81 ] = tNeighbor->get_child(  91 );
                       tChildren[  82 ] = tNeighbor->get_child(  92 );
                       tChildren[  83 ] = tNeighbor->get_child(  93 );
                       tChildren[  84 ] = tNeighbor->get_child(  94 );
                       tChildren[  85 ] = tNeighbor->get_child(  95 );
                       tChildren[  86 ] = tNeighbor->get_child(  96 );
                       tChildren[  87 ] = tNeighbor->get_child(  97 );
                       tChildren[  88 ] = tNeighbor->get_child(  98 );
                       tChildren[  89 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 110 );
                       tChildren[ 101 ] = tNeighbor->get_child( 111 );
                       tChildren[ 102 ] = tNeighbor->get_child( 112 );
                       tChildren[ 103 ] = tNeighbor->get_child( 113 );
                       tChildren[ 104 ] = tNeighbor->get_child( 114 );
                       tChildren[ 105 ] = tNeighbor->get_child( 115 );
                       tChildren[ 106 ] = tNeighbor->get_child( 116 );
                       tChildren[ 107 ] = tNeighbor->get_child( 117 );
                       tChildren[ 108 ] = tNeighbor->get_child( 118 );
                       tChildren[ 109 ] = tNeighbor->get_child( 119 );
                       tChildren[ 110 ] = tNeighbor->get_child( 120 );
                       tChildren[ 111 ] = tNeighbor->get_child( 121 );
                       tChildren[ 112 ] = tNeighbor->get_child( 122 );
                       tChildren[ 113 ] = tNeighbor->get_child( 123 );
                       tChildren[ 114 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 1 exists
                if( tNeighbors[ 1 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 1 ];

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
                       tChildren[  27 ] = tNeighbor->get_child(  25 );
                       tChildren[  28 ] = tNeighbor->get_child(  26 );
                       tChildren[  29 ] = tNeighbor->get_child(  27 );
                       tChildren[  32 ] = tNeighbor->get_child(  30 );
                       tChildren[  33 ] = tNeighbor->get_child(  31 );
                       tChildren[  34 ] = tNeighbor->get_child(  32 );
                       tChildren[  37 ] = tNeighbor->get_child(  35 );
                       tChildren[  38 ] = tNeighbor->get_child(  36 );
                       tChildren[  39 ] = tNeighbor->get_child(  37 );
                       tChildren[  42 ] = tNeighbor->get_child(  40 );
                       tChildren[  43 ] = tNeighbor->get_child(  41 );
                       tChildren[  44 ] = tNeighbor->get_child(  42 );
                       tChildren[  47 ] = tNeighbor->get_child(  45 );
                       tChildren[  48 ] = tNeighbor->get_child(  46 );
                       tChildren[  49 ] = tNeighbor->get_child(  47 );
                       tChildren[  52 ] = tNeighbor->get_child(  50 );
                       tChildren[  53 ] = tNeighbor->get_child(  51 );
                       tChildren[  54 ] = tNeighbor->get_child(  52 );
                       tChildren[  57 ] = tNeighbor->get_child(  55 );
                       tChildren[  58 ] = tNeighbor->get_child(  56 );
                       tChildren[  59 ] = tNeighbor->get_child(  57 );
                       tChildren[  62 ] = tNeighbor->get_child(  60 );
                       tChildren[  63 ] = tNeighbor->get_child(  61 );
                       tChildren[  64 ] = tNeighbor->get_child(  62 );
                       tChildren[  67 ] = tNeighbor->get_child(  65 );
                       tChildren[  68 ] = tNeighbor->get_child(  66 );
                       tChildren[  69 ] = tNeighbor->get_child(  67 );
                       tChildren[  72 ] = tNeighbor->get_child(  70 );
                       tChildren[  73 ] = tNeighbor->get_child(  71 );
                       tChildren[  74 ] = tNeighbor->get_child(  72 );
                       tChildren[  77 ] = tNeighbor->get_child(  75 );
                       tChildren[  78 ] = tNeighbor->get_child(  76 );
                       tChildren[  79 ] = tNeighbor->get_child(  77 );
                       tChildren[  82 ] = tNeighbor->get_child(  80 );
                       tChildren[  83 ] = tNeighbor->get_child(  81 );
                       tChildren[  84 ] = tNeighbor->get_child(  82 );
                       tChildren[  87 ] = tNeighbor->get_child(  85 );
                       tChildren[  88 ] = tNeighbor->get_child(  86 );
                       tChildren[  89 ] = tNeighbor->get_child(  87 );
                       tChildren[  92 ] = tNeighbor->get_child(  90 );
                       tChildren[  93 ] = tNeighbor->get_child(  91 );
                       tChildren[  94 ] = tNeighbor->get_child(  92 );
                       tChildren[  97 ] = tNeighbor->get_child(  95 );
                       tChildren[  98 ] = tNeighbor->get_child(  96 );
                       tChildren[  99 ] = tNeighbor->get_child(  97 );
                       tChildren[ 102 ] = tNeighbor->get_child( 100 );
                       tChildren[ 103 ] = tNeighbor->get_child( 101 );
                       tChildren[ 104 ] = tNeighbor->get_child( 102 );
                       tChildren[ 107 ] = tNeighbor->get_child( 105 );
                       tChildren[ 108 ] = tNeighbor->get_child( 106 );
                       tChildren[ 109 ] = tNeighbor->get_child( 107 );
                       tChildren[ 112 ] = tNeighbor->get_child( 110 );
                       tChildren[ 113 ] = tNeighbor->get_child( 111 );
                       tChildren[ 114 ] = tNeighbor->get_child( 112 );
                       tChildren[ 117 ] = tNeighbor->get_child( 115 );
                       tChildren[ 118 ] = tNeighbor->get_child( 116 );
                       tChildren[ 119 ] = tNeighbor->get_child( 117 );
                       tChildren[ 122 ] = tNeighbor->get_child( 120 );
                       tChildren[ 123 ] = tNeighbor->get_child( 121 );
                       tChildren[ 124 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 2 exists
                if( tNeighbors[ 2 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 2 ];

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
                       tChildren[  35 ] = tNeighbor->get_child(  25 );
                       tChildren[  36 ] = tNeighbor->get_child(  26 );
                       tChildren[  37 ] = tNeighbor->get_child(  27 );
                       tChildren[  38 ] = tNeighbor->get_child(  28 );
                       tChildren[  39 ] = tNeighbor->get_child(  29 );
                       tChildren[  40 ] = tNeighbor->get_child(  30 );
                       tChildren[  41 ] = tNeighbor->get_child(  31 );
                       tChildren[  42 ] = tNeighbor->get_child(  32 );
                       tChildren[  43 ] = tNeighbor->get_child(  33 );
                       tChildren[  44 ] = tNeighbor->get_child(  34 );
                       tChildren[  45 ] = tNeighbor->get_child(  35 );
                       tChildren[  46 ] = tNeighbor->get_child(  36 );
                       tChildren[  47 ] = tNeighbor->get_child(  37 );
                       tChildren[  48 ] = tNeighbor->get_child(  38 );
                       tChildren[  49 ] = tNeighbor->get_child(  39 );
                       tChildren[  60 ] = tNeighbor->get_child(  50 );
                       tChildren[  61 ] = tNeighbor->get_child(  51 );
                       tChildren[  62 ] = tNeighbor->get_child(  52 );
                       tChildren[  63 ] = tNeighbor->get_child(  53 );
                       tChildren[  64 ] = tNeighbor->get_child(  54 );
                       tChildren[  65 ] = tNeighbor->get_child(  55 );
                       tChildren[  66 ] = tNeighbor->get_child(  56 );
                       tChildren[  67 ] = tNeighbor->get_child(  57 );
                       tChildren[  68 ] = tNeighbor->get_child(  58 );
                       tChildren[  69 ] = tNeighbor->get_child(  59 );
                       tChildren[  70 ] = tNeighbor->get_child(  60 );
                       tChildren[  71 ] = tNeighbor->get_child(  61 );
                       tChildren[  72 ] = tNeighbor->get_child(  62 );
                       tChildren[  73 ] = tNeighbor->get_child(  63 );
                       tChildren[  74 ] = tNeighbor->get_child(  64 );
                       tChildren[  85 ] = tNeighbor->get_child(  75 );
                       tChildren[  86 ] = tNeighbor->get_child(  76 );
                       tChildren[  87 ] = tNeighbor->get_child(  77 );
                       tChildren[  88 ] = tNeighbor->get_child(  78 );
                       tChildren[  89 ] = tNeighbor->get_child(  79 );
                       tChildren[  90 ] = tNeighbor->get_child(  80 );
                       tChildren[  91 ] = tNeighbor->get_child(  81 );
                       tChildren[  92 ] = tNeighbor->get_child(  82 );
                       tChildren[  93 ] = tNeighbor->get_child(  83 );
                       tChildren[  94 ] = tNeighbor->get_child(  84 );
                       tChildren[  95 ] = tNeighbor->get_child(  85 );
                       tChildren[  96 ] = tNeighbor->get_child(  86 );
                       tChildren[  97 ] = tNeighbor->get_child(  87 );
                       tChildren[  98 ] = tNeighbor->get_child(  88 );
                       tChildren[  99 ] = tNeighbor->get_child(  89 );
                       tChildren[ 110 ] = tNeighbor->get_child( 100 );
                       tChildren[ 111 ] = tNeighbor->get_child( 101 );
                       tChildren[ 112 ] = tNeighbor->get_child( 102 );
                       tChildren[ 113 ] = tNeighbor->get_child( 103 );
                       tChildren[ 114 ] = tNeighbor->get_child( 104 );
                       tChildren[ 115 ] = tNeighbor->get_child( 105 );
                       tChildren[ 116 ] = tNeighbor->get_child( 106 );
                       tChildren[ 117 ] = tNeighbor->get_child( 107 );
                       tChildren[ 118 ] = tNeighbor->get_child( 108 );
                       tChildren[ 119 ] = tNeighbor->get_child( 109 );
                       tChildren[ 120 ] = tNeighbor->get_child( 110 );
                       tChildren[ 121 ] = tNeighbor->get_child( 111 );
                       tChildren[ 122 ] = tNeighbor->get_child( 112 );
                       tChildren[ 123 ] = tNeighbor->get_child( 113 );
                       tChildren[ 124 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 3 exists
                if( tNeighbors[ 3 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 3 ];

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
                       tChildren[  25 ] = tNeighbor->get_child(  27 );
                       tChildren[  26 ] = tNeighbor->get_child(  28 );
                       tChildren[  27 ] = tNeighbor->get_child(  29 );
                       tChildren[  30 ] = tNeighbor->get_child(  32 );
                       tChildren[  31 ] = tNeighbor->get_child(  33 );
                       tChildren[  32 ] = tNeighbor->get_child(  34 );
                       tChildren[  35 ] = tNeighbor->get_child(  37 );
                       tChildren[  36 ] = tNeighbor->get_child(  38 );
                       tChildren[  37 ] = tNeighbor->get_child(  39 );
                       tChildren[  40 ] = tNeighbor->get_child(  42 );
                       tChildren[  41 ] = tNeighbor->get_child(  43 );
                       tChildren[  42 ] = tNeighbor->get_child(  44 );
                       tChildren[  45 ] = tNeighbor->get_child(  47 );
                       tChildren[  46 ] = tNeighbor->get_child(  48 );
                       tChildren[  47 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  52 );
                       tChildren[  51 ] = tNeighbor->get_child(  53 );
                       tChildren[  52 ] = tNeighbor->get_child(  54 );
                       tChildren[  55 ] = tNeighbor->get_child(  57 );
                       tChildren[  56 ] = tNeighbor->get_child(  58 );
                       tChildren[  57 ] = tNeighbor->get_child(  59 );
                       tChildren[  60 ] = tNeighbor->get_child(  62 );
                       tChildren[  61 ] = tNeighbor->get_child(  63 );
                       tChildren[  62 ] = tNeighbor->get_child(  64 );
                       tChildren[  65 ] = tNeighbor->get_child(  67 );
                       tChildren[  66 ] = tNeighbor->get_child(  68 );
                       tChildren[  67 ] = tNeighbor->get_child(  69 );
                       tChildren[  70 ] = tNeighbor->get_child(  72 );
                       tChildren[  71 ] = tNeighbor->get_child(  73 );
                       tChildren[  72 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  77 );
                       tChildren[  76 ] = tNeighbor->get_child(  78 );
                       tChildren[  77 ] = tNeighbor->get_child(  79 );
                       tChildren[  80 ] = tNeighbor->get_child(  82 );
                       tChildren[  81 ] = tNeighbor->get_child(  83 );
                       tChildren[  82 ] = tNeighbor->get_child(  84 );
                       tChildren[  85 ] = tNeighbor->get_child(  87 );
                       tChildren[  86 ] = tNeighbor->get_child(  88 );
                       tChildren[  87 ] = tNeighbor->get_child(  89 );
                       tChildren[  90 ] = tNeighbor->get_child(  92 );
                       tChildren[  91 ] = tNeighbor->get_child(  93 );
                       tChildren[  92 ] = tNeighbor->get_child(  94 );
                       tChildren[  95 ] = tNeighbor->get_child(  97 );
                       tChildren[  96 ] = tNeighbor->get_child(  98 );
                       tChildren[  97 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 102 );
                       tChildren[ 101 ] = tNeighbor->get_child( 103 );
                       tChildren[ 102 ] = tNeighbor->get_child( 104 );
                       tChildren[ 105 ] = tNeighbor->get_child( 107 );
                       tChildren[ 106 ] = tNeighbor->get_child( 108 );
                       tChildren[ 107 ] = tNeighbor->get_child( 109 );
                       tChildren[ 110 ] = tNeighbor->get_child( 112 );
                       tChildren[ 111 ] = tNeighbor->get_child( 113 );
                       tChildren[ 112 ] = tNeighbor->get_child( 114 );
                       tChildren[ 115 ] = tNeighbor->get_child( 117 );
                       tChildren[ 116 ] = tNeighbor->get_child( 118 );
                       tChildren[ 117 ] = tNeighbor->get_child( 119 );
                       tChildren[ 120 ] = tNeighbor->get_child( 122 );
                       tChildren[ 121 ] = tNeighbor->get_child( 123 );
                       tChildren[ 122 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 4 exists
                if( tNeighbors[ 4 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 4 ];

                    // test if neighbor 4 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  50 );
                       tChildren[   1 ] = tNeighbor->get_child(  51 );
                       tChildren[   2 ] = tNeighbor->get_child(  52 );
                       tChildren[   3 ] = tNeighbor->get_child(  53 );
                       tChildren[   4 ] = tNeighbor->get_child(  54 );
                       tChildren[   5 ] = tNeighbor->get_child(  55 );
                       tChildren[   6 ] = tNeighbor->get_child(  56 );
                       tChildren[   7 ] = tNeighbor->get_child(  57 );
                       tChildren[   8 ] = tNeighbor->get_child(  58 );
                       tChildren[   9 ] = tNeighbor->get_child(  59 );
                       tChildren[  10 ] = tNeighbor->get_child(  60 );
                       tChildren[  11 ] = tNeighbor->get_child(  61 );
                       tChildren[  12 ] = tNeighbor->get_child(  62 );
                       tChildren[  13 ] = tNeighbor->get_child(  63 );
                       tChildren[  14 ] = tNeighbor->get_child(  64 );
                       tChildren[  15 ] = tNeighbor->get_child(  65 );
                       tChildren[  16 ] = tNeighbor->get_child(  66 );
                       tChildren[  17 ] = tNeighbor->get_child(  67 );
                       tChildren[  18 ] = tNeighbor->get_child(  68 );
                       tChildren[  19 ] = tNeighbor->get_child(  69 );
                       tChildren[  20 ] = tNeighbor->get_child(  70 );
                       tChildren[  21 ] = tNeighbor->get_child(  71 );
                       tChildren[  22 ] = tNeighbor->get_child(  72 );
                       tChildren[  23 ] = tNeighbor->get_child(  73 );
                       tChildren[  24 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  75 );
                       tChildren[  26 ] = tNeighbor->get_child(  76 );
                       tChildren[  27 ] = tNeighbor->get_child(  77 );
                       tChildren[  28 ] = tNeighbor->get_child(  78 );
                       tChildren[  29 ] = tNeighbor->get_child(  79 );
                       tChildren[  30 ] = tNeighbor->get_child(  80 );
                       tChildren[  31 ] = tNeighbor->get_child(  81 );
                       tChildren[  32 ] = tNeighbor->get_child(  82 );
                       tChildren[  33 ] = tNeighbor->get_child(  83 );
                       tChildren[  34 ] = tNeighbor->get_child(  84 );
                       tChildren[  35 ] = tNeighbor->get_child(  85 );
                       tChildren[  36 ] = tNeighbor->get_child(  86 );
                       tChildren[  37 ] = tNeighbor->get_child(  87 );
                       tChildren[  38 ] = tNeighbor->get_child(  88 );
                       tChildren[  39 ] = tNeighbor->get_child(  89 );
                       tChildren[  40 ] = tNeighbor->get_child(  90 );
                       tChildren[  41 ] = tNeighbor->get_child(  91 );
                       tChildren[  42 ] = tNeighbor->get_child(  92 );
                       tChildren[  43 ] = tNeighbor->get_child(  93 );
                       tChildren[  44 ] = tNeighbor->get_child(  94 );
                       tChildren[  45 ] = tNeighbor->get_child(  95 );
                       tChildren[  46 ] = tNeighbor->get_child(  96 );
                       tChildren[  47 ] = tNeighbor->get_child(  97 );
                       tChildren[  48 ] = tNeighbor->get_child(  98 );
                       tChildren[  49 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 100 );
                       tChildren[  51 ] = tNeighbor->get_child( 101 );
                       tChildren[  52 ] = tNeighbor->get_child( 102 );
                       tChildren[  53 ] = tNeighbor->get_child( 103 );
                       tChildren[  54 ] = tNeighbor->get_child( 104 );
                       tChildren[  55 ] = tNeighbor->get_child( 105 );
                       tChildren[  56 ] = tNeighbor->get_child( 106 );
                       tChildren[  57 ] = tNeighbor->get_child( 107 );
                       tChildren[  58 ] = tNeighbor->get_child( 108 );
                       tChildren[  59 ] = tNeighbor->get_child( 109 );
                       tChildren[  60 ] = tNeighbor->get_child( 110 );
                       tChildren[  61 ] = tNeighbor->get_child( 111 );
                       tChildren[  62 ] = tNeighbor->get_child( 112 );
                       tChildren[  63 ] = tNeighbor->get_child( 113 );
                       tChildren[  64 ] = tNeighbor->get_child( 114 );
                       tChildren[  65 ] = tNeighbor->get_child( 115 );
                       tChildren[  66 ] = tNeighbor->get_child( 116 );
                       tChildren[  67 ] = tNeighbor->get_child( 117 );
                       tChildren[  68 ] = tNeighbor->get_child( 118 );
                       tChildren[  69 ] = tNeighbor->get_child( 119 );
                       tChildren[  70 ] = tNeighbor->get_child( 120 );
                       tChildren[  71 ] = tNeighbor->get_child( 121 );
                       tChildren[  72 ] = tNeighbor->get_child( 122 );
                       tChildren[  73 ] = tNeighbor->get_child( 123 );
                       tChildren[  74 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 5 exists
                if( tNeighbors[ 5 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 5 ];

                    // test if neighbor 5 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(   0 );
                       tChildren[  51 ] = tNeighbor->get_child(   1 );
                       tChildren[  52 ] = tNeighbor->get_child(   2 );
                       tChildren[  53 ] = tNeighbor->get_child(   3 );
                       tChildren[  54 ] = tNeighbor->get_child(   4 );
                       tChildren[  55 ] = tNeighbor->get_child(   5 );
                       tChildren[  56 ] = tNeighbor->get_child(   6 );
                       tChildren[  57 ] = tNeighbor->get_child(   7 );
                       tChildren[  58 ] = tNeighbor->get_child(   8 );
                       tChildren[  59 ] = tNeighbor->get_child(   9 );
                       tChildren[  60 ] = tNeighbor->get_child(  10 );
                       tChildren[  61 ] = tNeighbor->get_child(  11 );
                       tChildren[  62 ] = tNeighbor->get_child(  12 );
                       tChildren[  63 ] = tNeighbor->get_child(  13 );
                       tChildren[  64 ] = tNeighbor->get_child(  14 );
                       tChildren[  65 ] = tNeighbor->get_child(  15 );
                       tChildren[  66 ] = tNeighbor->get_child(  16 );
                       tChildren[  67 ] = tNeighbor->get_child(  17 );
                       tChildren[  68 ] = tNeighbor->get_child(  18 );
                       tChildren[  69 ] = tNeighbor->get_child(  19 );
                       tChildren[  70 ] = tNeighbor->get_child(  20 );
                       tChildren[  71 ] = tNeighbor->get_child(  21 );
                       tChildren[  72 ] = tNeighbor->get_child(  22 );
                       tChildren[  73 ] = tNeighbor->get_child(  23 );
                       tChildren[  74 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  25 );
                       tChildren[  76 ] = tNeighbor->get_child(  26 );
                       tChildren[  77 ] = tNeighbor->get_child(  27 );
                       tChildren[  78 ] = tNeighbor->get_child(  28 );
                       tChildren[  79 ] = tNeighbor->get_child(  29 );
                       tChildren[  80 ] = tNeighbor->get_child(  30 );
                       tChildren[  81 ] = tNeighbor->get_child(  31 );
                       tChildren[  82 ] = tNeighbor->get_child(  32 );
                       tChildren[  83 ] = tNeighbor->get_child(  33 );
                       tChildren[  84 ] = tNeighbor->get_child(  34 );
                       tChildren[  85 ] = tNeighbor->get_child(  35 );
                       tChildren[  86 ] = tNeighbor->get_child(  36 );
                       tChildren[  87 ] = tNeighbor->get_child(  37 );
                       tChildren[  88 ] = tNeighbor->get_child(  38 );
                       tChildren[  89 ] = tNeighbor->get_child(  39 );
                       tChildren[  90 ] = tNeighbor->get_child(  40 );
                       tChildren[  91 ] = tNeighbor->get_child(  41 );
                       tChildren[  92 ] = tNeighbor->get_child(  42 );
                       tChildren[  93 ] = tNeighbor->get_child(  43 );
                       tChildren[  94 ] = tNeighbor->get_child(  44 );
                       tChildren[  95 ] = tNeighbor->get_child(  45 );
                       tChildren[  96 ] = tNeighbor->get_child(  46 );
                       tChildren[  97 ] = tNeighbor->get_child(  47 );
                       tChildren[  98 ] = tNeighbor->get_child(  48 );
                       tChildren[  99 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  50 );
                       tChildren[ 101 ] = tNeighbor->get_child(  51 );
                       tChildren[ 102 ] = tNeighbor->get_child(  52 );
                       tChildren[ 103 ] = tNeighbor->get_child(  53 );
                       tChildren[ 104 ] = tNeighbor->get_child(  54 );
                       tChildren[ 105 ] = tNeighbor->get_child(  55 );
                       tChildren[ 106 ] = tNeighbor->get_child(  56 );
                       tChildren[ 107 ] = tNeighbor->get_child(  57 );
                       tChildren[ 108 ] = tNeighbor->get_child(  58 );
                       tChildren[ 109 ] = tNeighbor->get_child(  59 );
                       tChildren[ 110 ] = tNeighbor->get_child(  60 );
                       tChildren[ 111 ] = tNeighbor->get_child(  61 );
                       tChildren[ 112 ] = tNeighbor->get_child(  62 );
                       tChildren[ 113 ] = tNeighbor->get_child(  63 );
                       tChildren[ 114 ] = tNeighbor->get_child(  64 );
                       tChildren[ 115 ] = tNeighbor->get_child(  65 );
                       tChildren[ 116 ] = tNeighbor->get_child(  66 );
                       tChildren[ 117 ] = tNeighbor->get_child(  67 );
                       tChildren[ 118 ] = tNeighbor->get_child(  68 );
                       tChildren[ 119 ] = tNeighbor->get_child(  69 );
                       tChildren[ 120 ] = tNeighbor->get_child(  70 );
                       tChildren[ 121 ] = tNeighbor->get_child(  71 );
                       tChildren[ 122 ] = tNeighbor->get_child(  72 );
                       tChildren[ 123 ] = tNeighbor->get_child(  73 );
                       tChildren[ 124 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 6 exists
                if( tNeighbors[ 6 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 6 ];

                    // test if neighbor 6 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  60 );
                       tChildren[   1 ] = tNeighbor->get_child(  61 );
                       tChildren[   2 ] = tNeighbor->get_child(  62 );
                       tChildren[   3 ] = tNeighbor->get_child(  63 );
                       tChildren[   4 ] = tNeighbor->get_child(  64 );
                       tChildren[   5 ] = tNeighbor->get_child(  65 );
                       tChildren[   6 ] = tNeighbor->get_child(  66 );
                       tChildren[   7 ] = tNeighbor->get_child(  67 );
                       tChildren[   8 ] = tNeighbor->get_child(  68 );
                       tChildren[   9 ] = tNeighbor->get_child(  69 );
                       tChildren[  10 ] = tNeighbor->get_child(  70 );
                       tChildren[  11 ] = tNeighbor->get_child(  71 );
                       tChildren[  12 ] = tNeighbor->get_child(  72 );
                       tChildren[  13 ] = tNeighbor->get_child(  73 );
                       tChildren[  14 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  85 );
                       tChildren[  26 ] = tNeighbor->get_child(  86 );
                       tChildren[  27 ] = tNeighbor->get_child(  87 );
                       tChildren[  28 ] = tNeighbor->get_child(  88 );
                       tChildren[  29 ] = tNeighbor->get_child(  89 );
                       tChildren[  30 ] = tNeighbor->get_child(  90 );
                       tChildren[  31 ] = tNeighbor->get_child(  91 );
                       tChildren[  32 ] = tNeighbor->get_child(  92 );
                       tChildren[  33 ] = tNeighbor->get_child(  93 );
                       tChildren[  34 ] = tNeighbor->get_child(  94 );
                       tChildren[  35 ] = tNeighbor->get_child(  95 );
                       tChildren[  36 ] = tNeighbor->get_child(  96 );
                       tChildren[  37 ] = tNeighbor->get_child(  97 );
                       tChildren[  38 ] = tNeighbor->get_child(  98 );
                       tChildren[  39 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 110 );
                       tChildren[  51 ] = tNeighbor->get_child( 111 );
                       tChildren[  52 ] = tNeighbor->get_child( 112 );
                       tChildren[  53 ] = tNeighbor->get_child( 113 );
                       tChildren[  54 ] = tNeighbor->get_child( 114 );
                       tChildren[  55 ] = tNeighbor->get_child( 115 );
                       tChildren[  56 ] = tNeighbor->get_child( 116 );
                       tChildren[  57 ] = tNeighbor->get_child( 117 );
                       tChildren[  58 ] = tNeighbor->get_child( 118 );
                       tChildren[  59 ] = tNeighbor->get_child( 119 );
                       tChildren[  60 ] = tNeighbor->get_child( 120 );
                       tChildren[  61 ] = tNeighbor->get_child( 121 );
                       tChildren[  62 ] = tNeighbor->get_child( 122 );
                       tChildren[  63 ] = tNeighbor->get_child( 123 );
                       tChildren[  64 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 7 exists
                if( tNeighbors[ 7 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 7 ];

                    // test if neighbor 7 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  50 );
                       tChildren[   3 ] = tNeighbor->get_child(  51 );
                       tChildren[   4 ] = tNeighbor->get_child(  52 );
                       tChildren[   7 ] = tNeighbor->get_child(  55 );
                       tChildren[   8 ] = tNeighbor->get_child(  56 );
                       tChildren[   9 ] = tNeighbor->get_child(  57 );
                       tChildren[  12 ] = tNeighbor->get_child(  60 );
                       tChildren[  13 ] = tNeighbor->get_child(  61 );
                       tChildren[  14 ] = tNeighbor->get_child(  62 );
                       tChildren[  17 ] = tNeighbor->get_child(  65 );
                       tChildren[  18 ] = tNeighbor->get_child(  66 );
                       tChildren[  19 ] = tNeighbor->get_child(  67 );
                       tChildren[  22 ] = tNeighbor->get_child(  70 );
                       tChildren[  23 ] = tNeighbor->get_child(  71 );
                       tChildren[  24 ] = tNeighbor->get_child(  72 );
                       tChildren[  27 ] = tNeighbor->get_child(  75 );
                       tChildren[  28 ] = tNeighbor->get_child(  76 );
                       tChildren[  29 ] = tNeighbor->get_child(  77 );
                       tChildren[  32 ] = tNeighbor->get_child(  80 );
                       tChildren[  33 ] = tNeighbor->get_child(  81 );
                       tChildren[  34 ] = tNeighbor->get_child(  82 );
                       tChildren[  37 ] = tNeighbor->get_child(  85 );
                       tChildren[  38 ] = tNeighbor->get_child(  86 );
                       tChildren[  39 ] = tNeighbor->get_child(  87 );
                       tChildren[  42 ] = tNeighbor->get_child(  90 );
                       tChildren[  43 ] = tNeighbor->get_child(  91 );
                       tChildren[  44 ] = tNeighbor->get_child(  92 );
                       tChildren[  47 ] = tNeighbor->get_child(  95 );
                       tChildren[  48 ] = tNeighbor->get_child(  96 );
                       tChildren[  49 ] = tNeighbor->get_child(  97 );
                       tChildren[  52 ] = tNeighbor->get_child( 100 );
                       tChildren[  53 ] = tNeighbor->get_child( 101 );
                       tChildren[  54 ] = tNeighbor->get_child( 102 );
                       tChildren[  57 ] = tNeighbor->get_child( 105 );
                       tChildren[  58 ] = tNeighbor->get_child( 106 );
                       tChildren[  59 ] = tNeighbor->get_child( 107 );
                       tChildren[  62 ] = tNeighbor->get_child( 110 );
                       tChildren[  63 ] = tNeighbor->get_child( 111 );
                       tChildren[  64 ] = tNeighbor->get_child( 112 );
                       tChildren[  67 ] = tNeighbor->get_child( 115 );
                       tChildren[  68 ] = tNeighbor->get_child( 116 );
                       tChildren[  69 ] = tNeighbor->get_child( 117 );
                       tChildren[  72 ] = tNeighbor->get_child( 120 );
                       tChildren[  73 ] = tNeighbor->get_child( 121 );
                       tChildren[  74 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 8 exists
                if( tNeighbors[ 8 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 8 ];

                    // test if neighbor 8 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(  50 );
                       tChildren[  11 ] = tNeighbor->get_child(  51 );
                       tChildren[  12 ] = tNeighbor->get_child(  52 );
                       tChildren[  13 ] = tNeighbor->get_child(  53 );
                       tChildren[  14 ] = tNeighbor->get_child(  54 );
                       tChildren[  15 ] = tNeighbor->get_child(  55 );
                       tChildren[  16 ] = tNeighbor->get_child(  56 );
                       tChildren[  17 ] = tNeighbor->get_child(  57 );
                       tChildren[  18 ] = tNeighbor->get_child(  58 );
                       tChildren[  19 ] = tNeighbor->get_child(  59 );
                       tChildren[  20 ] = tNeighbor->get_child(  60 );
                       tChildren[  21 ] = tNeighbor->get_child(  61 );
                       tChildren[  22 ] = tNeighbor->get_child(  62 );
                       tChildren[  23 ] = tNeighbor->get_child(  63 );
                       tChildren[  24 ] = tNeighbor->get_child(  64 );
                       tChildren[  35 ] = tNeighbor->get_child(  75 );
                       tChildren[  36 ] = tNeighbor->get_child(  76 );
                       tChildren[  37 ] = tNeighbor->get_child(  77 );
                       tChildren[  38 ] = tNeighbor->get_child(  78 );
                       tChildren[  39 ] = tNeighbor->get_child(  79 );
                       tChildren[  40 ] = tNeighbor->get_child(  80 );
                       tChildren[  41 ] = tNeighbor->get_child(  81 );
                       tChildren[  42 ] = tNeighbor->get_child(  82 );
                       tChildren[  43 ] = tNeighbor->get_child(  83 );
                       tChildren[  44 ] = tNeighbor->get_child(  84 );
                       tChildren[  45 ] = tNeighbor->get_child(  85 );
                       tChildren[  46 ] = tNeighbor->get_child(  86 );
                       tChildren[  47 ] = tNeighbor->get_child(  87 );
                       tChildren[  48 ] = tNeighbor->get_child(  88 );
                       tChildren[  49 ] = tNeighbor->get_child(  89 );
                       tChildren[  60 ] = tNeighbor->get_child( 100 );
                       tChildren[  61 ] = tNeighbor->get_child( 101 );
                       tChildren[  62 ] = tNeighbor->get_child( 102 );
                       tChildren[  63 ] = tNeighbor->get_child( 103 );
                       tChildren[  64 ] = tNeighbor->get_child( 104 );
                       tChildren[  65 ] = tNeighbor->get_child( 105 );
                       tChildren[  66 ] = tNeighbor->get_child( 106 );
                       tChildren[  67 ] = tNeighbor->get_child( 107 );
                       tChildren[  68 ] = tNeighbor->get_child( 108 );
                       tChildren[  69 ] = tNeighbor->get_child( 109 );
                       tChildren[  70 ] = tNeighbor->get_child( 110 );
                       tChildren[  71 ] = tNeighbor->get_child( 111 );
                       tChildren[  72 ] = tNeighbor->get_child( 112 );
                       tChildren[  73 ] = tNeighbor->get_child( 113 );
                       tChildren[  74 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 9 exists
                if( tNeighbors[ 9 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 9 ];

                    // test if neighbor 9 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  52 );
                       tChildren[   1 ] = tNeighbor->get_child(  53 );
                       tChildren[   2 ] = tNeighbor->get_child(  54 );
                       tChildren[   5 ] = tNeighbor->get_child(  57 );
                       tChildren[   6 ] = tNeighbor->get_child(  58 );
                       tChildren[   7 ] = tNeighbor->get_child(  59 );
                       tChildren[  10 ] = tNeighbor->get_child(  62 );
                       tChildren[  11 ] = tNeighbor->get_child(  63 );
                       tChildren[  12 ] = tNeighbor->get_child(  64 );
                       tChildren[  15 ] = tNeighbor->get_child(  67 );
                       tChildren[  16 ] = tNeighbor->get_child(  68 );
                       tChildren[  17 ] = tNeighbor->get_child(  69 );
                       tChildren[  20 ] = tNeighbor->get_child(  72 );
                       tChildren[  21 ] = tNeighbor->get_child(  73 );
                       tChildren[  22 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  77 );
                       tChildren[  26 ] = tNeighbor->get_child(  78 );
                       tChildren[  27 ] = tNeighbor->get_child(  79 );
                       tChildren[  30 ] = tNeighbor->get_child(  82 );
                       tChildren[  31 ] = tNeighbor->get_child(  83 );
                       tChildren[  32 ] = tNeighbor->get_child(  84 );
                       tChildren[  35 ] = tNeighbor->get_child(  87 );
                       tChildren[  36 ] = tNeighbor->get_child(  88 );
                       tChildren[  37 ] = tNeighbor->get_child(  89 );
                       tChildren[  40 ] = tNeighbor->get_child(  92 );
                       tChildren[  41 ] = tNeighbor->get_child(  93 );
                       tChildren[  42 ] = tNeighbor->get_child(  94 );
                       tChildren[  45 ] = tNeighbor->get_child(  97 );
                       tChildren[  46 ] = tNeighbor->get_child(  98 );
                       tChildren[  47 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 102 );
                       tChildren[  51 ] = tNeighbor->get_child( 103 );
                       tChildren[  52 ] = tNeighbor->get_child( 104 );
                       tChildren[  55 ] = tNeighbor->get_child( 107 );
                       tChildren[  56 ] = tNeighbor->get_child( 108 );
                       tChildren[  57 ] = tNeighbor->get_child( 109 );
                       tChildren[  60 ] = tNeighbor->get_child( 112 );
                       tChildren[  61 ] = tNeighbor->get_child( 113 );
                       tChildren[  62 ] = tNeighbor->get_child( 114 );
                       tChildren[  65 ] = tNeighbor->get_child( 117 );
                       tChildren[  66 ] = tNeighbor->get_child( 118 );
                       tChildren[  67 ] = tNeighbor->get_child( 119 );
                       tChildren[  70 ] = tNeighbor->get_child( 122 );
                       tChildren[  71 ] = tNeighbor->get_child( 123 );
                       tChildren[  72 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 10 exists
                if( tNeighbors[ 10 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 10 ];

                    // test if neighbor 10 has children and copy them if so
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
                       tChildren[  25 ] = tNeighbor->get_child(  37 );
                       tChildren[  26 ] = tNeighbor->get_child(  38 );
                       tChildren[  27 ] = tNeighbor->get_child(  39 );
                       tChildren[  30 ] = tNeighbor->get_child(  42 );
                       tChildren[  31 ] = tNeighbor->get_child(  43 );
                       tChildren[  32 ] = tNeighbor->get_child(  44 );
                       tChildren[  35 ] = tNeighbor->get_child(  47 );
                       tChildren[  36 ] = tNeighbor->get_child(  48 );
                       tChildren[  37 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  62 );
                       tChildren[  51 ] = tNeighbor->get_child(  63 );
                       tChildren[  52 ] = tNeighbor->get_child(  64 );
                       tChildren[  55 ] = tNeighbor->get_child(  67 );
                       tChildren[  56 ] = tNeighbor->get_child(  68 );
                       tChildren[  57 ] = tNeighbor->get_child(  69 );
                       tChildren[  60 ] = tNeighbor->get_child(  72 );
                       tChildren[  61 ] = tNeighbor->get_child(  73 );
                       tChildren[  62 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  87 );
                       tChildren[  76 ] = tNeighbor->get_child(  88 );
                       tChildren[  77 ] = tNeighbor->get_child(  89 );
                       tChildren[  80 ] = tNeighbor->get_child(  92 );
                       tChildren[  81 ] = tNeighbor->get_child(  93 );
                       tChildren[  82 ] = tNeighbor->get_child(  94 );
                       tChildren[  85 ] = tNeighbor->get_child(  97 );
                       tChildren[  86 ] = tNeighbor->get_child(  98 );
                       tChildren[  87 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 112 );
                       tChildren[ 101 ] = tNeighbor->get_child( 113 );
                       tChildren[ 102 ] = tNeighbor->get_child( 114 );
                       tChildren[ 105 ] = tNeighbor->get_child( 117 );
                       tChildren[ 106 ] = tNeighbor->get_child( 118 );
                       tChildren[ 107 ] = tNeighbor->get_child( 119 );
                       tChildren[ 110 ] = tNeighbor->get_child( 122 );
                       tChildren[ 111 ] = tNeighbor->get_child( 123 );
                       tChildren[ 112 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 11 exists
                if( tNeighbors[ 11 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 11 ];

                    // test if neighbor 11 has children and copy them if so
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
                       tChildren[  27 ] = tNeighbor->get_child(  35 );
                       tChildren[  28 ] = tNeighbor->get_child(  36 );
                       tChildren[  29 ] = tNeighbor->get_child(  37 );
                       tChildren[  32 ] = tNeighbor->get_child(  40 );
                       tChildren[  33 ] = tNeighbor->get_child(  41 );
                       tChildren[  34 ] = tNeighbor->get_child(  42 );
                       tChildren[  37 ] = tNeighbor->get_child(  45 );
                       tChildren[  38 ] = tNeighbor->get_child(  46 );
                       tChildren[  39 ] = tNeighbor->get_child(  47 );
                       tChildren[  52 ] = tNeighbor->get_child(  60 );
                       tChildren[  53 ] = tNeighbor->get_child(  61 );
                       tChildren[  54 ] = tNeighbor->get_child(  62 );
                       tChildren[  57 ] = tNeighbor->get_child(  65 );
                       tChildren[  58 ] = tNeighbor->get_child(  66 );
                       tChildren[  59 ] = tNeighbor->get_child(  67 );
                       tChildren[  62 ] = tNeighbor->get_child(  70 );
                       tChildren[  63 ] = tNeighbor->get_child(  71 );
                       tChildren[  64 ] = tNeighbor->get_child(  72 );
                       tChildren[  77 ] = tNeighbor->get_child(  85 );
                       tChildren[  78 ] = tNeighbor->get_child(  86 );
                       tChildren[  79 ] = tNeighbor->get_child(  87 );
                       tChildren[  82 ] = tNeighbor->get_child(  90 );
                       tChildren[  83 ] = tNeighbor->get_child(  91 );
                       tChildren[  84 ] = tNeighbor->get_child(  92 );
                       tChildren[  87 ] = tNeighbor->get_child(  95 );
                       tChildren[  88 ] = tNeighbor->get_child(  96 );
                       tChildren[  89 ] = tNeighbor->get_child(  97 );
                       tChildren[ 102 ] = tNeighbor->get_child( 110 );
                       tChildren[ 103 ] = tNeighbor->get_child( 111 );
                       tChildren[ 104 ] = tNeighbor->get_child( 112 );
                       tChildren[ 107 ] = tNeighbor->get_child( 115 );
                       tChildren[ 108 ] = tNeighbor->get_child( 116 );
                       tChildren[ 109 ] = tNeighbor->get_child( 117 );
                       tChildren[ 112 ] = tNeighbor->get_child( 120 );
                       tChildren[ 113 ] = tNeighbor->get_child( 121 );
                       tChildren[ 114 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 12 exists
                if( tNeighbors[ 12 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 12 ];

                    // test if neighbor 12 has children and copy them if so
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
                       tChildren[  37 ] = tNeighbor->get_child(  25 );
                       tChildren[  38 ] = tNeighbor->get_child(  26 );
                       tChildren[  39 ] = tNeighbor->get_child(  27 );
                       tChildren[  42 ] = tNeighbor->get_child(  30 );
                       tChildren[  43 ] = tNeighbor->get_child(  31 );
                       tChildren[  44 ] = tNeighbor->get_child(  32 );
                       tChildren[  47 ] = tNeighbor->get_child(  35 );
                       tChildren[  48 ] = tNeighbor->get_child(  36 );
                       tChildren[  49 ] = tNeighbor->get_child(  37 );
                       tChildren[  62 ] = tNeighbor->get_child(  50 );
                       tChildren[  63 ] = tNeighbor->get_child(  51 );
                       tChildren[  64 ] = tNeighbor->get_child(  52 );
                       tChildren[  67 ] = tNeighbor->get_child(  55 );
                       tChildren[  68 ] = tNeighbor->get_child(  56 );
                       tChildren[  69 ] = tNeighbor->get_child(  57 );
                       tChildren[  72 ] = tNeighbor->get_child(  60 );
                       tChildren[  73 ] = tNeighbor->get_child(  61 );
                       tChildren[  74 ] = tNeighbor->get_child(  62 );
                       tChildren[  87 ] = tNeighbor->get_child(  75 );
                       tChildren[  88 ] = tNeighbor->get_child(  76 );
                       tChildren[  89 ] = tNeighbor->get_child(  77 );
                       tChildren[  92 ] = tNeighbor->get_child(  80 );
                       tChildren[  93 ] = tNeighbor->get_child(  81 );
                       tChildren[  94 ] = tNeighbor->get_child(  82 );
                       tChildren[  97 ] = tNeighbor->get_child(  85 );
                       tChildren[  98 ] = tNeighbor->get_child(  86 );
                       tChildren[  99 ] = tNeighbor->get_child(  87 );
                       tChildren[ 112 ] = tNeighbor->get_child( 100 );
                       tChildren[ 113 ] = tNeighbor->get_child( 101 );
                       tChildren[ 114 ] = tNeighbor->get_child( 102 );
                       tChildren[ 117 ] = tNeighbor->get_child( 105 );
                       tChildren[ 118 ] = tNeighbor->get_child( 106 );
                       tChildren[ 119 ] = tNeighbor->get_child( 107 );
                       tChildren[ 122 ] = tNeighbor->get_child( 110 );
                       tChildren[ 123 ] = tNeighbor->get_child( 111 );
                       tChildren[ 124 ] = tNeighbor->get_child( 112 );
                    }
                }

                // test if neighbor 13 exists
                if( tNeighbors[ 13 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 13 ];

                    // test if neighbor 13 has children and copy them if so
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
                       tChildren[  35 ] = tNeighbor->get_child(  27 );
                       tChildren[  36 ] = tNeighbor->get_child(  28 );
                       tChildren[  37 ] = tNeighbor->get_child(  29 );
                       tChildren[  40 ] = tNeighbor->get_child(  32 );
                       tChildren[  41 ] = tNeighbor->get_child(  33 );
                       tChildren[  42 ] = tNeighbor->get_child(  34 );
                       tChildren[  45 ] = tNeighbor->get_child(  37 );
                       tChildren[  46 ] = tNeighbor->get_child(  38 );
                       tChildren[  47 ] = tNeighbor->get_child(  39 );
                       tChildren[  60 ] = tNeighbor->get_child(  52 );
                       tChildren[  61 ] = tNeighbor->get_child(  53 );
                       tChildren[  62 ] = tNeighbor->get_child(  54 );
                       tChildren[  65 ] = tNeighbor->get_child(  57 );
                       tChildren[  66 ] = tNeighbor->get_child(  58 );
                       tChildren[  67 ] = tNeighbor->get_child(  59 );
                       tChildren[  70 ] = tNeighbor->get_child(  62 );
                       tChildren[  71 ] = tNeighbor->get_child(  63 );
                       tChildren[  72 ] = tNeighbor->get_child(  64 );
                       tChildren[  85 ] = tNeighbor->get_child(  77 );
                       tChildren[  86 ] = tNeighbor->get_child(  78 );
                       tChildren[  87 ] = tNeighbor->get_child(  79 );
                       tChildren[  90 ] = tNeighbor->get_child(  82 );
                       tChildren[  91 ] = tNeighbor->get_child(  83 );
                       tChildren[  92 ] = tNeighbor->get_child(  84 );
                       tChildren[  95 ] = tNeighbor->get_child(  87 );
                       tChildren[  96 ] = tNeighbor->get_child(  88 );
                       tChildren[  97 ] = tNeighbor->get_child(  89 );
                       tChildren[ 110 ] = tNeighbor->get_child( 102 );
                       tChildren[ 111 ] = tNeighbor->get_child( 103 );
                       tChildren[ 112 ] = tNeighbor->get_child( 104 );
                       tChildren[ 115 ] = tNeighbor->get_child( 107 );
                       tChildren[ 116 ] = tNeighbor->get_child( 108 );
                       tChildren[ 117 ] = tNeighbor->get_child( 109 );
                       tChildren[ 120 ] = tNeighbor->get_child( 112 );
                       tChildren[ 121 ] = tNeighbor->get_child( 113 );
                       tChildren[ 122 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 14 exists
                if( tNeighbors[ 14 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 14 ];

                    // test if neighbor 14 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  10 );
                       tChildren[  51 ] = tNeighbor->get_child(  11 );
                       tChildren[  52 ] = tNeighbor->get_child(  12 );
                       tChildren[  53 ] = tNeighbor->get_child(  13 );
                       tChildren[  54 ] = tNeighbor->get_child(  14 );
                       tChildren[  55 ] = tNeighbor->get_child(  15 );
                       tChildren[  56 ] = tNeighbor->get_child(  16 );
                       tChildren[  57 ] = tNeighbor->get_child(  17 );
                       tChildren[  58 ] = tNeighbor->get_child(  18 );
                       tChildren[  59 ] = tNeighbor->get_child(  19 );
                       tChildren[  60 ] = tNeighbor->get_child(  20 );
                       tChildren[  61 ] = tNeighbor->get_child(  21 );
                       tChildren[  62 ] = tNeighbor->get_child(  22 );
                       tChildren[  63 ] = tNeighbor->get_child(  23 );
                       tChildren[  64 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  35 );
                       tChildren[  76 ] = tNeighbor->get_child(  36 );
                       tChildren[  77 ] = tNeighbor->get_child(  37 );
                       tChildren[  78 ] = tNeighbor->get_child(  38 );
                       tChildren[  79 ] = tNeighbor->get_child(  39 );
                       tChildren[  80 ] = tNeighbor->get_child(  40 );
                       tChildren[  81 ] = tNeighbor->get_child(  41 );
                       tChildren[  82 ] = tNeighbor->get_child(  42 );
                       tChildren[  83 ] = tNeighbor->get_child(  43 );
                       tChildren[  84 ] = tNeighbor->get_child(  44 );
                       tChildren[  85 ] = tNeighbor->get_child(  45 );
                       tChildren[  86 ] = tNeighbor->get_child(  46 );
                       tChildren[  87 ] = tNeighbor->get_child(  47 );
                       tChildren[  88 ] = tNeighbor->get_child(  48 );
                       tChildren[  89 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  60 );
                       tChildren[ 101 ] = tNeighbor->get_child(  61 );
                       tChildren[ 102 ] = tNeighbor->get_child(  62 );
                       tChildren[ 103 ] = tNeighbor->get_child(  63 );
                       tChildren[ 104 ] = tNeighbor->get_child(  64 );
                       tChildren[ 105 ] = tNeighbor->get_child(  65 );
                       tChildren[ 106 ] = tNeighbor->get_child(  66 );
                       tChildren[ 107 ] = tNeighbor->get_child(  67 );
                       tChildren[ 108 ] = tNeighbor->get_child(  68 );
                       tChildren[ 109 ] = tNeighbor->get_child(  69 );
                       tChildren[ 110 ] = tNeighbor->get_child(  70 );
                       tChildren[ 111 ] = tNeighbor->get_child(  71 );
                       tChildren[ 112 ] = tNeighbor->get_child(  72 );
                       tChildren[ 113 ] = tNeighbor->get_child(  73 );
                       tChildren[ 114 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 15 exists
                if( tNeighbors[ 15 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 15 ];

                    // test if neighbor 15 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  52 ] = tNeighbor->get_child(   0 );
                       tChildren[  53 ] = tNeighbor->get_child(   1 );
                       tChildren[  54 ] = tNeighbor->get_child(   2 );
                       tChildren[  57 ] = tNeighbor->get_child(   5 );
                       tChildren[  58 ] = tNeighbor->get_child(   6 );
                       tChildren[  59 ] = tNeighbor->get_child(   7 );
                       tChildren[  62 ] = tNeighbor->get_child(  10 );
                       tChildren[  63 ] = tNeighbor->get_child(  11 );
                       tChildren[  64 ] = tNeighbor->get_child(  12 );
                       tChildren[  67 ] = tNeighbor->get_child(  15 );
                       tChildren[  68 ] = tNeighbor->get_child(  16 );
                       tChildren[  69 ] = tNeighbor->get_child(  17 );
                       tChildren[  72 ] = tNeighbor->get_child(  20 );
                       tChildren[  73 ] = tNeighbor->get_child(  21 );
                       tChildren[  74 ] = tNeighbor->get_child(  22 );
                       tChildren[  77 ] = tNeighbor->get_child(  25 );
                       tChildren[  78 ] = tNeighbor->get_child(  26 );
                       tChildren[  79 ] = tNeighbor->get_child(  27 );
                       tChildren[  82 ] = tNeighbor->get_child(  30 );
                       tChildren[  83 ] = tNeighbor->get_child(  31 );
                       tChildren[  84 ] = tNeighbor->get_child(  32 );
                       tChildren[  87 ] = tNeighbor->get_child(  35 );
                       tChildren[  88 ] = tNeighbor->get_child(  36 );
                       tChildren[  89 ] = tNeighbor->get_child(  37 );
                       tChildren[  92 ] = tNeighbor->get_child(  40 );
                       tChildren[  93 ] = tNeighbor->get_child(  41 );
                       tChildren[  94 ] = tNeighbor->get_child(  42 );
                       tChildren[  97 ] = tNeighbor->get_child(  45 );
                       tChildren[  98 ] = tNeighbor->get_child(  46 );
                       tChildren[  99 ] = tNeighbor->get_child(  47 );
                       tChildren[ 102 ] = tNeighbor->get_child(  50 );
                       tChildren[ 103 ] = tNeighbor->get_child(  51 );
                       tChildren[ 104 ] = tNeighbor->get_child(  52 );
                       tChildren[ 107 ] = tNeighbor->get_child(  55 );
                       tChildren[ 108 ] = tNeighbor->get_child(  56 );
                       tChildren[ 109 ] = tNeighbor->get_child(  57 );
                       tChildren[ 112 ] = tNeighbor->get_child(  60 );
                       tChildren[ 113 ] = tNeighbor->get_child(  61 );
                       tChildren[ 114 ] = tNeighbor->get_child(  62 );
                       tChildren[ 117 ] = tNeighbor->get_child(  65 );
                       tChildren[ 118 ] = tNeighbor->get_child(  66 );
                       tChildren[ 119 ] = tNeighbor->get_child(  67 );
                       tChildren[ 122 ] = tNeighbor->get_child(  70 );
                       tChildren[ 123 ] = tNeighbor->get_child(  71 );
                       tChildren[ 124 ] = tNeighbor->get_child(  72 );
                    }
                }

                // test if neighbor 16 exists
                if( tNeighbors[ 16 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 16 ];

                    // test if neighbor 16 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  60 ] = tNeighbor->get_child(   0 );
                       tChildren[  61 ] = tNeighbor->get_child(   1 );
                       tChildren[  62 ] = tNeighbor->get_child(   2 );
                       tChildren[  63 ] = tNeighbor->get_child(   3 );
                       tChildren[  64 ] = tNeighbor->get_child(   4 );
                       tChildren[  65 ] = tNeighbor->get_child(   5 );
                       tChildren[  66 ] = tNeighbor->get_child(   6 );
                       tChildren[  67 ] = tNeighbor->get_child(   7 );
                       tChildren[  68 ] = tNeighbor->get_child(   8 );
                       tChildren[  69 ] = tNeighbor->get_child(   9 );
                       tChildren[  70 ] = tNeighbor->get_child(  10 );
                       tChildren[  71 ] = tNeighbor->get_child(  11 );
                       tChildren[  72 ] = tNeighbor->get_child(  12 );
                       tChildren[  73 ] = tNeighbor->get_child(  13 );
                       tChildren[  74 ] = tNeighbor->get_child(  14 );
                       tChildren[  85 ] = tNeighbor->get_child(  25 );
                       tChildren[  86 ] = tNeighbor->get_child(  26 );
                       tChildren[  87 ] = tNeighbor->get_child(  27 );
                       tChildren[  88 ] = tNeighbor->get_child(  28 );
                       tChildren[  89 ] = tNeighbor->get_child(  29 );
                       tChildren[  90 ] = tNeighbor->get_child(  30 );
                       tChildren[  91 ] = tNeighbor->get_child(  31 );
                       tChildren[  92 ] = tNeighbor->get_child(  32 );
                       tChildren[  93 ] = tNeighbor->get_child(  33 );
                       tChildren[  94 ] = tNeighbor->get_child(  34 );
                       tChildren[  95 ] = tNeighbor->get_child(  35 );
                       tChildren[  96 ] = tNeighbor->get_child(  36 );
                       tChildren[  97 ] = tNeighbor->get_child(  37 );
                       tChildren[  98 ] = tNeighbor->get_child(  38 );
                       tChildren[  99 ] = tNeighbor->get_child(  39 );
                       tChildren[ 110 ] = tNeighbor->get_child(  50 );
                       tChildren[ 111 ] = tNeighbor->get_child(  51 );
                       tChildren[ 112 ] = tNeighbor->get_child(  52 );
                       tChildren[ 113 ] = tNeighbor->get_child(  53 );
                       tChildren[ 114 ] = tNeighbor->get_child(  54 );
                       tChildren[ 115 ] = tNeighbor->get_child(  55 );
                       tChildren[ 116 ] = tNeighbor->get_child(  56 );
                       tChildren[ 117 ] = tNeighbor->get_child(  57 );
                       tChildren[ 118 ] = tNeighbor->get_child(  58 );
                       tChildren[ 119 ] = tNeighbor->get_child(  59 );
                       tChildren[ 120 ] = tNeighbor->get_child(  60 );
                       tChildren[ 121 ] = tNeighbor->get_child(  61 );
                       tChildren[ 122 ] = tNeighbor->get_child(  62 );
                       tChildren[ 123 ] = tNeighbor->get_child(  63 );
                       tChildren[ 124 ] = tNeighbor->get_child(  64 );
                    }
                }

                // test if neighbor 17 exists
                if( tNeighbors[ 17 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 17 ];

                    // test if neighbor 17 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(   2 );
                       tChildren[  51 ] = tNeighbor->get_child(   3 );
                       tChildren[  52 ] = tNeighbor->get_child(   4 );
                       tChildren[  55 ] = tNeighbor->get_child(   7 );
                       tChildren[  56 ] = tNeighbor->get_child(   8 );
                       tChildren[  57 ] = tNeighbor->get_child(   9 );
                       tChildren[  60 ] = tNeighbor->get_child(  12 );
                       tChildren[  61 ] = tNeighbor->get_child(  13 );
                       tChildren[  62 ] = tNeighbor->get_child(  14 );
                       tChildren[  65 ] = tNeighbor->get_child(  17 );
                       tChildren[  66 ] = tNeighbor->get_child(  18 );
                       tChildren[  67 ] = tNeighbor->get_child(  19 );
                       tChildren[  70 ] = tNeighbor->get_child(  22 );
                       tChildren[  71 ] = tNeighbor->get_child(  23 );
                       tChildren[  72 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  27 );
                       tChildren[  76 ] = tNeighbor->get_child(  28 );
                       tChildren[  77 ] = tNeighbor->get_child(  29 );
                       tChildren[  80 ] = tNeighbor->get_child(  32 );
                       tChildren[  81 ] = tNeighbor->get_child(  33 );
                       tChildren[  82 ] = tNeighbor->get_child(  34 );
                       tChildren[  85 ] = tNeighbor->get_child(  37 );
                       tChildren[  86 ] = tNeighbor->get_child(  38 );
                       tChildren[  87 ] = tNeighbor->get_child(  39 );
                       tChildren[  90 ] = tNeighbor->get_child(  42 );
                       tChildren[  91 ] = tNeighbor->get_child(  43 );
                       tChildren[  92 ] = tNeighbor->get_child(  44 );
                       tChildren[  95 ] = tNeighbor->get_child(  47 );
                       tChildren[  96 ] = tNeighbor->get_child(  48 );
                       tChildren[  97 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  52 );
                       tChildren[ 101 ] = tNeighbor->get_child(  53 );
                       tChildren[ 102 ] = tNeighbor->get_child(  54 );
                       tChildren[ 105 ] = tNeighbor->get_child(  57 );
                       tChildren[ 106 ] = tNeighbor->get_child(  58 );
                       tChildren[ 107 ] = tNeighbor->get_child(  59 );
                       tChildren[ 110 ] = tNeighbor->get_child(  62 );
                       tChildren[ 111 ] = tNeighbor->get_child(  63 );
                       tChildren[ 112 ] = tNeighbor->get_child(  64 );
                       tChildren[ 115 ] = tNeighbor->get_child(  67 );
                       tChildren[ 116 ] = tNeighbor->get_child(  68 );
                       tChildren[ 117 ] = tNeighbor->get_child(  69 );
                       tChildren[ 120 ] = tNeighbor->get_child(  72 );
                       tChildren[ 121 ] = tNeighbor->get_child(  73 );
                       tChildren[ 122 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 18 exists
                if( tNeighbors[ 18 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 18 ];

                    // test if neighbor 18 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  62 );
                       tChildren[   1 ] = tNeighbor->get_child(  63 );
                       tChildren[   2 ] = tNeighbor->get_child(  64 );
                       tChildren[   5 ] = tNeighbor->get_child(  67 );
                       tChildren[   6 ] = tNeighbor->get_child(  68 );
                       tChildren[   7 ] = tNeighbor->get_child(  69 );
                       tChildren[  10 ] = tNeighbor->get_child(  72 );
                       tChildren[  11 ] = tNeighbor->get_child(  73 );
                       tChildren[  12 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  87 );
                       tChildren[  26 ] = tNeighbor->get_child(  88 );
                       tChildren[  27 ] = tNeighbor->get_child(  89 );
                       tChildren[  30 ] = tNeighbor->get_child(  92 );
                       tChildren[  31 ] = tNeighbor->get_child(  93 );
                       tChildren[  32 ] = tNeighbor->get_child(  94 );
                       tChildren[  35 ] = tNeighbor->get_child(  97 );
                       tChildren[  36 ] = tNeighbor->get_child(  98 );
                       tChildren[  37 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 112 );
                       tChildren[  51 ] = tNeighbor->get_child( 113 );
                       tChildren[  52 ] = tNeighbor->get_child( 114 );
                       tChildren[  55 ] = tNeighbor->get_child( 117 );
                       tChildren[  56 ] = tNeighbor->get_child( 118 );
                       tChildren[  57 ] = tNeighbor->get_child( 119 );
                       tChildren[  60 ] = tNeighbor->get_child( 122 );
                       tChildren[  61 ] = tNeighbor->get_child( 123 );
                       tChildren[  62 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 19 exists
                if( tNeighbors[ 19 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 19 ];

                    // test if neighbor 19 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  60 );
                       tChildren[   3 ] = tNeighbor->get_child(  61 );
                       tChildren[   4 ] = tNeighbor->get_child(  62 );
                       tChildren[   7 ] = tNeighbor->get_child(  65 );
                       tChildren[   8 ] = tNeighbor->get_child(  66 );
                       tChildren[   9 ] = tNeighbor->get_child(  67 );
                       tChildren[  12 ] = tNeighbor->get_child(  70 );
                       tChildren[  13 ] = tNeighbor->get_child(  71 );
                       tChildren[  14 ] = tNeighbor->get_child(  72 );
                       tChildren[  27 ] = tNeighbor->get_child(  85 );
                       tChildren[  28 ] = tNeighbor->get_child(  86 );
                       tChildren[  29 ] = tNeighbor->get_child(  87 );
                       tChildren[  32 ] = tNeighbor->get_child(  90 );
                       tChildren[  33 ] = tNeighbor->get_child(  91 );
                       tChildren[  34 ] = tNeighbor->get_child(  92 );
                       tChildren[  37 ] = tNeighbor->get_child(  95 );
                       tChildren[  38 ] = tNeighbor->get_child(  96 );
                       tChildren[  39 ] = tNeighbor->get_child(  97 );
                       tChildren[  52 ] = tNeighbor->get_child( 110 );
                       tChildren[  53 ] = tNeighbor->get_child( 111 );
                       tChildren[  54 ] = tNeighbor->get_child( 112 );
                       tChildren[  57 ] = tNeighbor->get_child( 115 );
                       tChildren[  58 ] = tNeighbor->get_child( 116 );
                       tChildren[  59 ] = tNeighbor->get_child( 117 );
                       tChildren[  62 ] = tNeighbor->get_child( 120 );
                       tChildren[  63 ] = tNeighbor->get_child( 121 );
                       tChildren[  64 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 20 exists
                if( tNeighbors[ 20 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 20 ];

                    // test if neighbor 20 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  12 ] = tNeighbor->get_child(  50 );
                       tChildren[  13 ] = tNeighbor->get_child(  51 );
                       tChildren[  14 ] = tNeighbor->get_child(  52 );
                       tChildren[  17 ] = tNeighbor->get_child(  55 );
                       tChildren[  18 ] = tNeighbor->get_child(  56 );
                       tChildren[  19 ] = tNeighbor->get_child(  57 );
                       tChildren[  22 ] = tNeighbor->get_child(  60 );
                       tChildren[  23 ] = tNeighbor->get_child(  61 );
                       tChildren[  24 ] = tNeighbor->get_child(  62 );
                       tChildren[  37 ] = tNeighbor->get_child(  75 );
                       tChildren[  38 ] = tNeighbor->get_child(  76 );
                       tChildren[  39 ] = tNeighbor->get_child(  77 );
                       tChildren[  42 ] = tNeighbor->get_child(  80 );
                       tChildren[  43 ] = tNeighbor->get_child(  81 );
                       tChildren[  44 ] = tNeighbor->get_child(  82 );
                       tChildren[  47 ] = tNeighbor->get_child(  85 );
                       tChildren[  48 ] = tNeighbor->get_child(  86 );
                       tChildren[  49 ] = tNeighbor->get_child(  87 );
                       tChildren[  62 ] = tNeighbor->get_child( 100 );
                       tChildren[  63 ] = tNeighbor->get_child( 101 );
                       tChildren[  64 ] = tNeighbor->get_child( 102 );
                       tChildren[  67 ] = tNeighbor->get_child( 105 );
                       tChildren[  68 ] = tNeighbor->get_child( 106 );
                       tChildren[  69 ] = tNeighbor->get_child( 107 );
                       tChildren[  72 ] = tNeighbor->get_child( 110 );
                       tChildren[  73 ] = tNeighbor->get_child( 111 );
                       tChildren[  74 ] = tNeighbor->get_child( 112 );
                    }
                }

                // test if neighbor 21 exists
                if( tNeighbors[ 21 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 21 ];

                    // test if neighbor 21 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(  52 );
                       tChildren[  11 ] = tNeighbor->get_child(  53 );
                       tChildren[  12 ] = tNeighbor->get_child(  54 );
                       tChildren[  15 ] = tNeighbor->get_child(  57 );
                       tChildren[  16 ] = tNeighbor->get_child(  58 );
                       tChildren[  17 ] = tNeighbor->get_child(  59 );
                       tChildren[  20 ] = tNeighbor->get_child(  62 );
                       tChildren[  21 ] = tNeighbor->get_child(  63 );
                       tChildren[  22 ] = tNeighbor->get_child(  64 );
                       tChildren[  35 ] = tNeighbor->get_child(  77 );
                       tChildren[  36 ] = tNeighbor->get_child(  78 );
                       tChildren[  37 ] = tNeighbor->get_child(  79 );
                       tChildren[  40 ] = tNeighbor->get_child(  82 );
                       tChildren[  41 ] = tNeighbor->get_child(  83 );
                       tChildren[  42 ] = tNeighbor->get_child(  84 );
                       tChildren[  45 ] = tNeighbor->get_child(  87 );
                       tChildren[  46 ] = tNeighbor->get_child(  88 );
                       tChildren[  47 ] = tNeighbor->get_child(  89 );
                       tChildren[  60 ] = tNeighbor->get_child( 102 );
                       tChildren[  61 ] = tNeighbor->get_child( 103 );
                       tChildren[  62 ] = tNeighbor->get_child( 104 );
                       tChildren[  65 ] = tNeighbor->get_child( 107 );
                       tChildren[  66 ] = tNeighbor->get_child( 108 );
                       tChildren[  67 ] = tNeighbor->get_child( 109 );
                       tChildren[  70 ] = tNeighbor->get_child( 112 );
                       tChildren[  71 ] = tNeighbor->get_child( 113 );
                       tChildren[  72 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 22 exists
                if( tNeighbors[ 22 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 22 ];

                    // test if neighbor 22 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  12 );
                       tChildren[  51 ] = tNeighbor->get_child(  13 );
                       tChildren[  52 ] = tNeighbor->get_child(  14 );
                       tChildren[  55 ] = tNeighbor->get_child(  17 );
                       tChildren[  56 ] = tNeighbor->get_child(  18 );
                       tChildren[  57 ] = tNeighbor->get_child(  19 );
                       tChildren[  60 ] = tNeighbor->get_child(  22 );
                       tChildren[  61 ] = tNeighbor->get_child(  23 );
                       tChildren[  62 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  37 );
                       tChildren[  76 ] = tNeighbor->get_child(  38 );
                       tChildren[  77 ] = tNeighbor->get_child(  39 );
                       tChildren[  80 ] = tNeighbor->get_child(  42 );
                       tChildren[  81 ] = tNeighbor->get_child(  43 );
                       tChildren[  82 ] = tNeighbor->get_child(  44 );
                       tChildren[  85 ] = tNeighbor->get_child(  47 );
                       tChildren[  86 ] = tNeighbor->get_child(  48 );
                       tChildren[  87 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  62 );
                       tChildren[ 101 ] = tNeighbor->get_child(  63 );
                       tChildren[ 102 ] = tNeighbor->get_child(  64 );
                       tChildren[ 105 ] = tNeighbor->get_child(  67 );
                       tChildren[ 106 ] = tNeighbor->get_child(  68 );
                       tChildren[ 107 ] = tNeighbor->get_child(  69 );
                       tChildren[ 110 ] = tNeighbor->get_child(  72 );
                       tChildren[ 111 ] = tNeighbor->get_child(  73 );
                       tChildren[ 112 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 23 exists
                if( tNeighbors[ 23 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 23 ];

                    // test if neighbor 23 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  52 ] = tNeighbor->get_child(  10 );
                       tChildren[  53 ] = tNeighbor->get_child(  11 );
                       tChildren[  54 ] = tNeighbor->get_child(  12 );
                       tChildren[  57 ] = tNeighbor->get_child(  15 );
                       tChildren[  58 ] = tNeighbor->get_child(  16 );
                       tChildren[  59 ] = tNeighbor->get_child(  17 );
                       tChildren[  62 ] = tNeighbor->get_child(  20 );
                       tChildren[  63 ] = tNeighbor->get_child(  21 );
                       tChildren[  64 ] = tNeighbor->get_child(  22 );
                       tChildren[  77 ] = tNeighbor->get_child(  35 );
                       tChildren[  78 ] = tNeighbor->get_child(  36 );
                       tChildren[  79 ] = tNeighbor->get_child(  37 );
                       tChildren[  82 ] = tNeighbor->get_child(  40 );
                       tChildren[  83 ] = tNeighbor->get_child(  41 );
                       tChildren[  84 ] = tNeighbor->get_child(  42 );
                       tChildren[  87 ] = tNeighbor->get_child(  45 );
                       tChildren[  88 ] = tNeighbor->get_child(  46 );
                       tChildren[  89 ] = tNeighbor->get_child(  47 );
                       tChildren[ 102 ] = tNeighbor->get_child(  60 );
                       tChildren[ 103 ] = tNeighbor->get_child(  61 );
                       tChildren[ 104 ] = tNeighbor->get_child(  62 );
                       tChildren[ 107 ] = tNeighbor->get_child(  65 );
                       tChildren[ 108 ] = tNeighbor->get_child(  66 );
                       tChildren[ 109 ] = tNeighbor->get_child(  67 );
                       tChildren[ 112 ] = tNeighbor->get_child(  70 );
                       tChildren[ 113 ] = tNeighbor->get_child(  71 );
                       tChildren[ 114 ] = tNeighbor->get_child(  72 );
                    }
                }

                // test if neighbor 24 exists
                if( tNeighbors[ 24 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 24 ];

                    // test if neighbor 24 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  62 ] = tNeighbor->get_child(   0 );
                       tChildren[  63 ] = tNeighbor->get_child(   1 );
                       tChildren[  64 ] = tNeighbor->get_child(   2 );
                       tChildren[  67 ] = tNeighbor->get_child(   5 );
                       tChildren[  68 ] = tNeighbor->get_child(   6 );
                       tChildren[  69 ] = tNeighbor->get_child(   7 );
                       tChildren[  72 ] = tNeighbor->get_child(  10 );
                       tChildren[  73 ] = tNeighbor->get_child(  11 );
                       tChildren[  74 ] = tNeighbor->get_child(  12 );
                       tChildren[  87 ] = tNeighbor->get_child(  25 );
                       tChildren[  88 ] = tNeighbor->get_child(  26 );
                       tChildren[  89 ] = tNeighbor->get_child(  27 );
                       tChildren[  92 ] = tNeighbor->get_child(  30 );
                       tChildren[  93 ] = tNeighbor->get_child(  31 );
                       tChildren[  94 ] = tNeighbor->get_child(  32 );
                       tChildren[  97 ] = tNeighbor->get_child(  35 );
                       tChildren[  98 ] = tNeighbor->get_child(  36 );
                       tChildren[  99 ] = tNeighbor->get_child(  37 );
                       tChildren[ 112 ] = tNeighbor->get_child(  50 );
                       tChildren[ 113 ] = tNeighbor->get_child(  51 );
                       tChildren[ 114 ] = tNeighbor->get_child(  52 );
                       tChildren[ 117 ] = tNeighbor->get_child(  55 );
                       tChildren[ 118 ] = tNeighbor->get_child(  56 );
                       tChildren[ 119 ] = tNeighbor->get_child(  57 );
                       tChildren[ 122 ] = tNeighbor->get_child(  60 );
                       tChildren[ 123 ] = tNeighbor->get_child(  61 );
                       tChildren[ 124 ] = tNeighbor->get_child(  62 );
                    }
                }

                // test if neighbor 25 exists
                if( tNeighbors[ 25 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 25 ];

                    // test if neighbor 25 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  60 ] = tNeighbor->get_child(   2 );
                       tChildren[  61 ] = tNeighbor->get_child(   3 );
                       tChildren[  62 ] = tNeighbor->get_child(   4 );
                       tChildren[  65 ] = tNeighbor->get_child(   7 );
                       tChildren[  66 ] = tNeighbor->get_child(   8 );
                       tChildren[  67 ] = tNeighbor->get_child(   9 );
                       tChildren[  70 ] = tNeighbor->get_child(  12 );
                       tChildren[  71 ] = tNeighbor->get_child(  13 );
                       tChildren[  72 ] = tNeighbor->get_child(  14 );
                       tChildren[  85 ] = tNeighbor->get_child(  27 );
                       tChildren[  86 ] = tNeighbor->get_child(  28 );
                       tChildren[  87 ] = tNeighbor->get_child(  29 );
                       tChildren[  90 ] = tNeighbor->get_child(  32 );
                       tChildren[  91 ] = tNeighbor->get_child(  33 );
                       tChildren[  92 ] = tNeighbor->get_child(  34 );
                       tChildren[  95 ] = tNeighbor->get_child(  37 );
                       tChildren[  96 ] = tNeighbor->get_child(  38 );
                       tChildren[  97 ] = tNeighbor->get_child(  39 );
                       tChildren[ 110 ] = tNeighbor->get_child(  52 );
                       tChildren[ 111 ] = tNeighbor->get_child(  53 );
                       tChildren[ 112 ] = tNeighbor->get_child(  54 );
                       tChildren[ 115 ] = tNeighbor->get_child(  57 );
                       tChildren[ 116 ] = tNeighbor->get_child(  58 );
                       tChildren[ 117 ] = tNeighbor->get_child(  59 );
                       tChildren[ 120 ] = tNeighbor->get_child(  62 );
                       tChildren[ 121 ] = tNeighbor->get_child(  63 );
                       tChildren[ 122 ] = tNeighbor->get_child(  64 );
                    }
                }

                // test if neighbor 26 exists
                if( tNeighbors[ 26 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 26 ];

                    // test if neighbor 26 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 27 exists
                if( tNeighbors[ 27 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 27 ];

                    // test if neighbor 27 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 122 );
                       tChildren[   1 ] = tNeighbor->get_child( 123 );
                       tChildren[   2 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 28 exists
                if( tNeighbors[ 28 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 28 ];

                    // test if neighbor 28 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 120 );
                       tChildren[   1 ] = tNeighbor->get_child( 121 );
                       tChildren[   2 ] = tNeighbor->get_child( 122 );
                       tChildren[   3 ] = tNeighbor->get_child( 123 );
                       tChildren[   4 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 29 exists
                if( tNeighbors[ 29 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 29 ];

                    // test if neighbor 29 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child( 120 );
                       tChildren[   3 ] = tNeighbor->get_child( 121 );
                       tChildren[   4 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 30 exists
                if( tNeighbors[ 30 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 30 ];

                    // test if neighbor 30 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 31 exists
                if( tNeighbors[ 31 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 31 ];

                    // test if neighbor 31 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 114 );
                       tChildren[   5 ] = tNeighbor->get_child( 119 );
                       tChildren[  10 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 32 exists
                if( tNeighbors[ 32 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 32 ];

                    // test if neighbor 32 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 112 );
                       tChildren[   1 ] = tNeighbor->get_child( 113 );
                       tChildren[   2 ] = tNeighbor->get_child( 114 );
                       tChildren[   5 ] = tNeighbor->get_child( 117 );
                       tChildren[   6 ] = tNeighbor->get_child( 118 );
                       tChildren[   7 ] = tNeighbor->get_child( 119 );
                       tChildren[  10 ] = tNeighbor->get_child( 122 );
                       tChildren[  11 ] = tNeighbor->get_child( 123 );
                       tChildren[  12 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 33 exists
                if( tNeighbors[ 33 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 33 ];

                    // test if neighbor 33 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 110 );
                       tChildren[   1 ] = tNeighbor->get_child( 111 );
                       tChildren[   2 ] = tNeighbor->get_child( 112 );
                       tChildren[   3 ] = tNeighbor->get_child( 113 );
                       tChildren[   4 ] = tNeighbor->get_child( 114 );
                       tChildren[   5 ] = tNeighbor->get_child( 115 );
                       tChildren[   6 ] = tNeighbor->get_child( 116 );
                       tChildren[   7 ] = tNeighbor->get_child( 117 );
                       tChildren[   8 ] = tNeighbor->get_child( 118 );
                       tChildren[   9 ] = tNeighbor->get_child( 119 );
                       tChildren[  10 ] = tNeighbor->get_child( 120 );
                       tChildren[  11 ] = tNeighbor->get_child( 121 );
                       tChildren[  12 ] = tNeighbor->get_child( 122 );
                       tChildren[  13 ] = tNeighbor->get_child( 123 );
                       tChildren[  14 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 34 exists
                if( tNeighbors[ 34 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 34 ];

                    // test if neighbor 34 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child( 110 );
                       tChildren[   3 ] = tNeighbor->get_child( 111 );
                       tChildren[   4 ] = tNeighbor->get_child( 112 );
                       tChildren[   7 ] = tNeighbor->get_child( 115 );
                       tChildren[   8 ] = tNeighbor->get_child( 116 );
                       tChildren[   9 ] = tNeighbor->get_child( 117 );
                       tChildren[  12 ] = tNeighbor->get_child( 120 );
                       tChildren[  13 ] = tNeighbor->get_child( 121 );
                       tChildren[  14 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 35 exists
                if( tNeighbors[ 35 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 35 ];

                    // test if neighbor 35 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child( 110 );
                       tChildren[   9 ] = tNeighbor->get_child( 115 );
                       tChildren[  14 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 36 exists
                if( tNeighbors[ 36 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 36 ];

                    // test if neighbor 36 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 104 );
                       tChildren[   5 ] = tNeighbor->get_child( 109 );
                       tChildren[  10 ] = tNeighbor->get_child( 114 );
                       tChildren[  15 ] = tNeighbor->get_child( 119 );
                       tChildren[  20 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 37 exists
                if( tNeighbors[ 37 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 37 ];

                    // test if neighbor 37 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 102 );
                       tChildren[   1 ] = tNeighbor->get_child( 103 );
                       tChildren[   2 ] = tNeighbor->get_child( 104 );
                       tChildren[   5 ] = tNeighbor->get_child( 107 );
                       tChildren[   6 ] = tNeighbor->get_child( 108 );
                       tChildren[   7 ] = tNeighbor->get_child( 109 );
                       tChildren[  10 ] = tNeighbor->get_child( 112 );
                       tChildren[  11 ] = tNeighbor->get_child( 113 );
                       tChildren[  12 ] = tNeighbor->get_child( 114 );
                       tChildren[  15 ] = tNeighbor->get_child( 117 );
                       tChildren[  16 ] = tNeighbor->get_child( 118 );
                       tChildren[  17 ] = tNeighbor->get_child( 119 );
                       tChildren[  20 ] = tNeighbor->get_child( 122 );
                       tChildren[  21 ] = tNeighbor->get_child( 123 );
                       tChildren[  22 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 38 exists
                if( tNeighbors[ 38 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 38 ];

                    // test if neighbor 38 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child( 100 );
                       tChildren[   1 ] = tNeighbor->get_child( 101 );
                       tChildren[   2 ] = tNeighbor->get_child( 102 );
                       tChildren[   3 ] = tNeighbor->get_child( 103 );
                       tChildren[   4 ] = tNeighbor->get_child( 104 );
                       tChildren[   5 ] = tNeighbor->get_child( 105 );
                       tChildren[   6 ] = tNeighbor->get_child( 106 );
                       tChildren[   7 ] = tNeighbor->get_child( 107 );
                       tChildren[   8 ] = tNeighbor->get_child( 108 );
                       tChildren[   9 ] = tNeighbor->get_child( 109 );
                       tChildren[  10 ] = tNeighbor->get_child( 110 );
                       tChildren[  11 ] = tNeighbor->get_child( 111 );
                       tChildren[  12 ] = tNeighbor->get_child( 112 );
                       tChildren[  13 ] = tNeighbor->get_child( 113 );
                       tChildren[  14 ] = tNeighbor->get_child( 114 );
                       tChildren[  15 ] = tNeighbor->get_child( 115 );
                       tChildren[  16 ] = tNeighbor->get_child( 116 );
                       tChildren[  17 ] = tNeighbor->get_child( 117 );
                       tChildren[  18 ] = tNeighbor->get_child( 118 );
                       tChildren[  19 ] = tNeighbor->get_child( 119 );
                       tChildren[  20 ] = tNeighbor->get_child( 120 );
                       tChildren[  21 ] = tNeighbor->get_child( 121 );
                       tChildren[  22 ] = tNeighbor->get_child( 122 );
                       tChildren[  23 ] = tNeighbor->get_child( 123 );
                       tChildren[  24 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 39 exists
                if( tNeighbors[ 39 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 39 ];

                    // test if neighbor 39 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child( 100 );
                       tChildren[   3 ] = tNeighbor->get_child( 101 );
                       tChildren[   4 ] = tNeighbor->get_child( 102 );
                       tChildren[   7 ] = tNeighbor->get_child( 105 );
                       tChildren[   8 ] = tNeighbor->get_child( 106 );
                       tChildren[   9 ] = tNeighbor->get_child( 107 );
                       tChildren[  12 ] = tNeighbor->get_child( 110 );
                       tChildren[  13 ] = tNeighbor->get_child( 111 );
                       tChildren[  14 ] = tNeighbor->get_child( 112 );
                       tChildren[  17 ] = tNeighbor->get_child( 115 );
                       tChildren[  18 ] = tNeighbor->get_child( 116 );
                       tChildren[  19 ] = tNeighbor->get_child( 117 );
                       tChildren[  22 ] = tNeighbor->get_child( 120 );
                       tChildren[  23 ] = tNeighbor->get_child( 121 );
                       tChildren[  24 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 40 exists
                if( tNeighbors[ 40 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 40 ];

                    // test if neighbor 40 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child( 100 );
                       tChildren[   9 ] = tNeighbor->get_child( 105 );
                       tChildren[  14 ] = tNeighbor->get_child( 110 );
                       tChildren[  19 ] = tNeighbor->get_child( 115 );
                       tChildren[  24 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 41 exists
                if( tNeighbors[ 41 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 41 ];

                    // test if neighbor 41 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child( 104 );
                       tChildren[  15 ] = tNeighbor->get_child( 109 );
                       tChildren[  20 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 42 exists
                if( tNeighbors[ 42 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 42 ];

                    // test if neighbor 42 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child( 102 );
                       tChildren[  11 ] = tNeighbor->get_child( 103 );
                       tChildren[  12 ] = tNeighbor->get_child( 104 );
                       tChildren[  15 ] = tNeighbor->get_child( 107 );
                       tChildren[  16 ] = tNeighbor->get_child( 108 );
                       tChildren[  17 ] = tNeighbor->get_child( 109 );
                       tChildren[  20 ] = tNeighbor->get_child( 112 );
                       tChildren[  21 ] = tNeighbor->get_child( 113 );
                       tChildren[  22 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 43 exists
                if( tNeighbors[ 43 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 43 ];

                    // test if neighbor 43 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child( 100 );
                       tChildren[  11 ] = tNeighbor->get_child( 101 );
                       tChildren[  12 ] = tNeighbor->get_child( 102 );
                       tChildren[  13 ] = tNeighbor->get_child( 103 );
                       tChildren[  14 ] = tNeighbor->get_child( 104 );
                       tChildren[  15 ] = tNeighbor->get_child( 105 );
                       tChildren[  16 ] = tNeighbor->get_child( 106 );
                       tChildren[  17 ] = tNeighbor->get_child( 107 );
                       tChildren[  18 ] = tNeighbor->get_child( 108 );
                       tChildren[  19 ] = tNeighbor->get_child( 109 );
                       tChildren[  20 ] = tNeighbor->get_child( 110 );
                       tChildren[  21 ] = tNeighbor->get_child( 111 );
                       tChildren[  22 ] = tNeighbor->get_child( 112 );
                       tChildren[  23 ] = tNeighbor->get_child( 113 );
                       tChildren[  24 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 44 exists
                if( tNeighbors[ 44 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 44 ];

                    // test if neighbor 44 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  12 ] = tNeighbor->get_child( 100 );
                       tChildren[  13 ] = tNeighbor->get_child( 101 );
                       tChildren[  14 ] = tNeighbor->get_child( 102 );
                       tChildren[  17 ] = tNeighbor->get_child( 105 );
                       tChildren[  18 ] = tNeighbor->get_child( 106 );
                       tChildren[  19 ] = tNeighbor->get_child( 107 );
                       tChildren[  22 ] = tNeighbor->get_child( 110 );
                       tChildren[  23 ] = tNeighbor->get_child( 111 );
                       tChildren[  24 ] = tNeighbor->get_child( 112 );
                    }
                }

                // test if neighbor 45 exists
                if( tNeighbors[ 45 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 45 ];

                    // test if neighbor 45 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  14 ] = tNeighbor->get_child( 100 );
                       tChildren[  19 ] = tNeighbor->get_child( 105 );
                       tChildren[  24 ] = tNeighbor->get_child( 110 );
                    }
                }

                // test if neighbor 46 exists
                if( tNeighbors[ 46 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 46 ];

                    // test if neighbor 46 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 47 exists
                if( tNeighbors[ 47 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 47 ];

                    // test if neighbor 47 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child( 102 );
                       tChildren[  21 ] = tNeighbor->get_child( 103 );
                       tChildren[  22 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 48 exists
                if( tNeighbors[ 48 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 48 ];

                    // test if neighbor 48 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child( 100 );
                       tChildren[  21 ] = tNeighbor->get_child( 101 );
                       tChildren[  22 ] = tNeighbor->get_child( 102 );
                       tChildren[  23 ] = tNeighbor->get_child( 103 );
                       tChildren[  24 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 49 exists
                if( tNeighbors[ 49 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 49 ];

                    // test if neighbor 49 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  22 ] = tNeighbor->get_child( 100 );
                       tChildren[  23 ] = tNeighbor->get_child( 101 );
                       tChildren[  24 ] = tNeighbor->get_child( 102 );
                    }
                }

                // test if neighbor 50 exists
                if( tNeighbors[ 50 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 50 ];

                    // test if neighbor 50 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  24 ] = tNeighbor->get_child( 100 );
                    }
                }

                // test if neighbor 51 exists
                if( tNeighbors[ 51 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 51 ];

                    // test if neighbor 51 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 52 exists
                if( tNeighbors[ 52 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 52 ];

                    // test if neighbor 52 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  72 );
                       tChildren[   1 ] = tNeighbor->get_child(  73 );
                       tChildren[   2 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  97 );
                       tChildren[  26 ] = tNeighbor->get_child(  98 );
                       tChildren[  27 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 122 );
                       tChildren[  51 ] = tNeighbor->get_child( 123 );
                       tChildren[  52 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 53 exists
                if( tNeighbors[ 53 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 53 ];

                    // test if neighbor 53 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  70 );
                       tChildren[   1 ] = tNeighbor->get_child(  71 );
                       tChildren[   2 ] = tNeighbor->get_child(  72 );
                       tChildren[   3 ] = tNeighbor->get_child(  73 );
                       tChildren[   4 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  95 );
                       tChildren[  26 ] = tNeighbor->get_child(  96 );
                       tChildren[  27 ] = tNeighbor->get_child(  97 );
                       tChildren[  28 ] = tNeighbor->get_child(  98 );
                       tChildren[  29 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 120 );
                       tChildren[  51 ] = tNeighbor->get_child( 121 );
                       tChildren[  52 ] = tNeighbor->get_child( 122 );
                       tChildren[  53 ] = tNeighbor->get_child( 123 );
                       tChildren[  54 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 54 exists
                if( tNeighbors[ 54 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 54 ];

                    // test if neighbor 54 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  70 );
                       tChildren[   3 ] = tNeighbor->get_child(  71 );
                       tChildren[   4 ] = tNeighbor->get_child(  72 );
                       tChildren[  27 ] = tNeighbor->get_child(  95 );
                       tChildren[  28 ] = tNeighbor->get_child(  96 );
                       tChildren[  29 ] = tNeighbor->get_child(  97 );
                       tChildren[  52 ] = tNeighbor->get_child( 120 );
                       tChildren[  53 ] = tNeighbor->get_child( 121 );
                       tChildren[  54 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 55 exists
                if( tNeighbors[ 55 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 55 ];

                    // test if neighbor 55 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  70 );
                       tChildren[  29 ] = tNeighbor->get_child(  95 );
                       tChildren[  54 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 56 exists
                if( tNeighbors[ 56 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 56 ];

                    // test if neighbor 56 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  64 );
                       tChildren[   5 ] = tNeighbor->get_child(  69 );
                       tChildren[  10 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  89 );
                       tChildren[  30 ] = tNeighbor->get_child(  94 );
                       tChildren[  35 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 114 );
                       tChildren[  55 ] = tNeighbor->get_child( 119 );
                       tChildren[  60 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 57 exists
                if( tNeighbors[ 57 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 57 ];

                    // test if neighbor 57 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  60 );
                       tChildren[   9 ] = tNeighbor->get_child(  65 );
                       tChildren[  14 ] = tNeighbor->get_child(  70 );
                       tChildren[  29 ] = tNeighbor->get_child(  85 );
                       tChildren[  34 ] = tNeighbor->get_child(  90 );
                       tChildren[  39 ] = tNeighbor->get_child(  95 );
                       tChildren[  54 ] = tNeighbor->get_child( 110 );
                       tChildren[  59 ] = tNeighbor->get_child( 115 );
                       tChildren[  64 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 58 exists
                if( tNeighbors[ 58 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 58 ];

                    // test if neighbor 58 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  54 );
                       tChildren[   5 ] = tNeighbor->get_child(  59 );
                       tChildren[  10 ] = tNeighbor->get_child(  64 );
                       tChildren[  15 ] = tNeighbor->get_child(  69 );
                       tChildren[  20 ] = tNeighbor->get_child(  74 );
                       tChildren[  25 ] = tNeighbor->get_child(  79 );
                       tChildren[  30 ] = tNeighbor->get_child(  84 );
                       tChildren[  35 ] = tNeighbor->get_child(  89 );
                       tChildren[  40 ] = tNeighbor->get_child(  94 );
                       tChildren[  45 ] = tNeighbor->get_child(  99 );
                       tChildren[  50 ] = tNeighbor->get_child( 104 );
                       tChildren[  55 ] = tNeighbor->get_child( 109 );
                       tChildren[  60 ] = tNeighbor->get_child( 114 );
                       tChildren[  65 ] = tNeighbor->get_child( 119 );
                       tChildren[  70 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 59 exists
                if( tNeighbors[ 59 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 59 ];

                    // test if neighbor 59 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  50 );
                       tChildren[   9 ] = tNeighbor->get_child(  55 );
                       tChildren[  14 ] = tNeighbor->get_child(  60 );
                       tChildren[  19 ] = tNeighbor->get_child(  65 );
                       tChildren[  24 ] = tNeighbor->get_child(  70 );
                       tChildren[  29 ] = tNeighbor->get_child(  75 );
                       tChildren[  34 ] = tNeighbor->get_child(  80 );
                       tChildren[  39 ] = tNeighbor->get_child(  85 );
                       tChildren[  44 ] = tNeighbor->get_child(  90 );
                       tChildren[  49 ] = tNeighbor->get_child(  95 );
                       tChildren[  54 ] = tNeighbor->get_child( 100 );
                       tChildren[  59 ] = tNeighbor->get_child( 105 );
                       tChildren[  64 ] = tNeighbor->get_child( 110 );
                       tChildren[  69 ] = tNeighbor->get_child( 115 );
                       tChildren[  74 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 60 exists
                if( tNeighbors[ 60 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 60 ];

                    // test if neighbor 60 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(  54 );
                       tChildren[  15 ] = tNeighbor->get_child(  59 );
                       tChildren[  20 ] = tNeighbor->get_child(  64 );
                       tChildren[  35 ] = tNeighbor->get_child(  79 );
                       tChildren[  40 ] = tNeighbor->get_child(  84 );
                       tChildren[  45 ] = tNeighbor->get_child(  89 );
                       tChildren[  60 ] = tNeighbor->get_child( 104 );
                       tChildren[  65 ] = tNeighbor->get_child( 109 );
                       tChildren[  70 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 61 exists
                if( tNeighbors[ 61 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 61 ];

                    // test if neighbor 61 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  14 ] = tNeighbor->get_child(  50 );
                       tChildren[  19 ] = tNeighbor->get_child(  55 );
                       tChildren[  24 ] = tNeighbor->get_child(  60 );
                       tChildren[  39 ] = tNeighbor->get_child(  75 );
                       tChildren[  44 ] = tNeighbor->get_child(  80 );
                       tChildren[  49 ] = tNeighbor->get_child(  85 );
                       tChildren[  64 ] = tNeighbor->get_child( 100 );
                       tChildren[  69 ] = tNeighbor->get_child( 105 );
                       tChildren[  74 ] = tNeighbor->get_child( 110 );
                    }
                }

                // test if neighbor 62 exists
                if( tNeighbors[ 62 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 62 ];

                    // test if neighbor 62 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(  54 );
                       tChildren[  45 ] = tNeighbor->get_child(  79 );
                       tChildren[  70 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 63 exists
                if( tNeighbors[ 63 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 63 ];

                    // test if neighbor 63 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(  52 );
                       tChildren[  21 ] = tNeighbor->get_child(  53 );
                       tChildren[  22 ] = tNeighbor->get_child(  54 );
                       tChildren[  45 ] = tNeighbor->get_child(  77 );
                       tChildren[  46 ] = tNeighbor->get_child(  78 );
                       tChildren[  47 ] = tNeighbor->get_child(  79 );
                       tChildren[  70 ] = tNeighbor->get_child( 102 );
                       tChildren[  71 ] = tNeighbor->get_child( 103 );
                       tChildren[  72 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 64 exists
                if( tNeighbors[ 64 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 64 ];

                    // test if neighbor 64 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(  50 );
                       tChildren[  21 ] = tNeighbor->get_child(  51 );
                       tChildren[  22 ] = tNeighbor->get_child(  52 );
                       tChildren[  23 ] = tNeighbor->get_child(  53 );
                       tChildren[  24 ] = tNeighbor->get_child(  54 );
                       tChildren[  45 ] = tNeighbor->get_child(  75 );
                       tChildren[  46 ] = tNeighbor->get_child(  76 );
                       tChildren[  47 ] = tNeighbor->get_child(  77 );
                       tChildren[  48 ] = tNeighbor->get_child(  78 );
                       tChildren[  49 ] = tNeighbor->get_child(  79 );
                       tChildren[  70 ] = tNeighbor->get_child( 100 );
                       tChildren[  71 ] = tNeighbor->get_child( 101 );
                       tChildren[  72 ] = tNeighbor->get_child( 102 );
                       tChildren[  73 ] = tNeighbor->get_child( 103 );
                       tChildren[  74 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 65 exists
                if( tNeighbors[ 65 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 65 ];

                    // test if neighbor 65 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  22 ] = tNeighbor->get_child(  50 );
                       tChildren[  23 ] = tNeighbor->get_child(  51 );
                       tChildren[  24 ] = tNeighbor->get_child(  52 );
                       tChildren[  47 ] = tNeighbor->get_child(  75 );
                       tChildren[  48 ] = tNeighbor->get_child(  76 );
                       tChildren[  49 ] = tNeighbor->get_child(  77 );
                       tChildren[  72 ] = tNeighbor->get_child( 100 );
                       tChildren[  73 ] = tNeighbor->get_child( 101 );
                       tChildren[  74 ] = tNeighbor->get_child( 102 );
                    }
                }

                // test if neighbor 66 exists
                if( tNeighbors[ 66 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 66 ];

                    // test if neighbor 66 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  24 ] = tNeighbor->get_child(  50 );
                       tChildren[  49 ] = tNeighbor->get_child(  75 );
                       tChildren[  74 ] = tNeighbor->get_child( 100 );
                    }
                }

                // test if neighbor 67 exists
                if( tNeighbors[ 67 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 67 ];

                    // test if neighbor 67 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  24 );
                       tChildren[  25 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 68 exists
                if( tNeighbors[ 68 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 68 ];

                    // test if neighbor 68 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  22 );
                       tChildren[   1 ] = tNeighbor->get_child(  23 );
                       tChildren[   2 ] = tNeighbor->get_child(  24 );
                       tChildren[  25 ] = tNeighbor->get_child(  47 );
                       tChildren[  26 ] = tNeighbor->get_child(  48 );
                       tChildren[  27 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  72 );
                       tChildren[  51 ] = tNeighbor->get_child(  73 );
                       tChildren[  52 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  97 );
                       tChildren[  76 ] = tNeighbor->get_child(  98 );
                       tChildren[  77 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 122 );
                       tChildren[ 101 ] = tNeighbor->get_child( 123 );
                       tChildren[ 102 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 69 exists
                if( tNeighbors[ 69 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 69 ];

                    // test if neighbor 69 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  20 );
                       tChildren[   1 ] = tNeighbor->get_child(  21 );
                       tChildren[   2 ] = tNeighbor->get_child(  22 );
                       tChildren[   3 ] = tNeighbor->get_child(  23 );
                       tChildren[   4 ] = tNeighbor->get_child(  24 );
                       tChildren[  25 ] = tNeighbor->get_child(  45 );
                       tChildren[  26 ] = tNeighbor->get_child(  46 );
                       tChildren[  27 ] = tNeighbor->get_child(  47 );
                       tChildren[  28 ] = tNeighbor->get_child(  48 );
                       tChildren[  29 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  70 );
                       tChildren[  51 ] = tNeighbor->get_child(  71 );
                       tChildren[  52 ] = tNeighbor->get_child(  72 );
                       tChildren[  53 ] = tNeighbor->get_child(  73 );
                       tChildren[  54 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  95 );
                       tChildren[  76 ] = tNeighbor->get_child(  96 );
                       tChildren[  77 ] = tNeighbor->get_child(  97 );
                       tChildren[  78 ] = tNeighbor->get_child(  98 );
                       tChildren[  79 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 120 );
                       tChildren[ 101 ] = tNeighbor->get_child( 121 );
                       tChildren[ 102 ] = tNeighbor->get_child( 122 );
                       tChildren[ 103 ] = tNeighbor->get_child( 123 );
                       tChildren[ 104 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 70 exists
                if( tNeighbors[ 70 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 70 ];

                    // test if neighbor 70 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   2 ] = tNeighbor->get_child(  20 );
                       tChildren[   3 ] = tNeighbor->get_child(  21 );
                       tChildren[   4 ] = tNeighbor->get_child(  22 );
                       tChildren[  27 ] = tNeighbor->get_child(  45 );
                       tChildren[  28 ] = tNeighbor->get_child(  46 );
                       tChildren[  29 ] = tNeighbor->get_child(  47 );
                       tChildren[  52 ] = tNeighbor->get_child(  70 );
                       tChildren[  53 ] = tNeighbor->get_child(  71 );
                       tChildren[  54 ] = tNeighbor->get_child(  72 );
                       tChildren[  77 ] = tNeighbor->get_child(  95 );
                       tChildren[  78 ] = tNeighbor->get_child(  96 );
                       tChildren[  79 ] = tNeighbor->get_child(  97 );
                       tChildren[ 102 ] = tNeighbor->get_child( 120 );
                       tChildren[ 103 ] = tNeighbor->get_child( 121 );
                       tChildren[ 104 ] = tNeighbor->get_child( 122 );
                    }
                }

                // test if neighbor 71 exists
                if( tNeighbors[ 71 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 71 ];

                    // test if neighbor 71 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  20 );
                       tChildren[  29 ] = tNeighbor->get_child(  45 );
                       tChildren[  54 ] = tNeighbor->get_child(  70 );
                       tChildren[  79 ] = tNeighbor->get_child(  95 );
                       tChildren[ 104 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 72 exists
                if( tNeighbors[ 72 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 72 ];

                    // test if neighbor 72 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(  14 );
                       tChildren[   5 ] = tNeighbor->get_child(  19 );
                       tChildren[  10 ] = tNeighbor->get_child(  24 );
                       tChildren[  25 ] = tNeighbor->get_child(  39 );
                       tChildren[  30 ] = tNeighbor->get_child(  44 );
                       tChildren[  35 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  64 );
                       tChildren[  55 ] = tNeighbor->get_child(  69 );
                       tChildren[  60 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  89 );
                       tChildren[  80 ] = tNeighbor->get_child(  94 );
                       tChildren[  85 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 114 );
                       tChildren[ 105 ] = tNeighbor->get_child( 119 );
                       tChildren[ 110 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 73 exists
                if( tNeighbors[ 73 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 73 ];

                    // test if neighbor 73 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(  10 );
                       tChildren[   9 ] = tNeighbor->get_child(  15 );
                       tChildren[  14 ] = tNeighbor->get_child(  20 );
                       tChildren[  29 ] = tNeighbor->get_child(  35 );
                       tChildren[  34 ] = tNeighbor->get_child(  40 );
                       tChildren[  39 ] = tNeighbor->get_child(  45 );
                       tChildren[  54 ] = tNeighbor->get_child(  60 );
                       tChildren[  59 ] = tNeighbor->get_child(  65 );
                       tChildren[  64 ] = tNeighbor->get_child(  70 );
                       tChildren[  79 ] = tNeighbor->get_child(  85 );
                       tChildren[  84 ] = tNeighbor->get_child(  90 );
                       tChildren[  89 ] = tNeighbor->get_child(  95 );
                       tChildren[ 104 ] = tNeighbor->get_child( 110 );
                       tChildren[ 109 ] = tNeighbor->get_child( 115 );
                       tChildren[ 114 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 74 exists
                if( tNeighbors[ 74 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 74 ];

                    // test if neighbor 74 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   0 ] = tNeighbor->get_child(   4 );
                       tChildren[   5 ] = tNeighbor->get_child(   9 );
                       tChildren[  10 ] = tNeighbor->get_child(  14 );
                       tChildren[  15 ] = tNeighbor->get_child(  19 );
                       tChildren[  20 ] = tNeighbor->get_child(  24 );
                       tChildren[  25 ] = tNeighbor->get_child(  29 );
                       tChildren[  30 ] = tNeighbor->get_child(  34 );
                       tChildren[  35 ] = tNeighbor->get_child(  39 );
                       tChildren[  40 ] = tNeighbor->get_child(  44 );
                       tChildren[  45 ] = tNeighbor->get_child(  49 );
                       tChildren[  50 ] = tNeighbor->get_child(  54 );
                       tChildren[  55 ] = tNeighbor->get_child(  59 );
                       tChildren[  60 ] = tNeighbor->get_child(  64 );
                       tChildren[  65 ] = tNeighbor->get_child(  69 );
                       tChildren[  70 ] = tNeighbor->get_child(  74 );
                       tChildren[  75 ] = tNeighbor->get_child(  79 );
                       tChildren[  80 ] = tNeighbor->get_child(  84 );
                       tChildren[  85 ] = tNeighbor->get_child(  89 );
                       tChildren[  90 ] = tNeighbor->get_child(  94 );
                       tChildren[  95 ] = tNeighbor->get_child(  99 );
                       tChildren[ 100 ] = tNeighbor->get_child( 104 );
                       tChildren[ 105 ] = tNeighbor->get_child( 109 );
                       tChildren[ 110 ] = tNeighbor->get_child( 114 );
                       tChildren[ 115 ] = tNeighbor->get_child( 119 );
                       tChildren[ 120 ] = tNeighbor->get_child( 124 );
                    }
                }

                // test if neighbor 75 exists
                if( tNeighbors[ 75 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 75 ];

                    // test if neighbor 75 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[   4 ] = tNeighbor->get_child(   0 );
                       tChildren[   9 ] = tNeighbor->get_child(   5 );
                       tChildren[  14 ] = tNeighbor->get_child(  10 );
                       tChildren[  19 ] = tNeighbor->get_child(  15 );
                       tChildren[  24 ] = tNeighbor->get_child(  20 );
                       tChildren[  29 ] = tNeighbor->get_child(  25 );
                       tChildren[  34 ] = tNeighbor->get_child(  30 );
                       tChildren[  39 ] = tNeighbor->get_child(  35 );
                       tChildren[  44 ] = tNeighbor->get_child(  40 );
                       tChildren[  49 ] = tNeighbor->get_child(  45 );
                       tChildren[  54 ] = tNeighbor->get_child(  50 );
                       tChildren[  59 ] = tNeighbor->get_child(  55 );
                       tChildren[  64 ] = tNeighbor->get_child(  60 );
                       tChildren[  69 ] = tNeighbor->get_child(  65 );
                       tChildren[  74 ] = tNeighbor->get_child(  70 );
                       tChildren[  79 ] = tNeighbor->get_child(  75 );
                       tChildren[  84 ] = tNeighbor->get_child(  80 );
                       tChildren[  89 ] = tNeighbor->get_child(  85 );
                       tChildren[  94 ] = tNeighbor->get_child(  90 );
                       tChildren[  99 ] = tNeighbor->get_child(  95 );
                       tChildren[ 104 ] = tNeighbor->get_child( 100 );
                       tChildren[ 109 ] = tNeighbor->get_child( 105 );
                       tChildren[ 114 ] = tNeighbor->get_child( 110 );
                       tChildren[ 119 ] = tNeighbor->get_child( 115 );
                       tChildren[ 124 ] = tNeighbor->get_child( 120 );
                    }
                }

                // test if neighbor 76 exists
                if( tNeighbors[ 76 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 76 ];

                    // test if neighbor 76 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  10 ] = tNeighbor->get_child(   4 );
                       tChildren[  15 ] = tNeighbor->get_child(   9 );
                       tChildren[  20 ] = tNeighbor->get_child(  14 );
                       tChildren[  35 ] = tNeighbor->get_child(  29 );
                       tChildren[  40 ] = tNeighbor->get_child(  34 );
                       tChildren[  45 ] = tNeighbor->get_child(  39 );
                       tChildren[  60 ] = tNeighbor->get_child(  54 );
                       tChildren[  65 ] = tNeighbor->get_child(  59 );
                       tChildren[  70 ] = tNeighbor->get_child(  64 );
                       tChildren[  85 ] = tNeighbor->get_child(  79 );
                       tChildren[  90 ] = tNeighbor->get_child(  84 );
                       tChildren[  95 ] = tNeighbor->get_child(  89 );
                       tChildren[ 110 ] = tNeighbor->get_child( 104 );
                       tChildren[ 115 ] = tNeighbor->get_child( 109 );
                       tChildren[ 120 ] = tNeighbor->get_child( 114 );
                    }
                }

                // test if neighbor 77 exists
                if( tNeighbors[ 77 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 77 ];

                    // test if neighbor 77 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  14 ] = tNeighbor->get_child(   0 );
                       tChildren[  19 ] = tNeighbor->get_child(   5 );
                       tChildren[  24 ] = tNeighbor->get_child(  10 );
                       tChildren[  39 ] = tNeighbor->get_child(  25 );
                       tChildren[  44 ] = tNeighbor->get_child(  30 );
                       tChildren[  49 ] = tNeighbor->get_child(  35 );
                       tChildren[  64 ] = tNeighbor->get_child(  50 );
                       tChildren[  69 ] = tNeighbor->get_child(  55 );
                       tChildren[  74 ] = tNeighbor->get_child(  60 );
                       tChildren[  89 ] = tNeighbor->get_child(  75 );
                       tChildren[  94 ] = tNeighbor->get_child(  80 );
                       tChildren[  99 ] = tNeighbor->get_child(  85 );
                       tChildren[ 114 ] = tNeighbor->get_child( 100 );
                       tChildren[ 119 ] = tNeighbor->get_child( 105 );
                       tChildren[ 124 ] = tNeighbor->get_child( 110 );
                    }
                }

                // test if neighbor 78 exists
                if( tNeighbors[ 78 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 78 ];

                    // test if neighbor 78 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   4 );
                       tChildren[  45 ] = tNeighbor->get_child(  29 );
                       tChildren[  70 ] = tNeighbor->get_child(  54 );
                       tChildren[  95 ] = tNeighbor->get_child(  79 );
                       tChildren[ 120 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 79 exists
                if( tNeighbors[ 79 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 79 ];

                    // test if neighbor 79 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   2 );
                       tChildren[  21 ] = tNeighbor->get_child(   3 );
                       tChildren[  22 ] = tNeighbor->get_child(   4 );
                       tChildren[  45 ] = tNeighbor->get_child(  27 );
                       tChildren[  46 ] = tNeighbor->get_child(  28 );
                       tChildren[  47 ] = tNeighbor->get_child(  29 );
                       tChildren[  70 ] = tNeighbor->get_child(  52 );
                       tChildren[  71 ] = tNeighbor->get_child(  53 );
                       tChildren[  72 ] = tNeighbor->get_child(  54 );
                       tChildren[  95 ] = tNeighbor->get_child(  77 );
                       tChildren[  96 ] = tNeighbor->get_child(  78 );
                       tChildren[  97 ] = tNeighbor->get_child(  79 );
                       tChildren[ 120 ] = tNeighbor->get_child( 102 );
                       tChildren[ 121 ] = tNeighbor->get_child( 103 );
                       tChildren[ 122 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 80 exists
                if( tNeighbors[ 80 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 80 ];

                    // test if neighbor 80 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  20 ] = tNeighbor->get_child(   0 );
                       tChildren[  21 ] = tNeighbor->get_child(   1 );
                       tChildren[  22 ] = tNeighbor->get_child(   2 );
                       tChildren[  23 ] = tNeighbor->get_child(   3 );
                       tChildren[  24 ] = tNeighbor->get_child(   4 );
                       tChildren[  45 ] = tNeighbor->get_child(  25 );
                       tChildren[  46 ] = tNeighbor->get_child(  26 );
                       tChildren[  47 ] = tNeighbor->get_child(  27 );
                       tChildren[  48 ] = tNeighbor->get_child(  28 );
                       tChildren[  49 ] = tNeighbor->get_child(  29 );
                       tChildren[  70 ] = tNeighbor->get_child(  50 );
                       tChildren[  71 ] = tNeighbor->get_child(  51 );
                       tChildren[  72 ] = tNeighbor->get_child(  52 );
                       tChildren[  73 ] = tNeighbor->get_child(  53 );
                       tChildren[  74 ] = tNeighbor->get_child(  54 );
                       tChildren[  95 ] = tNeighbor->get_child(  75 );
                       tChildren[  96 ] = tNeighbor->get_child(  76 );
                       tChildren[  97 ] = tNeighbor->get_child(  77 );
                       tChildren[  98 ] = tNeighbor->get_child(  78 );
                       tChildren[  99 ] = tNeighbor->get_child(  79 );
                       tChildren[ 120 ] = tNeighbor->get_child( 100 );
                       tChildren[ 121 ] = tNeighbor->get_child( 101 );
                       tChildren[ 122 ] = tNeighbor->get_child( 102 );
                       tChildren[ 123 ] = tNeighbor->get_child( 103 );
                       tChildren[ 124 ] = tNeighbor->get_child( 104 );
                    }
                }

                // test if neighbor 81 exists
                if( tNeighbors[ 81 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 81 ];

                    // test if neighbor 81 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  22 ] = tNeighbor->get_child(   0 );
                       tChildren[  23 ] = tNeighbor->get_child(   1 );
                       tChildren[  24 ] = tNeighbor->get_child(   2 );
                       tChildren[  47 ] = tNeighbor->get_child(  25 );
                       tChildren[  48 ] = tNeighbor->get_child(  26 );
                       tChildren[  49 ] = tNeighbor->get_child(  27 );
                       tChildren[  72 ] = tNeighbor->get_child(  50 );
                       tChildren[  73 ] = tNeighbor->get_child(  51 );
                       tChildren[  74 ] = tNeighbor->get_child(  52 );
                       tChildren[  97 ] = tNeighbor->get_child(  75 );
                       tChildren[  98 ] = tNeighbor->get_child(  76 );
                       tChildren[  99 ] = tNeighbor->get_child(  77 );
                       tChildren[ 122 ] = tNeighbor->get_child( 100 );
                       tChildren[ 123 ] = tNeighbor->get_child( 101 );
                       tChildren[ 124 ] = tNeighbor->get_child( 102 );
                    }
                }

                // test if neighbor 82 exists
                if( tNeighbors[ 82 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 82 ];

                    // test if neighbor 82 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  24 ] = tNeighbor->get_child(   0 );
                       tChildren[  49 ] = tNeighbor->get_child(  25 );
                       tChildren[  74 ] = tNeighbor->get_child(  50 );
                       tChildren[  99 ] = tNeighbor->get_child(  75 );
                       tChildren[ 124 ] = tNeighbor->get_child( 100 );
                    }
                }

                // test if neighbor 83 exists
                if( tNeighbors[ 83 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 83 ];

                    // test if neighbor 83 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 84 exists
                if( tNeighbors[ 84 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 84 ];

                    // test if neighbor 84 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  22 );
                       tChildren[  51 ] = tNeighbor->get_child(  23 );
                       tChildren[  52 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  47 );
                       tChildren[  76 ] = tNeighbor->get_child(  48 );
                       tChildren[  77 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  72 );
                       tChildren[ 101 ] = tNeighbor->get_child(  73 );
                       tChildren[ 102 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 85 exists
                if( tNeighbors[ 85 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 85 ];

                    // test if neighbor 85 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  20 );
                       tChildren[  51 ] = tNeighbor->get_child(  21 );
                       tChildren[  52 ] = tNeighbor->get_child(  22 );
                       tChildren[  53 ] = tNeighbor->get_child(  23 );
                       tChildren[  54 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  45 );
                       tChildren[  76 ] = tNeighbor->get_child(  46 );
                       tChildren[  77 ] = tNeighbor->get_child(  47 );
                       tChildren[  78 ] = tNeighbor->get_child(  48 );
                       tChildren[  79 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  70 );
                       tChildren[ 101 ] = tNeighbor->get_child(  71 );
                       tChildren[ 102 ] = tNeighbor->get_child(  72 );
                       tChildren[ 103 ] = tNeighbor->get_child(  73 );
                       tChildren[ 104 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 86 exists
                if( tNeighbors[ 86 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 86 ];

                    // test if neighbor 86 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  52 ] = tNeighbor->get_child(  20 );
                       tChildren[  53 ] = tNeighbor->get_child(  21 );
                       tChildren[  54 ] = tNeighbor->get_child(  22 );
                       tChildren[  77 ] = tNeighbor->get_child(  45 );
                       tChildren[  78 ] = tNeighbor->get_child(  46 );
                       tChildren[  79 ] = tNeighbor->get_child(  47 );
                       tChildren[ 102 ] = tNeighbor->get_child(  70 );
                       tChildren[ 103 ] = tNeighbor->get_child(  71 );
                       tChildren[ 104 ] = tNeighbor->get_child(  72 );
                    }
                }

                // test if neighbor 87 exists
                if( tNeighbors[ 87 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 87 ];

                    // test if neighbor 87 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  54 ] = tNeighbor->get_child(  20 );
                       tChildren[  79 ] = tNeighbor->get_child(  45 );
                       tChildren[ 104 ] = tNeighbor->get_child(  70 );
                    }
                }

                // test if neighbor 88 exists
                if( tNeighbors[ 88 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 88 ];

                    // test if neighbor 88 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(  14 );
                       tChildren[  55 ] = tNeighbor->get_child(  19 );
                       tChildren[  60 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  39 );
                       tChildren[  80 ] = tNeighbor->get_child(  44 );
                       tChildren[  85 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  64 );
                       tChildren[ 105 ] = tNeighbor->get_child(  69 );
                       tChildren[ 110 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 89 exists
                if( tNeighbors[ 89 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 89 ];

                    // test if neighbor 89 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  54 ] = tNeighbor->get_child(  10 );
                       tChildren[  59 ] = tNeighbor->get_child(  15 );
                       tChildren[  64 ] = tNeighbor->get_child(  20 );
                       tChildren[  79 ] = tNeighbor->get_child(  35 );
                       tChildren[  84 ] = tNeighbor->get_child(  40 );
                       tChildren[  89 ] = tNeighbor->get_child(  45 );
                       tChildren[ 104 ] = tNeighbor->get_child(  60 );
                       tChildren[ 109 ] = tNeighbor->get_child(  65 );
                       tChildren[ 114 ] = tNeighbor->get_child(  70 );
                    }
                }

                // test if neighbor 90 exists
                if( tNeighbors[ 90 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 90 ];

                    // test if neighbor 90 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  50 ] = tNeighbor->get_child(   4 );
                       tChildren[  55 ] = tNeighbor->get_child(   9 );
                       tChildren[  60 ] = tNeighbor->get_child(  14 );
                       tChildren[  65 ] = tNeighbor->get_child(  19 );
                       tChildren[  70 ] = tNeighbor->get_child(  24 );
                       tChildren[  75 ] = tNeighbor->get_child(  29 );
                       tChildren[  80 ] = tNeighbor->get_child(  34 );
                       tChildren[  85 ] = tNeighbor->get_child(  39 );
                       tChildren[  90 ] = tNeighbor->get_child(  44 );
                       tChildren[  95 ] = tNeighbor->get_child(  49 );
                       tChildren[ 100 ] = tNeighbor->get_child(  54 );
                       tChildren[ 105 ] = tNeighbor->get_child(  59 );
                       tChildren[ 110 ] = tNeighbor->get_child(  64 );
                       tChildren[ 115 ] = tNeighbor->get_child(  69 );
                       tChildren[ 120 ] = tNeighbor->get_child(  74 );
                    }
                }

                // test if neighbor 91 exists
                if( tNeighbors[ 91 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 91 ];

                    // test if neighbor 91 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  54 ] = tNeighbor->get_child(   0 );
                       tChildren[  59 ] = tNeighbor->get_child(   5 );
                       tChildren[  64 ] = tNeighbor->get_child(  10 );
                       tChildren[  69 ] = tNeighbor->get_child(  15 );
                       tChildren[  74 ] = tNeighbor->get_child(  20 );
                       tChildren[  79 ] = tNeighbor->get_child(  25 );
                       tChildren[  84 ] = tNeighbor->get_child(  30 );
                       tChildren[  89 ] = tNeighbor->get_child(  35 );
                       tChildren[  94 ] = tNeighbor->get_child(  40 );
                       tChildren[  99 ] = tNeighbor->get_child(  45 );
                       tChildren[ 104 ] = tNeighbor->get_child(  50 );
                       tChildren[ 109 ] = tNeighbor->get_child(  55 );
                       tChildren[ 114 ] = tNeighbor->get_child(  60 );
                       tChildren[ 119 ] = tNeighbor->get_child(  65 );
                       tChildren[ 124 ] = tNeighbor->get_child(  70 );
                    }
                }

                // test if neighbor 92 exists
                if( tNeighbors[ 92 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 92 ];

                    // test if neighbor 92 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  60 ] = tNeighbor->get_child(   4 );
                       tChildren[  65 ] = tNeighbor->get_child(   9 );
                       tChildren[  70 ] = tNeighbor->get_child(  14 );
                       tChildren[  85 ] = tNeighbor->get_child(  29 );
                       tChildren[  90 ] = tNeighbor->get_child(  34 );
                       tChildren[  95 ] = tNeighbor->get_child(  39 );
                       tChildren[ 110 ] = tNeighbor->get_child(  54 );
                       tChildren[ 115 ] = tNeighbor->get_child(  59 );
                       tChildren[ 120 ] = tNeighbor->get_child(  64 );
                    }
                }

                // test if neighbor 93 exists
                if( tNeighbors[ 93 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 93 ];

                    // test if neighbor 93 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  64 ] = tNeighbor->get_child(   0 );
                       tChildren[  69 ] = tNeighbor->get_child(   5 );
                       tChildren[  74 ] = tNeighbor->get_child(  10 );
                       tChildren[  89 ] = tNeighbor->get_child(  25 );
                       tChildren[  94 ] = tNeighbor->get_child(  30 );
                       tChildren[  99 ] = tNeighbor->get_child(  35 );
                       tChildren[ 114 ] = tNeighbor->get_child(  50 );
                       tChildren[ 119 ] = tNeighbor->get_child(  55 );
                       tChildren[ 124 ] = tNeighbor->get_child(  60 );
                    }
                }

                // test if neighbor 94 exists
                if( tNeighbors[ 94 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 94 ];

                    // test if neighbor 94 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  70 ] = tNeighbor->get_child(   4 );
                       tChildren[  95 ] = tNeighbor->get_child(  29 );
                       tChildren[ 120 ] = tNeighbor->get_child(  54 );
                    }
                }

                // test if neighbor 95 exists
                if( tNeighbors[ 95 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 95 ];

                    // test if neighbor 95 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  70 ] = tNeighbor->get_child(   2 );
                       tChildren[  71 ] = tNeighbor->get_child(   3 );
                       tChildren[  72 ] = tNeighbor->get_child(   4 );
                       tChildren[  95 ] = tNeighbor->get_child(  27 );
                       tChildren[  96 ] = tNeighbor->get_child(  28 );
                       tChildren[  97 ] = tNeighbor->get_child(  29 );
                       tChildren[ 120 ] = tNeighbor->get_child(  52 );
                       tChildren[ 121 ] = tNeighbor->get_child(  53 );
                       tChildren[ 122 ] = tNeighbor->get_child(  54 );
                    }
                }

                // test if neighbor 96 exists
                if( tNeighbors[ 96 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 96 ];

                    // test if neighbor 96 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  70 ] = tNeighbor->get_child(   0 );
                       tChildren[  71 ] = tNeighbor->get_child(   1 );
                       tChildren[  72 ] = tNeighbor->get_child(   2 );
                       tChildren[  73 ] = tNeighbor->get_child(   3 );
                       tChildren[  74 ] = tNeighbor->get_child(   4 );
                       tChildren[  95 ] = tNeighbor->get_child(  25 );
                       tChildren[  96 ] = tNeighbor->get_child(  26 );
                       tChildren[  97 ] = tNeighbor->get_child(  27 );
                       tChildren[  98 ] = tNeighbor->get_child(  28 );
                       tChildren[  99 ] = tNeighbor->get_child(  29 );
                       tChildren[ 120 ] = tNeighbor->get_child(  50 );
                       tChildren[ 121 ] = tNeighbor->get_child(  51 );
                       tChildren[ 122 ] = tNeighbor->get_child(  52 );
                       tChildren[ 123 ] = tNeighbor->get_child(  53 );
                       tChildren[ 124 ] = tNeighbor->get_child(  54 );
                    }
                }

                // test if neighbor 97 exists
                if( tNeighbors[ 97 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 97 ];

                    // test if neighbor 97 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  72 ] = tNeighbor->get_child(   0 );
                       tChildren[  73 ] = tNeighbor->get_child(   1 );
                       tChildren[  74 ] = tNeighbor->get_child(   2 );
                       tChildren[  97 ] = tNeighbor->get_child(  25 );
                       tChildren[  98 ] = tNeighbor->get_child(  26 );
                       tChildren[  99 ] = tNeighbor->get_child(  27 );
                       tChildren[ 122 ] = tNeighbor->get_child(  50 );
                       tChildren[ 123 ] = tNeighbor->get_child(  51 );
                       tChildren[ 124 ] = tNeighbor->get_child(  52 );
                    }
                }

                // test if neighbor 98 exists
                if( tNeighbors[ 98 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 98 ];

                    // test if neighbor 98 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[  74 ] = tNeighbor->get_child(   0 );
                       tChildren[  99 ] = tNeighbor->get_child(  25 );
                       tChildren[ 124 ] = tNeighbor->get_child(  50 );
                    }
                }

                // test if neighbor 99 exists
                if( tNeighbors[ 99 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 99 ];

                    // test if neighbor 99 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 100 exists
                if( tNeighbors[ 100 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 100 ];

                    // test if neighbor 100 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  22 );
                       tChildren[ 101 ] = tNeighbor->get_child(  23 );
                       tChildren[ 102 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 101 exists
                if( tNeighbors[ 101 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 101 ];

                    // test if neighbor 101 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  20 );
                       tChildren[ 101 ] = tNeighbor->get_child(  21 );
                       tChildren[ 102 ] = tNeighbor->get_child(  22 );
                       tChildren[ 103 ] = tNeighbor->get_child(  23 );
                       tChildren[ 104 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 102 exists
                if( tNeighbors[ 102 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 102 ];

                    // test if neighbor 102 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 102 ] = tNeighbor->get_child(  20 );
                       tChildren[ 103 ] = tNeighbor->get_child(  21 );
                       tChildren[ 104 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 103 exists
                if( tNeighbors[ 103 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 103 ];

                    // test if neighbor 103 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 104 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 104 exists
                if( tNeighbors[ 104 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 104 ];

                    // test if neighbor 104 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  14 );
                       tChildren[ 105 ] = tNeighbor->get_child(  19 );
                       tChildren[ 110 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 105 exists
                if( tNeighbors[ 105 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 105 ];

                    // test if neighbor 105 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  12 );
                       tChildren[ 101 ] = tNeighbor->get_child(  13 );
                       tChildren[ 102 ] = tNeighbor->get_child(  14 );
                       tChildren[ 105 ] = tNeighbor->get_child(  17 );
                       tChildren[ 106 ] = tNeighbor->get_child(  18 );
                       tChildren[ 107 ] = tNeighbor->get_child(  19 );
                       tChildren[ 110 ] = tNeighbor->get_child(  22 );
                       tChildren[ 111 ] = tNeighbor->get_child(  23 );
                       tChildren[ 112 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 106 exists
                if( tNeighbors[ 106 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 106 ];

                    // test if neighbor 106 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(  10 );
                       tChildren[ 101 ] = tNeighbor->get_child(  11 );
                       tChildren[ 102 ] = tNeighbor->get_child(  12 );
                       tChildren[ 103 ] = tNeighbor->get_child(  13 );
                       tChildren[ 104 ] = tNeighbor->get_child(  14 );
                       tChildren[ 105 ] = tNeighbor->get_child(  15 );
                       tChildren[ 106 ] = tNeighbor->get_child(  16 );
                       tChildren[ 107 ] = tNeighbor->get_child(  17 );
                       tChildren[ 108 ] = tNeighbor->get_child(  18 );
                       tChildren[ 109 ] = tNeighbor->get_child(  19 );
                       tChildren[ 110 ] = tNeighbor->get_child(  20 );
                       tChildren[ 111 ] = tNeighbor->get_child(  21 );
                       tChildren[ 112 ] = tNeighbor->get_child(  22 );
                       tChildren[ 113 ] = tNeighbor->get_child(  23 );
                       tChildren[ 114 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 107 exists
                if( tNeighbors[ 107 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 107 ];

                    // test if neighbor 107 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 102 ] = tNeighbor->get_child(  10 );
                       tChildren[ 103 ] = tNeighbor->get_child(  11 );
                       tChildren[ 104 ] = tNeighbor->get_child(  12 );
                       tChildren[ 107 ] = tNeighbor->get_child(  15 );
                       tChildren[ 108 ] = tNeighbor->get_child(  16 );
                       tChildren[ 109 ] = tNeighbor->get_child(  17 );
                       tChildren[ 112 ] = tNeighbor->get_child(  20 );
                       tChildren[ 113 ] = tNeighbor->get_child(  21 );
                       tChildren[ 114 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 108 exists
                if( tNeighbors[ 108 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 108 ];

                    // test if neighbor 108 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 104 ] = tNeighbor->get_child(  10 );
                       tChildren[ 109 ] = tNeighbor->get_child(  15 );
                       tChildren[ 114 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 109 exists
                if( tNeighbors[ 109 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 109 ];

                    // test if neighbor 109 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(   4 );
                       tChildren[ 105 ] = tNeighbor->get_child(   9 );
                       tChildren[ 110 ] = tNeighbor->get_child(  14 );
                       tChildren[ 115 ] = tNeighbor->get_child(  19 );
                       tChildren[ 120 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 110 exists
                if( tNeighbors[ 110 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 110 ];

                    // test if neighbor 110 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(   2 );
                       tChildren[ 101 ] = tNeighbor->get_child(   3 );
                       tChildren[ 102 ] = tNeighbor->get_child(   4 );
                       tChildren[ 105 ] = tNeighbor->get_child(   7 );
                       tChildren[ 106 ] = tNeighbor->get_child(   8 );
                       tChildren[ 107 ] = tNeighbor->get_child(   9 );
                       tChildren[ 110 ] = tNeighbor->get_child(  12 );
                       tChildren[ 111 ] = tNeighbor->get_child(  13 );
                       tChildren[ 112 ] = tNeighbor->get_child(  14 );
                       tChildren[ 115 ] = tNeighbor->get_child(  17 );
                       tChildren[ 116 ] = tNeighbor->get_child(  18 );
                       tChildren[ 117 ] = tNeighbor->get_child(  19 );
                       tChildren[ 120 ] = tNeighbor->get_child(  22 );
                       tChildren[ 121 ] = tNeighbor->get_child(  23 );
                       tChildren[ 122 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 111 exists
                if( tNeighbors[ 111 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 111 ];

                    // test if neighbor 111 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 100 ] = tNeighbor->get_child(   0 );
                       tChildren[ 101 ] = tNeighbor->get_child(   1 );
                       tChildren[ 102 ] = tNeighbor->get_child(   2 );
                       tChildren[ 103 ] = tNeighbor->get_child(   3 );
                       tChildren[ 104 ] = tNeighbor->get_child(   4 );
                       tChildren[ 105 ] = tNeighbor->get_child(   5 );
                       tChildren[ 106 ] = tNeighbor->get_child(   6 );
                       tChildren[ 107 ] = tNeighbor->get_child(   7 );
                       tChildren[ 108 ] = tNeighbor->get_child(   8 );
                       tChildren[ 109 ] = tNeighbor->get_child(   9 );
                       tChildren[ 110 ] = tNeighbor->get_child(  10 );
                       tChildren[ 111 ] = tNeighbor->get_child(  11 );
                       tChildren[ 112 ] = tNeighbor->get_child(  12 );
                       tChildren[ 113 ] = tNeighbor->get_child(  13 );
                       tChildren[ 114 ] = tNeighbor->get_child(  14 );
                       tChildren[ 115 ] = tNeighbor->get_child(  15 );
                       tChildren[ 116 ] = tNeighbor->get_child(  16 );
                       tChildren[ 117 ] = tNeighbor->get_child(  17 );
                       tChildren[ 118 ] = tNeighbor->get_child(  18 );
                       tChildren[ 119 ] = tNeighbor->get_child(  19 );
                       tChildren[ 120 ] = tNeighbor->get_child(  20 );
                       tChildren[ 121 ] = tNeighbor->get_child(  21 );
                       tChildren[ 122 ] = tNeighbor->get_child(  22 );
                       tChildren[ 123 ] = tNeighbor->get_child(  23 );
                       tChildren[ 124 ] = tNeighbor->get_child(  24 );
                    }
                }

                // test if neighbor 112 exists
                if( tNeighbors[ 112 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 112 ];

                    // test if neighbor 112 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 102 ] = tNeighbor->get_child(   0 );
                       tChildren[ 103 ] = tNeighbor->get_child(   1 );
                       tChildren[ 104 ] = tNeighbor->get_child(   2 );
                       tChildren[ 107 ] = tNeighbor->get_child(   5 );
                       tChildren[ 108 ] = tNeighbor->get_child(   6 );
                       tChildren[ 109 ] = tNeighbor->get_child(   7 );
                       tChildren[ 112 ] = tNeighbor->get_child(  10 );
                       tChildren[ 113 ] = tNeighbor->get_child(  11 );
                       tChildren[ 114 ] = tNeighbor->get_child(  12 );
                       tChildren[ 117 ] = tNeighbor->get_child(  15 );
                       tChildren[ 118 ] = tNeighbor->get_child(  16 );
                       tChildren[ 119 ] = tNeighbor->get_child(  17 );
                       tChildren[ 122 ] = tNeighbor->get_child(  20 );
                       tChildren[ 123 ] = tNeighbor->get_child(  21 );
                       tChildren[ 124 ] = tNeighbor->get_child(  22 );
                    }
                }

                // test if neighbor 113 exists
                if( tNeighbors[ 113 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 113 ];

                    // test if neighbor 113 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 104 ] = tNeighbor->get_child(   0 );
                       tChildren[ 109 ] = tNeighbor->get_child(   5 );
                       tChildren[ 114 ] = tNeighbor->get_child(  10 );
                       tChildren[ 119 ] = tNeighbor->get_child(  15 );
                       tChildren[ 124 ] = tNeighbor->get_child(  20 );
                    }
                }

                // test if neighbor 114 exists
                if( tNeighbors[ 114 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 114 ];

                    // test if neighbor 114 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 110 ] = tNeighbor->get_child(   4 );
                       tChildren[ 115 ] = tNeighbor->get_child(   9 );
                       tChildren[ 120 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 115 exists
                if( tNeighbors[ 115 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 115 ];

                    // test if neighbor 115 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 110 ] = tNeighbor->get_child(   2 );
                       tChildren[ 111 ] = tNeighbor->get_child(   3 );
                       tChildren[ 112 ] = tNeighbor->get_child(   4 );
                       tChildren[ 115 ] = tNeighbor->get_child(   7 );
                       tChildren[ 116 ] = tNeighbor->get_child(   8 );
                       tChildren[ 117 ] = tNeighbor->get_child(   9 );
                       tChildren[ 120 ] = tNeighbor->get_child(  12 );
                       tChildren[ 121 ] = tNeighbor->get_child(  13 );
                       tChildren[ 122 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 116 exists
                if( tNeighbors[ 116 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 116 ];

                    // test if neighbor 116 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 110 ] = tNeighbor->get_child(   0 );
                       tChildren[ 111 ] = tNeighbor->get_child(   1 );
                       tChildren[ 112 ] = tNeighbor->get_child(   2 );
                       tChildren[ 113 ] = tNeighbor->get_child(   3 );
                       tChildren[ 114 ] = tNeighbor->get_child(   4 );
                       tChildren[ 115 ] = tNeighbor->get_child(   5 );
                       tChildren[ 116 ] = tNeighbor->get_child(   6 );
                       tChildren[ 117 ] = tNeighbor->get_child(   7 );
                       tChildren[ 118 ] = tNeighbor->get_child(   8 );
                       tChildren[ 119 ] = tNeighbor->get_child(   9 );
                       tChildren[ 120 ] = tNeighbor->get_child(  10 );
                       tChildren[ 121 ] = tNeighbor->get_child(  11 );
                       tChildren[ 122 ] = tNeighbor->get_child(  12 );
                       tChildren[ 123 ] = tNeighbor->get_child(  13 );
                       tChildren[ 124 ] = tNeighbor->get_child(  14 );
                    }
                }

                // test if neighbor 117 exists
                if( tNeighbors[ 117 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 117 ];

                    // test if neighbor 117 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 112 ] = tNeighbor->get_child(   0 );
                       tChildren[ 113 ] = tNeighbor->get_child(   1 );
                       tChildren[ 114 ] = tNeighbor->get_child(   2 );
                       tChildren[ 117 ] = tNeighbor->get_child(   5 );
                       tChildren[ 118 ] = tNeighbor->get_child(   6 );
                       tChildren[ 119 ] = tNeighbor->get_child(   7 );
                       tChildren[ 122 ] = tNeighbor->get_child(  10 );
                       tChildren[ 123 ] = tNeighbor->get_child(  11 );
                       tChildren[ 124 ] = tNeighbor->get_child(  12 );
                    }
                }

                // test if neighbor 118 exists
                if( tNeighbors[ 118 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 118 ];

                    // test if neighbor 118 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 114 ] = tNeighbor->get_child(   0 );
                       tChildren[ 119 ] = tNeighbor->get_child(   5 );
                       tChildren[ 124 ] = tNeighbor->get_child(  10 );
                    }
                }

                // test if neighbor 119 exists
                if( tNeighbors[ 119 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 119 ];

                    // test if neighbor 119 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 120 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 120 exists
                if( tNeighbors[ 120 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 120 ];

                    // test if neighbor 120 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 120 ] = tNeighbor->get_child(   2 );
                       tChildren[ 121 ] = tNeighbor->get_child(   3 );
                       tChildren[ 122 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 121 exists
                if( tNeighbors[ 121 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 121 ];

                    // test if neighbor 121 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 120 ] = tNeighbor->get_child(   0 );
                       tChildren[ 121 ] = tNeighbor->get_child(   1 );
                       tChildren[ 122 ] = tNeighbor->get_child(   2 );
                       tChildren[ 123 ] = tNeighbor->get_child(   3 );
                       tChildren[ 124 ] = tNeighbor->get_child(   4 );
                    }
                }

                // test if neighbor 122 exists
                if( tNeighbors[ 122 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 122 ];

                    // test if neighbor 122 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 122 ] = tNeighbor->get_child(   0 );
                       tChildren[ 123 ] = tNeighbor->get_child(   1 );
                       tChildren[ 124 ] = tNeighbor->get_child(   2 );
                    }
                }

                // test if neighbor 123 exists
                if( tNeighbors[ 123 ] != nullptr )
                {
                    // get pointer to neighbor
                    Basis_Function* tNeighbor = tNeighbors[ 123 ];

                    // test if neighbor 123 has children and copy them if so
                    if( tNeighbor->has_children() )
                    {
                       tChildren[ 124 ] = tNeighbor->get_child(   0 );
                    }
                }

                // level of child basis
                uint tLevel = tBasis->get_level() + 1;

                // create container for children
                tBasis->init_children_container();

                // position of basis
                const luint* tParentIJK  = tBasis->get_ijk();

                // minumum i-position
                luint tIMin = 2*tParentIJK[ 0 ];

                // minumum j-position
                luint tJMin = 2*tParentIJK[ 1 ];

                // minumum k-position
                luint tKMin = 2*tParentIJK[ 2 ];

                // maximum i-position
                luint tIMax = tIMin + 5;

                // maximum j-position
                luint tJMax = tJMin + 5;

                // maximum K-position
                luint tKMax = tKMin + 5;

                // initialize counter
                uint tChildIndex = 0;

                // loop over all positions
                for( luint k=tKMin; k<tKMax; ++k )
                {
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
                                 // calculate i-j-k position of child
                                 luint tIJK[ 3 ] = { i, j, k };

                                 // create child
                                 tBasis->insert_child( tChildIndex,
                                     new BSpline< 3, 3, 3 >( tIJK, tLevel, gNoProcOwner ) );

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
    luint BSpline_Element< 3, 3, 3 >::refine( moris::Cell< Element* > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // refine basis if they have not been refined already
        for( uint k=0; k<64; ++k )
        {
            tBasisCounter += this->refine_basis( k );
        }

        // initialize temporary basis pattern
        Basis_Function* tBasis[ 125 ] = { nullptr };

        // populate basis pattern
        if ( mBasis[   0 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   0 ]->get_child(  93 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   0 ]->get_child(  94 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   0 ]->get_child(  98 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   0 ]->get_child(  99 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[   0 ]->get_child( 118 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   0 ]->get_child( 119 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[   0 ]->get_child( 123 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[   0 ]->get_child( 124 );
            }
        }

        if ( mBasis[   1 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   1 ]->get_child(  90 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   1 ]->get_child(  91 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   1 ]->get_child(  95 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   1 ]->get_child(  96 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[   1 ]->get_child( 115 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[   1 ]->get_child( 116 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[   1 ]->get_child( 120 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[   1 ]->get_child( 121 );
            }
        }

        if ( mBasis[   2 ] != nullptr )
        {
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[   2 ]->get_child(  75 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[   2 ]->get_child(  76 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[   2 ]->get_child(  80 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[   2 ]->get_child(  81 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[   2 ]->get_child( 100 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[   2 ]->get_child( 101 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[   2 ]->get_child( 105 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[   2 ]->get_child( 106 );
            }
        }

        if ( mBasis[   3 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   3 ]->get_child(  78 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[   3 ]->get_child(  79 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[   3 ]->get_child(  83 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[   3 ]->get_child(  84 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[   3 ]->get_child( 103 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[   3 ]->get_child( 104 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[   3 ]->get_child( 108 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[   3 ]->get_child( 109 );
            }
        }

        if ( mBasis[   4 ] != nullptr )
        {
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[   4 ]->get_child(  18 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[   4 ]->get_child(  19 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[   4 ]->get_child(  23 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[   4 ]->get_child(  24 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[   4 ]->get_child(  43 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[   4 ]->get_child(  44 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[   4 ]->get_child(  48 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[   4 ]->get_child(  49 );
            }
        }

        if ( mBasis[   5 ] != nullptr )
        {
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[   5 ]->get_child(  15 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[   5 ]->get_child(  16 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[   5 ]->get_child(  20 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[   5 ]->get_child(  21 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[   5 ]->get_child(  40 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[   5 ]->get_child(  41 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[   5 ]->get_child(  45 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[   5 ]->get_child(  46 );
            }
        }

        if ( mBasis[   6 ] != nullptr )
        {
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[   6 ]->get_child(   0 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[   6 ]->get_child(   1 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[   6 ]->get_child(   5 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[   6 ]->get_child(   6 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[   6 ]->get_child(  25 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[   6 ]->get_child(  26 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[   6 ]->get_child(  30 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[   6 ]->get_child(  31 );
            }
        }

        if ( mBasis[   7 ] != nullptr )
        {
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[   7 ]->get_child(   3 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[   7 ]->get_child(   4 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[   7 ]->get_child(   8 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[   7 ]->get_child(   9 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[   7 ]->get_child(  28 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[   7 ]->get_child(  29 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[   7 ]->get_child(  33 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[   7 ]->get_child(  34 );
            }
        }

        if ( mBasis[   8 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   8 ]->get_child(  91 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   8 ]->get_child(  92 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   8 ]->get_child(  93 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   8 ]->get_child(  94 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   8 ]->get_child(  96 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   8 ]->get_child(  97 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   8 ]->get_child(  98 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   8 ]->get_child(  99 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[   8 ]->get_child( 116 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   8 ]->get_child( 117 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[   8 ]->get_child( 118 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[   8 ]->get_child( 119 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[   8 ]->get_child( 121 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[   8 ]->get_child( 122 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[   8 ]->get_child( 123 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[   8 ]->get_child( 124 );
            }
        }

        if ( mBasis[   9 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   9 ]->get_child(  90 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   9 ]->get_child(  91 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   9 ]->get_child(  92 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   9 ]->get_child(  93 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   9 ]->get_child(  95 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   9 ]->get_child(  96 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   9 ]->get_child(  97 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   9 ]->get_child(  98 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[   9 ]->get_child( 115 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[   9 ]->get_child( 116 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[   9 ]->get_child( 117 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[   9 ]->get_child( 118 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[   9 ]->get_child( 120 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[   9 ]->get_child( 121 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[   9 ]->get_child( 122 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[   9 ]->get_child( 123 );
            }
        }

        if ( mBasis[  10 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  10 ]->get_child(  83 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  10 ]->get_child(  84 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  10 ]->get_child(  88 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  10 ]->get_child(  89 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  10 ]->get_child(  93 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  10 ]->get_child(  94 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  10 ]->get_child(  98 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  10 ]->get_child(  99 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  10 ]->get_child( 108 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  10 ]->get_child( 109 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  10 ]->get_child( 113 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  10 ]->get_child( 114 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  10 ]->get_child( 118 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  10 ]->get_child( 119 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  10 ]->get_child( 123 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  10 ]->get_child( 124 );
            }
        }

        if ( mBasis[  11 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  11 ]->get_child(  78 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  11 ]->get_child(  79 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  11 ]->get_child(  83 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  11 ]->get_child(  84 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  11 ]->get_child(  88 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  11 ]->get_child(  89 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  11 ]->get_child(  93 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  11 ]->get_child(  94 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  11 ]->get_child( 103 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  11 ]->get_child( 104 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  11 ]->get_child( 108 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  11 ]->get_child( 109 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  11 ]->get_child( 113 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  11 ]->get_child( 114 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  11 ]->get_child( 118 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  11 ]->get_child( 119 );
            }
        }

        if ( mBasis[  12 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  12 ]->get_child(  43 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  12 ]->get_child(  44 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  12 ]->get_child(  48 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  12 ]->get_child(  49 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  12 ]->get_child(  68 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  12 ]->get_child(  69 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  12 ]->get_child(  73 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  12 ]->get_child(  74 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  12 ]->get_child(  93 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  12 ]->get_child(  94 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  12 ]->get_child(  98 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  12 ]->get_child(  99 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  12 ]->get_child( 118 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  12 ]->get_child( 119 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  12 ]->get_child( 123 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  12 ]->get_child( 124 );
            }
        }

        if ( mBasis[  13 ] != nullptr )
        {
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  13 ]->get_child(  18 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  13 ]->get_child(  19 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  13 ]->get_child(  23 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  13 ]->get_child(  24 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  13 ]->get_child(  43 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  13 ]->get_child(  44 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  13 ]->get_child(  48 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  13 ]->get_child(  49 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  13 ]->get_child(  68 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  13 ]->get_child(  69 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  13 ]->get_child(  73 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  13 ]->get_child(  74 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  13 ]->get_child(  93 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  13 ]->get_child(  94 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  13 ]->get_child(  98 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  13 ]->get_child(  99 );
            }
        }

        if ( mBasis[  14 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  14 ]->get_child(  80 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  14 ]->get_child(  81 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  14 ]->get_child(  85 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  14 ]->get_child(  86 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  14 ]->get_child(  90 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  14 ]->get_child(  91 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  14 ]->get_child(  95 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  14 ]->get_child(  96 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  14 ]->get_child( 105 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  14 ]->get_child( 106 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  14 ]->get_child( 110 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  14 ]->get_child( 111 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  14 ]->get_child( 115 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  14 ]->get_child( 116 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  14 ]->get_child( 120 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  14 ]->get_child( 121 );
            }
        }

        if ( mBasis[  15 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  15 ]->get_child(  75 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  15 ]->get_child(  76 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  15 ]->get_child(  80 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  15 ]->get_child(  81 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  15 ]->get_child(  85 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  15 ]->get_child(  86 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  15 ]->get_child(  90 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  15 ]->get_child(  91 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  15 ]->get_child( 100 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  15 ]->get_child( 101 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  15 ]->get_child( 105 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  15 ]->get_child( 106 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  15 ]->get_child( 110 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  15 ]->get_child( 111 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  15 ]->get_child( 115 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  15 ]->get_child( 116 );
            }
        }

        if ( mBasis[  16 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  16 ]->get_child(  40 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  16 ]->get_child(  41 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  16 ]->get_child(  45 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  16 ]->get_child(  46 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  16 ]->get_child(  65 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  16 ]->get_child(  66 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  16 ]->get_child(  70 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  16 ]->get_child(  71 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  16 ]->get_child(  90 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  16 ]->get_child(  91 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  16 ]->get_child(  95 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  16 ]->get_child(  96 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  16 ]->get_child( 115 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  16 ]->get_child( 116 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  16 ]->get_child( 120 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  16 ]->get_child( 121 );
            }
        }

        if ( mBasis[  17 ] != nullptr )
        {
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  17 ]->get_child(  15 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  17 ]->get_child(  16 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  17 ]->get_child(  20 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  17 ]->get_child(  21 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  17 ]->get_child(  40 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  17 ]->get_child(  41 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  17 ]->get_child(  45 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  17 ]->get_child(  46 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  17 ]->get_child(  65 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  17 ]->get_child(  66 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  17 ]->get_child(  70 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  17 ]->get_child(  71 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  17 ]->get_child(  90 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  17 ]->get_child(  91 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  17 ]->get_child(  95 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  17 ]->get_child(  96 );
            }
        }

        if ( mBasis[  18 ] != nullptr )
        {
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  18 ]->get_child(  75 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  18 ]->get_child(  76 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  18 ]->get_child(  77 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  18 ]->get_child(  78 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  18 ]->get_child(  80 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  18 ]->get_child(  81 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  18 ]->get_child(  82 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  18 ]->get_child(  83 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  18 ]->get_child( 100 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  18 ]->get_child( 101 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  18 ]->get_child( 102 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  18 ]->get_child( 103 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  18 ]->get_child( 105 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  18 ]->get_child( 106 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  18 ]->get_child( 107 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  18 ]->get_child( 108 );
            }
        }

        if ( mBasis[  19 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  19 ]->get_child(  76 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  19 ]->get_child(  77 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  19 ]->get_child(  78 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  19 ]->get_child(  79 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  19 ]->get_child(  81 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  19 ]->get_child(  82 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  19 ]->get_child(  83 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  19 ]->get_child(  84 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  19 ]->get_child( 101 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  19 ]->get_child( 102 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  19 ]->get_child( 103 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  19 ]->get_child( 104 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  19 ]->get_child( 106 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  19 ]->get_child( 107 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  19 ]->get_child( 108 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  19 ]->get_child( 109 );
            }
        }

        if ( mBasis[  20 ] != nullptr )
        {
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  20 ]->get_child(  25 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  20 ]->get_child(  26 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  20 ]->get_child(  30 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  20 ]->get_child(  31 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  20 ]->get_child(  50 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  20 ]->get_child(  51 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  20 ]->get_child(  55 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  20 ]->get_child(  56 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  20 ]->get_child(  75 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  20 ]->get_child(  76 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  20 ]->get_child(  80 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  20 ]->get_child(  81 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  20 ]->get_child( 100 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  20 ]->get_child( 101 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  20 ]->get_child( 105 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  20 ]->get_child( 106 );
            }
        }

        if ( mBasis[  21 ] != nullptr )
        {
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  21 ]->get_child(   0 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  21 ]->get_child(   1 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  21 ]->get_child(   5 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  21 ]->get_child(   6 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  21 ]->get_child(  25 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  21 ]->get_child(  26 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  21 ]->get_child(  30 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  21 ]->get_child(  31 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  21 ]->get_child(  50 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  21 ]->get_child(  51 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  21 ]->get_child(  55 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  21 ]->get_child(  56 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  21 ]->get_child(  75 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  21 ]->get_child(  76 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  21 ]->get_child(  80 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  21 ]->get_child(  81 );
            }
        }

        if ( mBasis[  22 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  22 ]->get_child(  28 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  22 ]->get_child(  29 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  22 ]->get_child(  33 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  22 ]->get_child(  34 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  22 ]->get_child(  53 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  22 ]->get_child(  54 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  22 ]->get_child(  58 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  22 ]->get_child(  59 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  22 ]->get_child(  78 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  22 ]->get_child(  79 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  22 ]->get_child(  83 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  22 ]->get_child(  84 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  22 ]->get_child( 103 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  22 ]->get_child( 104 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  22 ]->get_child( 108 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  22 ]->get_child( 109 );
            }
        }

        if ( mBasis[  23 ] != nullptr )
        {
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  23 ]->get_child(   3 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  23 ]->get_child(   4 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  23 ]->get_child(   8 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  23 ]->get_child(   9 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  23 ]->get_child(  28 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  23 ]->get_child(  29 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  23 ]->get_child(  33 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  23 ]->get_child(  34 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  23 ]->get_child(  53 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  23 ]->get_child(  54 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  23 ]->get_child(  58 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  23 ]->get_child(  59 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  23 ]->get_child(  78 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  23 ]->get_child(  79 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  23 ]->get_child(  83 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  23 ]->get_child(  84 );
            }
        }

        if ( mBasis[  24 ] != nullptr )
        {
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  24 ]->get_child(  16 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  24 ]->get_child(  17 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  24 ]->get_child(  18 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  24 ]->get_child(  19 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  24 ]->get_child(  21 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  24 ]->get_child(  22 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  24 ]->get_child(  23 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  24 ]->get_child(  24 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  24 ]->get_child(  41 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  24 ]->get_child(  42 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  24 ]->get_child(  43 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  24 ]->get_child(  44 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  24 ]->get_child(  46 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  24 ]->get_child(  47 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  24 ]->get_child(  48 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  24 ]->get_child(  49 );
            }
        }

        if ( mBasis[  25 ] != nullptr )
        {
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  25 ]->get_child(  15 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  25 ]->get_child(  16 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  25 ]->get_child(  17 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  25 ]->get_child(  18 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  25 ]->get_child(  20 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  25 ]->get_child(  21 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  25 ]->get_child(  22 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  25 ]->get_child(  23 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  25 ]->get_child(  40 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  25 ]->get_child(  41 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  25 ]->get_child(  42 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  25 ]->get_child(  43 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  25 ]->get_child(  45 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  25 ]->get_child(  46 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  25 ]->get_child(  47 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  25 ]->get_child(  48 );
            }
        }

        if ( mBasis[  26 ] != nullptr )
        {
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  26 ]->get_child(   8 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  26 ]->get_child(   9 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  26 ]->get_child(  13 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  26 ]->get_child(  14 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  26 ]->get_child(  18 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  26 ]->get_child(  19 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  26 ]->get_child(  23 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  26 ]->get_child(  24 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  26 ]->get_child(  33 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  26 ]->get_child(  34 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  26 ]->get_child(  38 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  26 ]->get_child(  39 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  26 ]->get_child(  43 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  26 ]->get_child(  44 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  26 ]->get_child(  48 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  26 ]->get_child(  49 );
            }
        }

        if ( mBasis[  27 ] != nullptr )
        {
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  27 ]->get_child(   3 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  27 ]->get_child(   4 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  27 ]->get_child(   8 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  27 ]->get_child(   9 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  27 ]->get_child(  13 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  27 ]->get_child(  14 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  27 ]->get_child(  18 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  27 ]->get_child(  19 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  27 ]->get_child(  28 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  27 ]->get_child(  29 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  27 ]->get_child(  33 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  27 ]->get_child(  34 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  27 ]->get_child(  38 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  27 ]->get_child(  39 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  27 ]->get_child(  43 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  27 ]->get_child(  44 );
            }
        }

        if ( mBasis[  28 ] != nullptr )
        {
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  28 ]->get_child(   5 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  28 ]->get_child(   6 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  28 ]->get_child(  10 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  28 ]->get_child(  11 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  28 ]->get_child(  15 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  28 ]->get_child(  16 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  28 ]->get_child(  20 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  28 ]->get_child(  21 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  28 ]->get_child(  30 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  28 ]->get_child(  31 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  28 ]->get_child(  35 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  28 ]->get_child(  36 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  28 ]->get_child(  40 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  28 ]->get_child(  41 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  28 ]->get_child(  45 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  28 ]->get_child(  46 );
            }
        }

        if ( mBasis[  29 ] != nullptr )
        {
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  29 ]->get_child(   0 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  29 ]->get_child(   1 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  29 ]->get_child(   5 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  29 ]->get_child(   6 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  29 ]->get_child(  10 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  29 ]->get_child(  11 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  29 ]->get_child(  15 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  29 ]->get_child(  16 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  29 ]->get_child(  25 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  29 ]->get_child(  26 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  29 ]->get_child(  30 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  29 ]->get_child(  31 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  29 ]->get_child(  35 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  29 ]->get_child(  36 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  29 ]->get_child(  40 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  29 ]->get_child(  41 );
            }
        }

        if ( mBasis[  30 ] != nullptr )
        {
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  30 ]->get_child(   0 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  30 ]->get_child(   1 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  30 ]->get_child(   2 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  30 ]->get_child(   3 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  30 ]->get_child(   5 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  30 ]->get_child(   6 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  30 ]->get_child(   7 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  30 ]->get_child(   8 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  30 ]->get_child(  25 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  30 ]->get_child(  26 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  30 ]->get_child(  27 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  30 ]->get_child(  28 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  30 ]->get_child(  30 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  30 ]->get_child(  31 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  30 ]->get_child(  32 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  30 ]->get_child(  33 );
            }
        }

        if ( mBasis[  31 ] != nullptr )
        {
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  31 ]->get_child(   1 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  31 ]->get_child(   2 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  31 ]->get_child(   3 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  31 ]->get_child(   4 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  31 ]->get_child(   6 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  31 ]->get_child(   7 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  31 ]->get_child(   8 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  31 ]->get_child(   9 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  31 ]->get_child(  26 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  31 ]->get_child(  27 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  31 ]->get_child(  28 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  31 ]->get_child(  29 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  31 ]->get_child(  31 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  31 ]->get_child(  32 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  31 ]->get_child(  33 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  31 ]->get_child(  34 );
            }
        }

        if ( mBasis[  32 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  32 ]->get_child(  81 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  32 ]->get_child(  82 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  32 ]->get_child(  83 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  32 ]->get_child(  84 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  32 ]->get_child(  86 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  32 ]->get_child(  87 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  32 ]->get_child(  88 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  32 ]->get_child(  89 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  32 ]->get_child(  91 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  32 ]->get_child(  92 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  32 ]->get_child(  93 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  32 ]->get_child(  94 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  32 ]->get_child(  96 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  32 ]->get_child(  97 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  32 ]->get_child(  98 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  32 ]->get_child(  99 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  32 ]->get_child( 106 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  32 ]->get_child( 107 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  32 ]->get_child( 108 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  32 ]->get_child( 109 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  32 ]->get_child( 111 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  32 ]->get_child( 112 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  32 ]->get_child( 113 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  32 ]->get_child( 114 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  32 ]->get_child( 116 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  32 ]->get_child( 117 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  32 ]->get_child( 118 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  32 ]->get_child( 119 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  32 ]->get_child( 121 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  32 ]->get_child( 122 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  32 ]->get_child( 123 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  32 ]->get_child( 124 );
            }
        }

        if ( mBasis[  33 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  33 ]->get_child(  76 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  33 ]->get_child(  77 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  33 ]->get_child(  78 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  33 ]->get_child(  79 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  33 ]->get_child(  81 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  33 ]->get_child(  82 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  33 ]->get_child(  83 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  33 ]->get_child(  84 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  33 ]->get_child(  86 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  33 ]->get_child(  87 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  33 ]->get_child(  88 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  33 ]->get_child(  89 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  33 ]->get_child(  91 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  33 ]->get_child(  92 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  33 ]->get_child(  93 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  33 ]->get_child(  94 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  33 ]->get_child( 101 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  33 ]->get_child( 102 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  33 ]->get_child( 103 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  33 ]->get_child( 104 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  33 ]->get_child( 106 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  33 ]->get_child( 107 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  33 ]->get_child( 108 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  33 ]->get_child( 109 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  33 ]->get_child( 111 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  33 ]->get_child( 112 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  33 ]->get_child( 113 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  33 ]->get_child( 114 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  33 ]->get_child( 116 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  33 ]->get_child( 117 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  33 ]->get_child( 118 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  33 ]->get_child( 119 );
            }
        }

        if ( mBasis[  34 ] != nullptr )
        {
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  34 ]->get_child(  75 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  34 ]->get_child(  76 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  34 ]->get_child(  77 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  34 ]->get_child(  78 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  34 ]->get_child(  80 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  34 ]->get_child(  81 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  34 ]->get_child(  82 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  34 ]->get_child(  83 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  34 ]->get_child(  85 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  34 ]->get_child(  86 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  34 ]->get_child(  87 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  34 ]->get_child(  88 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  34 ]->get_child(  90 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  34 ]->get_child(  91 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  34 ]->get_child(  92 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  34 ]->get_child(  93 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  34 ]->get_child( 100 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  34 ]->get_child( 101 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  34 ]->get_child( 102 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  34 ]->get_child( 103 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  34 ]->get_child( 105 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  34 ]->get_child( 106 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  34 ]->get_child( 107 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  34 ]->get_child( 108 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  34 ]->get_child( 110 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  34 ]->get_child( 111 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  34 ]->get_child( 112 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  34 ]->get_child( 113 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  34 ]->get_child( 115 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  34 ]->get_child( 116 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  34 ]->get_child( 117 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  34 ]->get_child( 118 );
            }
        }

        if ( mBasis[  35 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  35 ]->get_child(  80 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  35 ]->get_child(  81 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  35 ]->get_child(  82 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  35 ]->get_child(  83 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  35 ]->get_child(  85 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  35 ]->get_child(  86 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  35 ]->get_child(  87 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  35 ]->get_child(  88 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  35 ]->get_child(  90 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  35 ]->get_child(  91 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  35 ]->get_child(  92 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  35 ]->get_child(  93 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  35 ]->get_child(  95 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  35 ]->get_child(  96 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  35 ]->get_child(  97 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  35 ]->get_child(  98 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  35 ]->get_child( 105 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  35 ]->get_child( 106 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  35 ]->get_child( 107 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  35 ]->get_child( 108 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  35 ]->get_child( 110 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  35 ]->get_child( 111 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  35 ]->get_child( 112 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  35 ]->get_child( 113 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  35 ]->get_child( 115 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  35 ]->get_child( 116 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  35 ]->get_child( 117 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  35 ]->get_child( 118 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  35 ]->get_child( 120 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  35 ]->get_child( 121 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  35 ]->get_child( 122 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  35 ]->get_child( 123 );
            }
        }

        if ( mBasis[  36 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  36 ]->get_child(  41 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  36 ]->get_child(  42 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  36 ]->get_child(  43 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  36 ]->get_child(  44 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  36 ]->get_child(  46 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  36 ]->get_child(  47 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  36 ]->get_child(  48 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  36 ]->get_child(  49 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  36 ]->get_child(  66 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  36 ]->get_child(  67 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  36 ]->get_child(  68 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  36 ]->get_child(  69 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  36 ]->get_child(  71 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  36 ]->get_child(  72 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  36 ]->get_child(  73 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  36 ]->get_child(  74 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  36 ]->get_child(  91 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  36 ]->get_child(  92 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  36 ]->get_child(  93 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  36 ]->get_child(  94 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  36 ]->get_child(  96 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  36 ]->get_child(  97 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  36 ]->get_child(  98 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  36 ]->get_child(  99 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  36 ]->get_child( 116 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  36 ]->get_child( 117 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  36 ]->get_child( 118 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  36 ]->get_child( 119 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  36 ]->get_child( 121 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  36 ]->get_child( 122 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  36 ]->get_child( 123 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  36 ]->get_child( 124 );
            }
        }

        if ( mBasis[  37 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  37 ]->get_child(  40 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  37 ]->get_child(  41 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  37 ]->get_child(  42 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  37 ]->get_child(  43 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  37 ]->get_child(  45 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  37 ]->get_child(  46 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  37 ]->get_child(  47 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  37 ]->get_child(  48 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  37 ]->get_child(  65 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  37 ]->get_child(  66 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  37 ]->get_child(  67 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  37 ]->get_child(  68 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  37 ]->get_child(  70 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  37 ]->get_child(  71 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  37 ]->get_child(  72 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  37 ]->get_child(  73 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  37 ]->get_child(  90 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  37 ]->get_child(  91 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  37 ]->get_child(  92 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  37 ]->get_child(  93 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  37 ]->get_child(  95 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  37 ]->get_child(  96 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  37 ]->get_child(  97 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  37 ]->get_child(  98 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  37 ]->get_child( 115 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  37 ]->get_child( 116 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  37 ]->get_child( 117 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  37 ]->get_child( 118 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  37 ]->get_child( 120 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  37 ]->get_child( 121 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  37 ]->get_child( 122 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  37 ]->get_child( 123 );
            }
        }

        if ( mBasis[  38 ] != nullptr )
        {
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  38 ]->get_child(  15 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  38 ]->get_child(  16 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  38 ]->get_child(  17 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  38 ]->get_child(  18 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  38 ]->get_child(  20 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  38 ]->get_child(  21 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  38 ]->get_child(  22 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  38 ]->get_child(  23 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  38 ]->get_child(  40 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  38 ]->get_child(  41 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  38 ]->get_child(  42 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  38 ]->get_child(  43 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  38 ]->get_child(  45 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  38 ]->get_child(  46 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  38 ]->get_child(  47 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  38 ]->get_child(  48 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  38 ]->get_child(  65 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  38 ]->get_child(  66 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  38 ]->get_child(  67 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  38 ]->get_child(  68 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  38 ]->get_child(  70 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  38 ]->get_child(  71 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  38 ]->get_child(  72 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  38 ]->get_child(  73 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  38 ]->get_child(  90 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  38 ]->get_child(  91 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  38 ]->get_child(  92 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  38 ]->get_child(  93 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  38 ]->get_child(  95 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  38 ]->get_child(  96 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  38 ]->get_child(  97 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  38 ]->get_child(  98 );
            }
        }

        if ( mBasis[  39 ] != nullptr )
        {
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  39 ]->get_child(  16 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  39 ]->get_child(  17 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  39 ]->get_child(  18 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  39 ]->get_child(  19 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  39 ]->get_child(  21 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  39 ]->get_child(  22 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  39 ]->get_child(  23 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  39 ]->get_child(  24 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  39 ]->get_child(  41 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  39 ]->get_child(  42 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  39 ]->get_child(  43 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  39 ]->get_child(  44 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  39 ]->get_child(  46 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  39 ]->get_child(  47 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  39 ]->get_child(  48 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  39 ]->get_child(  49 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  39 ]->get_child(  66 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  39 ]->get_child(  67 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  39 ]->get_child(  68 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  39 ]->get_child(  69 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  39 ]->get_child(  71 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  39 ]->get_child(  72 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  39 ]->get_child(  73 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  39 ]->get_child(  74 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  39 ]->get_child(  91 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  39 ]->get_child(  92 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  39 ]->get_child(  93 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  39 ]->get_child(  94 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  39 ]->get_child(  96 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  39 ]->get_child(  97 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  39 ]->get_child(  98 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  39 ]->get_child(  99 );
            }
        }

        if ( mBasis[  40 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  40 ]->get_child(  33 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  40 ]->get_child(  34 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  40 ]->get_child(  38 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  40 ]->get_child(  39 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  40 ]->get_child(  43 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  40 ]->get_child(  44 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  40 ]->get_child(  48 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  40 ]->get_child(  49 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  40 ]->get_child(  58 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  40 ]->get_child(  59 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  40 ]->get_child(  63 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  40 ]->get_child(  64 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  40 ]->get_child(  68 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  40 ]->get_child(  69 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  40 ]->get_child(  73 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  40 ]->get_child(  74 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  40 ]->get_child(  83 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  40 ]->get_child(  84 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  40 ]->get_child(  88 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  40 ]->get_child(  89 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  40 ]->get_child(  93 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  40 ]->get_child(  94 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  40 ]->get_child(  98 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  40 ]->get_child(  99 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  40 ]->get_child( 108 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  40 ]->get_child( 109 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  40 ]->get_child( 113 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  40 ]->get_child( 114 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  40 ]->get_child( 118 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  40 ]->get_child( 119 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  40 ]->get_child( 123 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  40 ]->get_child( 124 );
            }
        }

        if ( mBasis[  41 ] != nullptr )
        {
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  41 ]->get_child(   8 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  41 ]->get_child(   9 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  41 ]->get_child(  13 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  41 ]->get_child(  14 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  41 ]->get_child(  18 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  41 ]->get_child(  19 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  41 ]->get_child(  23 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  41 ]->get_child(  24 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  41 ]->get_child(  33 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  41 ]->get_child(  34 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  41 ]->get_child(  38 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  41 ]->get_child(  39 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  41 ]->get_child(  43 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  41 ]->get_child(  44 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  41 ]->get_child(  48 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  41 ]->get_child(  49 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  41 ]->get_child(  58 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  41 ]->get_child(  59 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  41 ]->get_child(  63 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  41 ]->get_child(  64 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  41 ]->get_child(  68 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  41 ]->get_child(  69 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  41 ]->get_child(  73 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  41 ]->get_child(  74 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  41 ]->get_child(  83 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  41 ]->get_child(  84 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  41 ]->get_child(  88 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  41 ]->get_child(  89 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  41 ]->get_child(  93 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  41 ]->get_child(  94 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  41 ]->get_child(  98 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  41 ]->get_child(  99 );
            }
        }

        if ( mBasis[  42 ] != nullptr )
        {
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  42 ]->get_child(   3 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  42 ]->get_child(   4 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  42 ]->get_child(   8 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  42 ]->get_child(   9 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  42 ]->get_child(  13 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  42 ]->get_child(  14 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  42 ]->get_child(  18 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  42 ]->get_child(  19 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  42 ]->get_child(  28 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  42 ]->get_child(  29 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  42 ]->get_child(  33 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  42 ]->get_child(  34 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  42 ]->get_child(  38 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  42 ]->get_child(  39 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  42 ]->get_child(  43 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  42 ]->get_child(  44 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  42 ]->get_child(  53 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  42 ]->get_child(  54 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  42 ]->get_child(  58 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  42 ]->get_child(  59 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  42 ]->get_child(  63 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  42 ]->get_child(  64 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  42 ]->get_child(  68 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  42 ]->get_child(  69 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  42 ]->get_child(  78 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  42 ]->get_child(  79 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  42 ]->get_child(  83 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  42 ]->get_child(  84 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  42 ]->get_child(  88 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  42 ]->get_child(  89 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  42 ]->get_child(  93 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  42 ]->get_child(  94 );
            }
        }

        if ( mBasis[  43 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  43 ]->get_child(  28 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  43 ]->get_child(  29 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  43 ]->get_child(  33 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  43 ]->get_child(  34 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  43 ]->get_child(  38 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  43 ]->get_child(  39 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  43 ]->get_child(  43 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  43 ]->get_child(  44 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  43 ]->get_child(  53 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  43 ]->get_child(  54 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  43 ]->get_child(  58 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  43 ]->get_child(  59 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  43 ]->get_child(  63 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  43 ]->get_child(  64 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  43 ]->get_child(  68 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  43 ]->get_child(  69 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  43 ]->get_child(  78 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  43 ]->get_child(  79 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  43 ]->get_child(  83 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  43 ]->get_child(  84 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  43 ]->get_child(  88 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  43 ]->get_child(  89 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  43 ]->get_child(  93 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  43 ]->get_child(  94 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  43 ]->get_child( 103 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  43 ]->get_child( 104 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  43 ]->get_child( 108 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  43 ]->get_child( 109 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  43 ]->get_child( 113 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  43 ]->get_child( 114 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  43 ]->get_child( 118 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  43 ]->get_child( 119 );
            }
        }

        if ( mBasis[  44 ] != nullptr )
        {
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  44 ]->get_child(  30 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  44 ]->get_child(  31 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  44 ]->get_child(  35 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  44 ]->get_child(  36 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  44 ]->get_child(  40 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  44 ]->get_child(  41 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  44 ]->get_child(  45 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  44 ]->get_child(  46 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  44 ]->get_child(  55 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  44 ]->get_child(  56 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  44 ]->get_child(  60 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  44 ]->get_child(  61 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  44 ]->get_child(  65 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  44 ]->get_child(  66 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  44 ]->get_child(  70 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  44 ]->get_child(  71 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  44 ]->get_child(  80 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  44 ]->get_child(  81 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  44 ]->get_child(  85 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  44 ]->get_child(  86 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  44 ]->get_child(  90 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  44 ]->get_child(  91 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  44 ]->get_child(  95 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  44 ]->get_child(  96 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  44 ]->get_child( 105 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  44 ]->get_child( 106 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  44 ]->get_child( 110 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  44 ]->get_child( 111 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  44 ]->get_child( 115 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  44 ]->get_child( 116 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  44 ]->get_child( 120 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  44 ]->get_child( 121 );
            }
        }

        if ( mBasis[  45 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  45 ]->get_child(  25 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  45 ]->get_child(  26 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  45 ]->get_child(  30 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  45 ]->get_child(  31 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  45 ]->get_child(  35 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  45 ]->get_child(  36 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  45 ]->get_child(  40 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  45 ]->get_child(  41 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  45 ]->get_child(  50 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  45 ]->get_child(  51 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  45 ]->get_child(  55 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  45 ]->get_child(  56 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  45 ]->get_child(  60 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  45 ]->get_child(  61 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  45 ]->get_child(  65 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  45 ]->get_child(  66 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  45 ]->get_child(  75 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  45 ]->get_child(  76 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  45 ]->get_child(  80 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  45 ]->get_child(  81 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  45 ]->get_child(  85 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  45 ]->get_child(  86 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  45 ]->get_child(  90 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  45 ]->get_child(  91 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  45 ]->get_child( 100 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  45 ]->get_child( 101 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  45 ]->get_child( 105 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  45 ]->get_child( 106 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  45 ]->get_child( 110 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  45 ]->get_child( 111 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  45 ]->get_child( 115 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  45 ]->get_child( 116 );
            }
        }

        if ( mBasis[  46 ] != nullptr )
        {
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  46 ]->get_child(   0 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  46 ]->get_child(   1 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  46 ]->get_child(   5 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  46 ]->get_child(   6 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  46 ]->get_child(  10 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  46 ]->get_child(  11 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  46 ]->get_child(  15 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  46 ]->get_child(  16 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  46 ]->get_child(  25 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  46 ]->get_child(  26 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  46 ]->get_child(  30 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  46 ]->get_child(  31 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  46 ]->get_child(  35 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  46 ]->get_child(  36 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  46 ]->get_child(  40 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  46 ]->get_child(  41 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  46 ]->get_child(  50 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  46 ]->get_child(  51 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  46 ]->get_child(  55 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  46 ]->get_child(  56 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  46 ]->get_child(  60 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  46 ]->get_child(  61 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  46 ]->get_child(  65 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  46 ]->get_child(  66 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  46 ]->get_child(  75 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  46 ]->get_child(  76 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  46 ]->get_child(  80 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  46 ]->get_child(  81 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  46 ]->get_child(  85 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  46 ]->get_child(  86 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  46 ]->get_child(  90 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  46 ]->get_child(  91 );
            }
        }

        if ( mBasis[  47 ] != nullptr )
        {
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  47 ]->get_child(   5 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  47 ]->get_child(   6 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  47 ]->get_child(  10 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  47 ]->get_child(  11 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  47 ]->get_child(  15 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  47 ]->get_child(  16 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  47 ]->get_child(  20 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  47 ]->get_child(  21 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  47 ]->get_child(  30 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  47 ]->get_child(  31 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  47 ]->get_child(  35 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  47 ]->get_child(  36 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  47 ]->get_child(  40 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  47 ]->get_child(  41 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  47 ]->get_child(  45 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  47 ]->get_child(  46 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  47 ]->get_child(  55 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  47 ]->get_child(  56 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  47 ]->get_child(  60 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  47 ]->get_child(  61 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  47 ]->get_child(  65 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  47 ]->get_child(  66 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  47 ]->get_child(  70 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  47 ]->get_child(  71 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  47 ]->get_child(  80 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  47 ]->get_child(  81 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  47 ]->get_child(  85 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  47 ]->get_child(  86 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  47 ]->get_child(  90 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  47 ]->get_child(  91 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  47 ]->get_child(  95 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  47 ]->get_child(  96 );
            }
        }

        if ( mBasis[  48 ] != nullptr )
        {
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  48 ]->get_child(  25 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  48 ]->get_child(  26 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  48 ]->get_child(  27 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  48 ]->get_child(  28 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  48 ]->get_child(  30 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  48 ]->get_child(  31 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  48 ]->get_child(  32 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  48 ]->get_child(  33 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  48 ]->get_child(  50 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  48 ]->get_child(  51 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  48 ]->get_child(  52 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  48 ]->get_child(  53 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  48 ]->get_child(  55 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  48 ]->get_child(  56 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  48 ]->get_child(  57 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  48 ]->get_child(  58 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  48 ]->get_child(  75 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  48 ]->get_child(  76 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  48 ]->get_child(  77 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  48 ]->get_child(  78 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  48 ]->get_child(  80 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  48 ]->get_child(  81 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  48 ]->get_child(  82 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  48 ]->get_child(  83 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  48 ]->get_child( 100 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  48 ]->get_child( 101 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  48 ]->get_child( 102 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  48 ]->get_child( 103 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  48 ]->get_child( 105 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  48 ]->get_child( 106 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  48 ]->get_child( 107 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  48 ]->get_child( 108 );
            }
        }

        if ( mBasis[  49 ] != nullptr )
        {
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  49 ]->get_child(  26 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  49 ]->get_child(  27 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  49 ]->get_child(  28 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  49 ]->get_child(  29 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  49 ]->get_child(  31 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  49 ]->get_child(  32 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  49 ]->get_child(  33 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  49 ]->get_child(  34 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  49 ]->get_child(  51 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  49 ]->get_child(  52 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  49 ]->get_child(  53 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  49 ]->get_child(  54 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  49 ]->get_child(  56 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  49 ]->get_child(  57 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  49 ]->get_child(  58 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  49 ]->get_child(  59 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  49 ]->get_child(  76 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  49 ]->get_child(  77 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  49 ]->get_child(  78 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  49 ]->get_child(  79 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  49 ]->get_child(  81 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  49 ]->get_child(  82 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  49 ]->get_child(  83 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  49 ]->get_child(  84 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  49 ]->get_child( 101 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  49 ]->get_child( 102 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  49 ]->get_child( 103 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  49 ]->get_child( 104 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  49 ]->get_child( 106 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  49 ]->get_child( 107 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  49 ]->get_child( 108 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  49 ]->get_child( 109 );
            }
        }

        if ( mBasis[  50 ] != nullptr )
        {
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  50 ]->get_child(   1 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  50 ]->get_child(   2 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  50 ]->get_child(   3 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  50 ]->get_child(   4 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  50 ]->get_child(   6 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  50 ]->get_child(   7 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  50 ]->get_child(   8 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  50 ]->get_child(   9 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  50 ]->get_child(  26 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  50 ]->get_child(  27 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  50 ]->get_child(  28 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  50 ]->get_child(  29 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  50 ]->get_child(  31 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  50 ]->get_child(  32 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  50 ]->get_child(  33 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  50 ]->get_child(  34 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  50 ]->get_child(  51 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  50 ]->get_child(  52 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  50 ]->get_child(  53 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  50 ]->get_child(  54 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  50 ]->get_child(  56 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  50 ]->get_child(  57 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  50 ]->get_child(  58 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  50 ]->get_child(  59 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  50 ]->get_child(  76 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  50 ]->get_child(  77 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  50 ]->get_child(  78 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  50 ]->get_child(  79 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  50 ]->get_child(  81 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  50 ]->get_child(  82 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  50 ]->get_child(  83 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  50 ]->get_child(  84 );
            }
        }

        if ( mBasis[  51 ] != nullptr )
        {
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  51 ]->get_child(   0 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  51 ]->get_child(   1 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  51 ]->get_child(   2 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  51 ]->get_child(   3 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  51 ]->get_child(   5 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  51 ]->get_child(   6 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  51 ]->get_child(   7 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  51 ]->get_child(   8 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  51 ]->get_child(  25 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  51 ]->get_child(  26 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  51 ]->get_child(  27 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  51 ]->get_child(  28 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  51 ]->get_child(  30 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  51 ]->get_child(  31 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  51 ]->get_child(  32 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  51 ]->get_child(  33 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  51 ]->get_child(  50 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  51 ]->get_child(  51 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  51 ]->get_child(  52 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  51 ]->get_child(  53 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  51 ]->get_child(  55 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  51 ]->get_child(  56 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  51 ]->get_child(  57 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  51 ]->get_child(  58 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  51 ]->get_child(  75 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  51 ]->get_child(  76 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  51 ]->get_child(  77 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  51 ]->get_child(  78 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  51 ]->get_child(  80 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  51 ]->get_child(  81 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  51 ]->get_child(  82 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  51 ]->get_child(  83 );
            }
        }

        if ( mBasis[  52 ] != nullptr )
        {
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  52 ]->get_child(   6 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  52 ]->get_child(   7 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  52 ]->get_child(   8 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  52 ]->get_child(   9 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  52 ]->get_child(  11 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  52 ]->get_child(  12 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  52 ]->get_child(  13 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  52 ]->get_child(  14 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  52 ]->get_child(  16 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  52 ]->get_child(  17 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  52 ]->get_child(  18 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  52 ]->get_child(  19 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  52 ]->get_child(  21 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  52 ]->get_child(  22 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  52 ]->get_child(  23 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  52 ]->get_child(  24 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  52 ]->get_child(  31 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  52 ]->get_child(  32 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  52 ]->get_child(  33 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  52 ]->get_child(  34 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  52 ]->get_child(  36 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  52 ]->get_child(  37 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  52 ]->get_child(  38 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  52 ]->get_child(  39 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  52 ]->get_child(  41 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  52 ]->get_child(  42 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  52 ]->get_child(  43 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  52 ]->get_child(  44 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  52 ]->get_child(  46 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  52 ]->get_child(  47 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  52 ]->get_child(  48 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  52 ]->get_child(  49 );
            }
        }

        if ( mBasis[  53 ] != nullptr )
        {
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  53 ]->get_child(   5 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  53 ]->get_child(   6 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  53 ]->get_child(   7 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  53 ]->get_child(   8 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  53 ]->get_child(  10 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  53 ]->get_child(  11 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  53 ]->get_child(  12 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  53 ]->get_child(  13 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  53 ]->get_child(  15 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  53 ]->get_child(  16 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  53 ]->get_child(  17 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  53 ]->get_child(  18 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  53 ]->get_child(  20 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  53 ]->get_child(  21 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  53 ]->get_child(  22 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  53 ]->get_child(  23 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  53 ]->get_child(  30 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  53 ]->get_child(  31 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  53 ]->get_child(  32 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  53 ]->get_child(  33 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  53 ]->get_child(  35 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  53 ]->get_child(  36 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  53 ]->get_child(  37 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  53 ]->get_child(  38 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  53 ]->get_child(  40 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  53 ]->get_child(  41 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  53 ]->get_child(  42 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  53 ]->get_child(  43 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  53 ]->get_child(  45 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  53 ]->get_child(  46 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  53 ]->get_child(  47 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  53 ]->get_child(  48 );
            }
        }

        if ( mBasis[  54 ] != nullptr )
        {
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  54 ]->get_child(   0 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  54 ]->get_child(   1 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  54 ]->get_child(   2 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  54 ]->get_child(   3 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  54 ]->get_child(   5 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  54 ]->get_child(   6 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  54 ]->get_child(   7 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  54 ]->get_child(   8 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  54 ]->get_child(  10 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  54 ]->get_child(  11 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  54 ]->get_child(  12 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  54 ]->get_child(  13 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  54 ]->get_child(  15 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  54 ]->get_child(  16 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  54 ]->get_child(  17 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  54 ]->get_child(  18 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  54 ]->get_child(  25 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  54 ]->get_child(  26 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  54 ]->get_child(  27 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  54 ]->get_child(  28 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  54 ]->get_child(  30 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  54 ]->get_child(  31 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  54 ]->get_child(  32 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  54 ]->get_child(  33 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  54 ]->get_child(  35 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  54 ]->get_child(  36 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  54 ]->get_child(  37 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  54 ]->get_child(  38 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  54 ]->get_child(  40 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  54 ]->get_child(  41 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  54 ]->get_child(  42 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  54 ]->get_child(  43 );
            }
        }

        if ( mBasis[  55 ] != nullptr )
        {
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  55 ]->get_child(   1 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  55 ]->get_child(   2 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  55 ]->get_child(   3 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  55 ]->get_child(   4 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  55 ]->get_child(   6 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  55 ]->get_child(   7 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  55 ]->get_child(   8 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  55 ]->get_child(   9 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  55 ]->get_child(  11 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  55 ]->get_child(  12 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  55 ]->get_child(  13 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  55 ]->get_child(  14 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  55 ]->get_child(  16 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  55 ]->get_child(  17 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  55 ]->get_child(  18 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  55 ]->get_child(  19 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  55 ]->get_child(  26 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  55 ]->get_child(  27 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  55 ]->get_child(  28 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  55 ]->get_child(  29 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  55 ]->get_child(  31 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  55 ]->get_child(  32 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  55 ]->get_child(  33 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  55 ]->get_child(  34 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  55 ]->get_child(  36 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  55 ]->get_child(  37 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  55 ]->get_child(  38 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  55 ]->get_child(  39 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  55 ]->get_child(  41 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  55 ]->get_child(  42 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  55 ]->get_child(  43 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  55 ]->get_child(  44 );
            }
        }

        if ( mBasis[  56 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[  56 ]->get_child(  31 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  56 ]->get_child(  32 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  56 ]->get_child(  33 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  56 ]->get_child(  34 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  56 ]->get_child(  36 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  56 ]->get_child(  37 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  56 ]->get_child(  38 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  56 ]->get_child(  39 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  56 ]->get_child(  41 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  56 ]->get_child(  42 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  56 ]->get_child(  43 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  56 ]->get_child(  44 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  56 ]->get_child(  46 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  56 ]->get_child(  47 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  56 ]->get_child(  48 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  56 ]->get_child(  49 );
            }
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  56 ]->get_child(  56 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  56 ]->get_child(  57 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  56 ]->get_child(  58 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  56 ]->get_child(  59 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  56 ]->get_child(  61 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  56 ]->get_child(  62 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  56 ]->get_child(  63 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  56 ]->get_child(  64 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  56 ]->get_child(  66 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  56 ]->get_child(  67 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  56 ]->get_child(  68 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  56 ]->get_child(  69 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  56 ]->get_child(  71 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  56 ]->get_child(  72 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  56 ]->get_child(  73 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  56 ]->get_child(  74 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  56 ]->get_child(  81 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  56 ]->get_child(  82 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  56 ]->get_child(  83 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  56 ]->get_child(  84 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  56 ]->get_child(  86 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  56 ]->get_child(  87 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  56 ]->get_child(  88 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  56 ]->get_child(  89 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  56 ]->get_child(  91 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  56 ]->get_child(  92 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  56 ]->get_child(  93 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  56 ]->get_child(  94 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  56 ]->get_child(  96 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  56 ]->get_child(  97 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  56 ]->get_child(  98 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  56 ]->get_child(  99 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  56 ]->get_child( 106 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  56 ]->get_child( 107 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  56 ]->get_child( 108 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  56 ]->get_child( 109 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  56 ]->get_child( 111 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  56 ]->get_child( 112 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  56 ]->get_child( 113 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  56 ]->get_child( 114 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  56 ]->get_child( 116 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  56 ]->get_child( 117 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  56 ]->get_child( 118 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  56 ]->get_child( 119 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  56 ]->get_child( 121 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  56 ]->get_child( 122 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  56 ]->get_child( 123 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  56 ]->get_child( 124 );
            }
        }

        if ( mBasis[  57 ] != nullptr )
        {
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[  57 ]->get_child(  30 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[  57 ]->get_child(  31 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[  57 ]->get_child(  32 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[  57 ]->get_child(  33 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  57 ]->get_child(  35 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  57 ]->get_child(  36 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  57 ]->get_child(  37 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  57 ]->get_child(  38 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  57 ]->get_child(  40 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  57 ]->get_child(  41 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  57 ]->get_child(  42 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  57 ]->get_child(  43 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  57 ]->get_child(  45 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  57 ]->get_child(  46 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  57 ]->get_child(  47 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  57 ]->get_child(  48 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  57 ]->get_child(  55 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  57 ]->get_child(  56 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  57 ]->get_child(  57 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  57 ]->get_child(  58 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  57 ]->get_child(  60 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  57 ]->get_child(  61 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  57 ]->get_child(  62 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  57 ]->get_child(  63 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  57 ]->get_child(  65 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  57 ]->get_child(  66 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  57 ]->get_child(  67 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  57 ]->get_child(  68 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  57 ]->get_child(  70 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  57 ]->get_child(  71 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  57 ]->get_child(  72 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  57 ]->get_child(  73 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  57 ]->get_child(  80 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  57 ]->get_child(  81 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  57 ]->get_child(  82 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  57 ]->get_child(  83 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  57 ]->get_child(  85 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  57 ]->get_child(  86 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  57 ]->get_child(  87 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  57 ]->get_child(  88 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  57 ]->get_child(  90 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  57 ]->get_child(  91 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  57 ]->get_child(  92 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  57 ]->get_child(  93 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  57 ]->get_child(  95 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  57 ]->get_child(  96 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  57 ]->get_child(  97 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  57 ]->get_child(  98 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  57 ]->get_child( 105 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  57 ]->get_child( 106 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  57 ]->get_child( 107 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  57 ]->get_child( 108 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  57 ]->get_child( 110 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  57 ]->get_child( 111 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  57 ]->get_child( 112 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  57 ]->get_child( 113 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  57 ]->get_child( 115 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  57 ]->get_child( 116 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  57 ]->get_child( 117 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  57 ]->get_child( 118 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  57 ]->get_child( 120 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  57 ]->get_child( 121 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  57 ]->get_child( 122 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  57 ]->get_child( 123 );
            }
        }

        if ( mBasis[  58 ] != nullptr )
        {
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  58 ]->get_child(  25 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  58 ]->get_child(  26 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  58 ]->get_child(  27 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[  58 ]->get_child(  28 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  58 ]->get_child(  30 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  58 ]->get_child(  31 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  58 ]->get_child(  32 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[  58 ]->get_child(  33 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  58 ]->get_child(  35 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  58 ]->get_child(  36 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  58 ]->get_child(  37 );
            }
            if ( tBasis[  19 ] == nullptr )
            {
                tBasis[  19 ] = mBasis[  58 ]->get_child(  38 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  58 ]->get_child(  40 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  58 ]->get_child(  41 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  58 ]->get_child(  42 );
            }
            if ( tBasis[  24 ] == nullptr )
            {
                tBasis[  24 ] = mBasis[  58 ]->get_child(  43 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  58 ]->get_child(  50 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  58 ]->get_child(  51 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  58 ]->get_child(  52 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  58 ]->get_child(  53 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  58 ]->get_child(  55 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  58 ]->get_child(  56 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  58 ]->get_child(  57 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  58 ]->get_child(  58 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  58 ]->get_child(  60 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  58 ]->get_child(  61 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  58 ]->get_child(  62 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  58 ]->get_child(  63 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  58 ]->get_child(  65 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  58 ]->get_child(  66 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  58 ]->get_child(  67 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  58 ]->get_child(  68 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  58 ]->get_child(  75 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  58 ]->get_child(  76 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  58 ]->get_child(  77 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  58 ]->get_child(  78 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  58 ]->get_child(  80 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  58 ]->get_child(  81 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  58 ]->get_child(  82 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  58 ]->get_child(  83 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  58 ]->get_child(  85 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  58 ]->get_child(  86 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  58 ]->get_child(  87 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  58 ]->get_child(  88 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  58 ]->get_child(  90 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  58 ]->get_child(  91 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  58 ]->get_child(  92 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  58 ]->get_child(  93 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  58 ]->get_child( 100 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  58 ]->get_child( 101 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  58 ]->get_child( 102 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  58 ]->get_child( 103 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  58 ]->get_child( 105 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  58 ]->get_child( 106 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  58 ]->get_child( 107 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  58 ]->get_child( 108 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  58 ]->get_child( 110 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  58 ]->get_child( 111 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  58 ]->get_child( 112 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  58 ]->get_child( 113 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  58 ]->get_child( 115 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  58 ]->get_child( 116 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  58 ]->get_child( 117 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  58 ]->get_child( 118 );
            }
        }

        if ( mBasis[  59 ] != nullptr )
        {
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[  59 ]->get_child(  26 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[  59 ]->get_child(  27 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[  59 ]->get_child(  28 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[  59 ]->get_child(  29 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[  59 ]->get_child(  31 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[  59 ]->get_child(  32 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[  59 ]->get_child(  33 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[  59 ]->get_child(  34 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[  59 ]->get_child(  36 );
            }
            if ( tBasis[  16 ] == nullptr )
            {
                tBasis[  16 ] = mBasis[  59 ]->get_child(  37 );
            }
            if ( tBasis[  17 ] == nullptr )
            {
                tBasis[  17 ] = mBasis[  59 ]->get_child(  38 );
            }
            if ( tBasis[  18 ] == nullptr )
            {
                tBasis[  18 ] = mBasis[  59 ]->get_child(  39 );
            }
            if ( tBasis[  20 ] == nullptr )
            {
                tBasis[  20 ] = mBasis[  59 ]->get_child(  41 );
            }
            if ( tBasis[  21 ] == nullptr )
            {
                tBasis[  21 ] = mBasis[  59 ]->get_child(  42 );
            }
            if ( tBasis[  22 ] == nullptr )
            {
                tBasis[  22 ] = mBasis[  59 ]->get_child(  43 );
            }
            if ( tBasis[  23 ] == nullptr )
            {
                tBasis[  23 ] = mBasis[  59 ]->get_child(  44 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  59 ]->get_child(  51 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  59 ]->get_child(  52 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  59 ]->get_child(  53 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  59 ]->get_child(  54 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  59 ]->get_child(  56 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  59 ]->get_child(  57 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  59 ]->get_child(  58 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  59 ]->get_child(  59 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  59 ]->get_child(  61 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  59 ]->get_child(  62 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  59 ]->get_child(  63 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  59 ]->get_child(  64 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  59 ]->get_child(  66 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  59 ]->get_child(  67 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  59 ]->get_child(  68 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  59 ]->get_child(  69 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  59 ]->get_child(  76 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  59 ]->get_child(  77 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  59 ]->get_child(  78 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  59 ]->get_child(  79 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  59 ]->get_child(  81 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  59 ]->get_child(  82 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  59 ]->get_child(  83 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  59 ]->get_child(  84 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  59 ]->get_child(  86 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  59 ]->get_child(  87 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  59 ]->get_child(  88 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  59 ]->get_child(  89 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  59 ]->get_child(  91 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  59 ]->get_child(  92 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  59 ]->get_child(  93 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  59 ]->get_child(  94 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  59 ]->get_child( 101 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  59 ]->get_child( 102 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  59 ]->get_child( 103 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  59 ]->get_child( 104 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  59 ]->get_child( 106 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  59 ]->get_child( 107 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  59 ]->get_child( 108 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  59 ]->get_child( 109 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  59 ]->get_child( 111 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  59 ]->get_child( 112 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  59 ]->get_child( 113 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  59 ]->get_child( 114 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  59 ]->get_child( 116 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  59 ]->get_child( 117 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  59 ]->get_child( 118 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  59 ]->get_child( 119 );
            }
        }

        if ( mBasis[  60 ] != nullptr )
        {
            if ( tBasis[  25 ] == nullptr )
            {
                tBasis[  25 ] = mBasis[  60 ]->get_child(   6 );
            }
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  60 ]->get_child(   7 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  60 ]->get_child(   8 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  60 ]->get_child(   9 );
            }
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  60 ]->get_child(  11 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  60 ]->get_child(  12 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  60 ]->get_child(  13 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  60 ]->get_child(  14 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  60 ]->get_child(  16 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  60 ]->get_child(  17 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  60 ]->get_child(  18 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  60 ]->get_child(  19 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  60 ]->get_child(  21 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  60 ]->get_child(  22 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  60 ]->get_child(  23 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  60 ]->get_child(  24 );
            }
            if ( tBasis[  50 ] == nullptr )
            {
                tBasis[  50 ] = mBasis[  60 ]->get_child(  31 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  60 ]->get_child(  32 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  60 ]->get_child(  33 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  60 ]->get_child(  34 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  60 ]->get_child(  36 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  60 ]->get_child(  37 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  60 ]->get_child(  38 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  60 ]->get_child(  39 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  60 ]->get_child(  41 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  60 ]->get_child(  42 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  60 ]->get_child(  43 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  60 ]->get_child(  44 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  60 ]->get_child(  46 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  60 ]->get_child(  47 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  60 ]->get_child(  48 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  60 ]->get_child(  49 );
            }
            if ( tBasis[  75 ] == nullptr )
            {
                tBasis[  75 ] = mBasis[  60 ]->get_child(  56 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  60 ]->get_child(  57 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  60 ]->get_child(  58 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  60 ]->get_child(  59 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  60 ]->get_child(  61 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  60 ]->get_child(  62 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  60 ]->get_child(  63 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  60 ]->get_child(  64 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  60 ]->get_child(  66 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  60 ]->get_child(  67 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  60 ]->get_child(  68 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  60 ]->get_child(  69 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  60 ]->get_child(  71 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  60 ]->get_child(  72 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  60 ]->get_child(  73 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  60 ]->get_child(  74 );
            }
            if ( tBasis[ 100 ] == nullptr )
            {
                tBasis[ 100 ] = mBasis[  60 ]->get_child(  81 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  60 ]->get_child(  82 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  60 ]->get_child(  83 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  60 ]->get_child(  84 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  60 ]->get_child(  86 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  60 ]->get_child(  87 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  60 ]->get_child(  88 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  60 ]->get_child(  89 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  60 ]->get_child(  91 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  60 ]->get_child(  92 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  60 ]->get_child(  93 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  60 ]->get_child(  94 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  60 ]->get_child(  96 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  60 ]->get_child(  97 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  60 ]->get_child(  98 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  60 ]->get_child(  99 );
            }
        }

        if ( mBasis[  61 ] != nullptr )
        {
            if ( tBasis[  26 ] == nullptr )
            {
                tBasis[  26 ] = mBasis[  61 ]->get_child(   5 );
            }
            if ( tBasis[  27 ] == nullptr )
            {
                tBasis[  27 ] = mBasis[  61 ]->get_child(   6 );
            }
            if ( tBasis[  28 ] == nullptr )
            {
                tBasis[  28 ] = mBasis[  61 ]->get_child(   7 );
            }
            if ( tBasis[  29 ] == nullptr )
            {
                tBasis[  29 ] = mBasis[  61 ]->get_child(   8 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  61 ]->get_child(  10 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  61 ]->get_child(  11 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  61 ]->get_child(  12 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  61 ]->get_child(  13 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  61 ]->get_child(  15 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  61 ]->get_child(  16 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  61 ]->get_child(  17 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  61 ]->get_child(  18 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  61 ]->get_child(  20 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  61 ]->get_child(  21 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  61 ]->get_child(  22 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  61 ]->get_child(  23 );
            }
            if ( tBasis[  51 ] == nullptr )
            {
                tBasis[  51 ] = mBasis[  61 ]->get_child(  30 );
            }
            if ( tBasis[  52 ] == nullptr )
            {
                tBasis[  52 ] = mBasis[  61 ]->get_child(  31 );
            }
            if ( tBasis[  53 ] == nullptr )
            {
                tBasis[  53 ] = mBasis[  61 ]->get_child(  32 );
            }
            if ( tBasis[  54 ] == nullptr )
            {
                tBasis[  54 ] = mBasis[  61 ]->get_child(  33 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  61 ]->get_child(  35 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  61 ]->get_child(  36 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  61 ]->get_child(  37 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  61 ]->get_child(  38 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  61 ]->get_child(  40 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  61 ]->get_child(  41 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  61 ]->get_child(  42 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  61 ]->get_child(  43 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  61 ]->get_child(  45 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  61 ]->get_child(  46 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  61 ]->get_child(  47 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  61 ]->get_child(  48 );
            }
            if ( tBasis[  76 ] == nullptr )
            {
                tBasis[  76 ] = mBasis[  61 ]->get_child(  55 );
            }
            if ( tBasis[  77 ] == nullptr )
            {
                tBasis[  77 ] = mBasis[  61 ]->get_child(  56 );
            }
            if ( tBasis[  78 ] == nullptr )
            {
                tBasis[  78 ] = mBasis[  61 ]->get_child(  57 );
            }
            if ( tBasis[  79 ] == nullptr )
            {
                tBasis[  79 ] = mBasis[  61 ]->get_child(  58 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  61 ]->get_child(  60 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  61 ]->get_child(  61 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  61 ]->get_child(  62 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  61 ]->get_child(  63 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  61 ]->get_child(  65 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  61 ]->get_child(  66 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  61 ]->get_child(  67 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  61 ]->get_child(  68 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  61 ]->get_child(  70 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  61 ]->get_child(  71 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  61 ]->get_child(  72 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  61 ]->get_child(  73 );
            }
            if ( tBasis[ 101 ] == nullptr )
            {
                tBasis[ 101 ] = mBasis[  61 ]->get_child(  80 );
            }
            if ( tBasis[ 102 ] == nullptr )
            {
                tBasis[ 102 ] = mBasis[  61 ]->get_child(  81 );
            }
            if ( tBasis[ 103 ] == nullptr )
            {
                tBasis[ 103 ] = mBasis[  61 ]->get_child(  82 );
            }
            if ( tBasis[ 104 ] == nullptr )
            {
                tBasis[ 104 ] = mBasis[  61 ]->get_child(  83 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  61 ]->get_child(  85 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  61 ]->get_child(  86 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  61 ]->get_child(  87 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  61 ]->get_child(  88 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  61 ]->get_child(  90 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  61 ]->get_child(  91 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  61 ]->get_child(  92 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  61 ]->get_child(  93 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  61 ]->get_child(  95 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  61 ]->get_child(  96 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  61 ]->get_child(  97 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  61 ]->get_child(  98 );
            }
        }

        if ( mBasis[  62 ] != nullptr )
        {
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  62 ]->get_child(   0 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  62 ]->get_child(   1 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  62 ]->get_child(   2 );
            }
            if ( tBasis[  34 ] == nullptr )
            {
                tBasis[  34 ] = mBasis[  62 ]->get_child(   3 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  62 ]->get_child(   5 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  62 ]->get_child(   6 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  62 ]->get_child(   7 );
            }
            if ( tBasis[  39 ] == nullptr )
            {
                tBasis[  39 ] = mBasis[  62 ]->get_child(   8 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  62 ]->get_child(  10 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  62 ]->get_child(  11 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  62 ]->get_child(  12 );
            }
            if ( tBasis[  44 ] == nullptr )
            {
                tBasis[  44 ] = mBasis[  62 ]->get_child(  13 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  62 ]->get_child(  15 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  62 ]->get_child(  16 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  62 ]->get_child(  17 );
            }
            if ( tBasis[  49 ] == nullptr )
            {
                tBasis[  49 ] = mBasis[  62 ]->get_child(  18 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  62 ]->get_child(  25 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  62 ]->get_child(  26 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  62 ]->get_child(  27 );
            }
            if ( tBasis[  59 ] == nullptr )
            {
                tBasis[  59 ] = mBasis[  62 ]->get_child(  28 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  62 ]->get_child(  30 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  62 ]->get_child(  31 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  62 ]->get_child(  32 );
            }
            if ( tBasis[  64 ] == nullptr )
            {
                tBasis[  64 ] = mBasis[  62 ]->get_child(  33 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  62 ]->get_child(  35 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  62 ]->get_child(  36 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  62 ]->get_child(  37 );
            }
            if ( tBasis[  69 ] == nullptr )
            {
                tBasis[  69 ] = mBasis[  62 ]->get_child(  38 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  62 ]->get_child(  40 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  62 ]->get_child(  41 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  62 ]->get_child(  42 );
            }
            if ( tBasis[  74 ] == nullptr )
            {
                tBasis[  74 ] = mBasis[  62 ]->get_child(  43 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  62 ]->get_child(  50 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  62 ]->get_child(  51 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  62 ]->get_child(  52 );
            }
            if ( tBasis[  84 ] == nullptr )
            {
                tBasis[  84 ] = mBasis[  62 ]->get_child(  53 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  62 ]->get_child(  55 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  62 ]->get_child(  56 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  62 ]->get_child(  57 );
            }
            if ( tBasis[  89 ] == nullptr )
            {
                tBasis[  89 ] = mBasis[  62 ]->get_child(  58 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  62 ]->get_child(  60 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  62 ]->get_child(  61 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  62 ]->get_child(  62 );
            }
            if ( tBasis[  94 ] == nullptr )
            {
                tBasis[  94 ] = mBasis[  62 ]->get_child(  63 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  62 ]->get_child(  65 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  62 ]->get_child(  66 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  62 ]->get_child(  67 );
            }
            if ( tBasis[  99 ] == nullptr )
            {
                tBasis[  99 ] = mBasis[  62 ]->get_child(  68 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  62 ]->get_child(  75 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  62 ]->get_child(  76 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  62 ]->get_child(  77 );
            }
            if ( tBasis[ 109 ] == nullptr )
            {
                tBasis[ 109 ] = mBasis[  62 ]->get_child(  78 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  62 ]->get_child(  80 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  62 ]->get_child(  81 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  62 ]->get_child(  82 );
            }
            if ( tBasis[ 114 ] == nullptr )
            {
                tBasis[ 114 ] = mBasis[  62 ]->get_child(  83 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  62 ]->get_child(  85 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  62 ]->get_child(  86 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  62 ]->get_child(  87 );
            }
            if ( tBasis[ 119 ] == nullptr )
            {
                tBasis[ 119 ] = mBasis[  62 ]->get_child(  88 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  62 ]->get_child(  90 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  62 ]->get_child(  91 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  62 ]->get_child(  92 );
            }
            if ( tBasis[ 124 ] == nullptr )
            {
                tBasis[ 124 ] = mBasis[  62 ]->get_child(  93 );
            }
        }

        if ( mBasis[  63 ] != nullptr )
        {
            if ( tBasis[  30 ] == nullptr )
            {
                tBasis[  30 ] = mBasis[  63 ]->get_child(   1 );
            }
            if ( tBasis[  31 ] == nullptr )
            {
                tBasis[  31 ] = mBasis[  63 ]->get_child(   2 );
            }
            if ( tBasis[  32 ] == nullptr )
            {
                tBasis[  32 ] = mBasis[  63 ]->get_child(   3 );
            }
            if ( tBasis[  33 ] == nullptr )
            {
                tBasis[  33 ] = mBasis[  63 ]->get_child(   4 );
            }
            if ( tBasis[  35 ] == nullptr )
            {
                tBasis[  35 ] = mBasis[  63 ]->get_child(   6 );
            }
            if ( tBasis[  36 ] == nullptr )
            {
                tBasis[  36 ] = mBasis[  63 ]->get_child(   7 );
            }
            if ( tBasis[  37 ] == nullptr )
            {
                tBasis[  37 ] = mBasis[  63 ]->get_child(   8 );
            }
            if ( tBasis[  38 ] == nullptr )
            {
                tBasis[  38 ] = mBasis[  63 ]->get_child(   9 );
            }
            if ( tBasis[  40 ] == nullptr )
            {
                tBasis[  40 ] = mBasis[  63 ]->get_child(  11 );
            }
            if ( tBasis[  41 ] == nullptr )
            {
                tBasis[  41 ] = mBasis[  63 ]->get_child(  12 );
            }
            if ( tBasis[  42 ] == nullptr )
            {
                tBasis[  42 ] = mBasis[  63 ]->get_child(  13 );
            }
            if ( tBasis[  43 ] == nullptr )
            {
                tBasis[  43 ] = mBasis[  63 ]->get_child(  14 );
            }
            if ( tBasis[  45 ] == nullptr )
            {
                tBasis[  45 ] = mBasis[  63 ]->get_child(  16 );
            }
            if ( tBasis[  46 ] == nullptr )
            {
                tBasis[  46 ] = mBasis[  63 ]->get_child(  17 );
            }
            if ( tBasis[  47 ] == nullptr )
            {
                tBasis[  47 ] = mBasis[  63 ]->get_child(  18 );
            }
            if ( tBasis[  48 ] == nullptr )
            {
                tBasis[  48 ] = mBasis[  63 ]->get_child(  19 );
            }
            if ( tBasis[  55 ] == nullptr )
            {
                tBasis[  55 ] = mBasis[  63 ]->get_child(  26 );
            }
            if ( tBasis[  56 ] == nullptr )
            {
                tBasis[  56 ] = mBasis[  63 ]->get_child(  27 );
            }
            if ( tBasis[  57 ] == nullptr )
            {
                tBasis[  57 ] = mBasis[  63 ]->get_child(  28 );
            }
            if ( tBasis[  58 ] == nullptr )
            {
                tBasis[  58 ] = mBasis[  63 ]->get_child(  29 );
            }
            if ( tBasis[  60 ] == nullptr )
            {
                tBasis[  60 ] = mBasis[  63 ]->get_child(  31 );
            }
            if ( tBasis[  61 ] == nullptr )
            {
                tBasis[  61 ] = mBasis[  63 ]->get_child(  32 );
            }
            if ( tBasis[  62 ] == nullptr )
            {
                tBasis[  62 ] = mBasis[  63 ]->get_child(  33 );
            }
            if ( tBasis[  63 ] == nullptr )
            {
                tBasis[  63 ] = mBasis[  63 ]->get_child(  34 );
            }
            if ( tBasis[  65 ] == nullptr )
            {
                tBasis[  65 ] = mBasis[  63 ]->get_child(  36 );
            }
            if ( tBasis[  66 ] == nullptr )
            {
                tBasis[  66 ] = mBasis[  63 ]->get_child(  37 );
            }
            if ( tBasis[  67 ] == nullptr )
            {
                tBasis[  67 ] = mBasis[  63 ]->get_child(  38 );
            }
            if ( tBasis[  68 ] == nullptr )
            {
                tBasis[  68 ] = mBasis[  63 ]->get_child(  39 );
            }
            if ( tBasis[  70 ] == nullptr )
            {
                tBasis[  70 ] = mBasis[  63 ]->get_child(  41 );
            }
            if ( tBasis[  71 ] == nullptr )
            {
                tBasis[  71 ] = mBasis[  63 ]->get_child(  42 );
            }
            if ( tBasis[  72 ] == nullptr )
            {
                tBasis[  72 ] = mBasis[  63 ]->get_child(  43 );
            }
            if ( tBasis[  73 ] == nullptr )
            {
                tBasis[  73 ] = mBasis[  63 ]->get_child(  44 );
            }
            if ( tBasis[  80 ] == nullptr )
            {
                tBasis[  80 ] = mBasis[  63 ]->get_child(  51 );
            }
            if ( tBasis[  81 ] == nullptr )
            {
                tBasis[  81 ] = mBasis[  63 ]->get_child(  52 );
            }
            if ( tBasis[  82 ] == nullptr )
            {
                tBasis[  82 ] = mBasis[  63 ]->get_child(  53 );
            }
            if ( tBasis[  83 ] == nullptr )
            {
                tBasis[  83 ] = mBasis[  63 ]->get_child(  54 );
            }
            if ( tBasis[  85 ] == nullptr )
            {
                tBasis[  85 ] = mBasis[  63 ]->get_child(  56 );
            }
            if ( tBasis[  86 ] == nullptr )
            {
                tBasis[  86 ] = mBasis[  63 ]->get_child(  57 );
            }
            if ( tBasis[  87 ] == nullptr )
            {
                tBasis[  87 ] = mBasis[  63 ]->get_child(  58 );
            }
            if ( tBasis[  88 ] == nullptr )
            {
                tBasis[  88 ] = mBasis[  63 ]->get_child(  59 );
            }
            if ( tBasis[  90 ] == nullptr )
            {
                tBasis[  90 ] = mBasis[  63 ]->get_child(  61 );
            }
            if ( tBasis[  91 ] == nullptr )
            {
                tBasis[  91 ] = mBasis[  63 ]->get_child(  62 );
            }
            if ( tBasis[  92 ] == nullptr )
            {
                tBasis[  92 ] = mBasis[  63 ]->get_child(  63 );
            }
            if ( tBasis[  93 ] == nullptr )
            {
                tBasis[  93 ] = mBasis[  63 ]->get_child(  64 );
            }
            if ( tBasis[  95 ] == nullptr )
            {
                tBasis[  95 ] = mBasis[  63 ]->get_child(  66 );
            }
            if ( tBasis[  96 ] == nullptr )
            {
                tBasis[  96 ] = mBasis[  63 ]->get_child(  67 );
            }
            if ( tBasis[  97 ] == nullptr )
            {
                tBasis[  97 ] = mBasis[  63 ]->get_child(  68 );
            }
            if ( tBasis[  98 ] == nullptr )
            {
                tBasis[  98 ] = mBasis[  63 ]->get_child(  69 );
            }
            if ( tBasis[ 105 ] == nullptr )
            {
                tBasis[ 105 ] = mBasis[  63 ]->get_child(  76 );
            }
            if ( tBasis[ 106 ] == nullptr )
            {
                tBasis[ 106 ] = mBasis[  63 ]->get_child(  77 );
            }
            if ( tBasis[ 107 ] == nullptr )
            {
                tBasis[ 107 ] = mBasis[  63 ]->get_child(  78 );
            }
            if ( tBasis[ 108 ] == nullptr )
            {
                tBasis[ 108 ] = mBasis[  63 ]->get_child(  79 );
            }
            if ( tBasis[ 110 ] == nullptr )
            {
                tBasis[ 110 ] = mBasis[  63 ]->get_child(  81 );
            }
            if ( tBasis[ 111 ] == nullptr )
            {
                tBasis[ 111 ] = mBasis[  63 ]->get_child(  82 );
            }
            if ( tBasis[ 112 ] == nullptr )
            {
                tBasis[ 112 ] = mBasis[  63 ]->get_child(  83 );
            }
            if ( tBasis[ 113 ] == nullptr )
            {
                tBasis[ 113 ] = mBasis[  63 ]->get_child(  84 );
            }
            if ( tBasis[ 115 ] == nullptr )
            {
                tBasis[ 115 ] = mBasis[  63 ]->get_child(  86 );
            }
            if ( tBasis[ 116 ] == nullptr )
            {
                tBasis[ 116 ] = mBasis[  63 ]->get_child(  87 );
            }
            if ( tBasis[ 117 ] == nullptr )
            {
                tBasis[ 117 ] = mBasis[  63 ]->get_child(  88 );
            }
            if ( tBasis[ 118 ] == nullptr )
            {
                tBasis[ 118 ] = mBasis[  63 ]->get_child(  89 );
            }
            if ( tBasis[ 120 ] == nullptr )
            {
                tBasis[ 120 ] = mBasis[  63 ]->get_child(  91 );
            }
            if ( tBasis[ 121 ] == nullptr )
            {
                tBasis[ 121 ] = mBasis[  63 ]->get_child(  92 );
            }
            if ( tBasis[ 122 ] == nullptr )
            {
                tBasis[ 122 ] = mBasis[  63 ]->get_child(  93 );
            }
            if ( tBasis[ 123 ] == nullptr )
            {
                tBasis[ 123 ] = mBasis[  63 ]->get_child(  94 );
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
        tChildren[ 0 ]->insert_basis(  1, tBasis[   3 ] );
        tChildren[ 0 ]->insert_basis(  2, tBasis[  18 ] );
        tChildren[ 0 ]->insert_basis(  3, tBasis[  15 ] );
        tChildren[ 0 ]->insert_basis(  4, tBasis[  75 ] );
        tChildren[ 0 ]->insert_basis(  5, tBasis[  78 ] );
        tChildren[ 0 ]->insert_basis(  6, tBasis[  93 ] );
        tChildren[ 0 ]->insert_basis(  7, tBasis[  90 ] );
        tChildren[ 0 ]->insert_basis(  8, tBasis[   1 ] );
        tChildren[ 0 ]->insert_basis(  9, tBasis[   2 ] );
        tChildren[ 0 ]->insert_basis( 10, tBasis[   5 ] );
        tChildren[ 0 ]->insert_basis( 11, tBasis[  10 ] );
        tChildren[ 0 ]->insert_basis( 12, tBasis[  25 ] );
        tChildren[ 0 ]->insert_basis( 13, tBasis[  50 ] );
        tChildren[ 0 ]->insert_basis( 14, tBasis[   8 ] );
        tChildren[ 0 ]->insert_basis( 15, tBasis[  13 ] );
        tChildren[ 0 ]->insert_basis( 16, tBasis[  28 ] );
        tChildren[ 0 ]->insert_basis( 17, tBasis[  53 ] );
        tChildren[ 0 ]->insert_basis( 18, tBasis[  17 ] );
        tChildren[ 0 ]->insert_basis( 19, tBasis[  16 ] );
        tChildren[ 0 ]->insert_basis( 20, tBasis[  43 ] );
        tChildren[ 0 ]->insert_basis( 21, tBasis[  68 ] );
        tChildren[ 0 ]->insert_basis( 22, tBasis[  40 ] );
        tChildren[ 0 ]->insert_basis( 23, tBasis[  65 ] );
        tChildren[ 0 ]->insert_basis( 24, tBasis[  76 ] );
        tChildren[ 0 ]->insert_basis( 25, tBasis[  77 ] );
        tChildren[ 0 ]->insert_basis( 26, tBasis[  80 ] );
        tChildren[ 0 ]->insert_basis( 27, tBasis[  85 ] );
        tChildren[ 0 ]->insert_basis( 28, tBasis[  83 ] );
        tChildren[ 0 ]->insert_basis( 29, tBasis[  88 ] );
        tChildren[ 0 ]->insert_basis( 30, tBasis[  92 ] );
        tChildren[ 0 ]->insert_basis( 31, tBasis[  91 ] );
        tChildren[ 0 ]->insert_basis( 32, tBasis[   6 ] );
        tChildren[ 0 ]->insert_basis( 33, tBasis[  11 ] );
        tChildren[ 0 ]->insert_basis( 34, tBasis[  12 ] );
        tChildren[ 0 ]->insert_basis( 35, tBasis[   7 ] );
        tChildren[ 0 ]->insert_basis( 36, tBasis[  26 ] );
        tChildren[ 0 ]->insert_basis( 37, tBasis[  27 ] );
        tChildren[ 0 ]->insert_basis( 38, tBasis[  52 ] );
        tChildren[ 0 ]->insert_basis( 39, tBasis[  51 ] );
        tChildren[ 0 ]->insert_basis( 40, tBasis[  30 ] );
        tChildren[ 0 ]->insert_basis( 41, tBasis[  55 ] );
        tChildren[ 0 ]->insert_basis( 42, tBasis[  60 ] );
        tChildren[ 0 ]->insert_basis( 43, tBasis[  35 ] );
        tChildren[ 0 ]->insert_basis( 44, tBasis[  33 ] );
        tChildren[ 0 ]->insert_basis( 45, tBasis[  38 ] );
        tChildren[ 0 ]->insert_basis( 46, tBasis[  63 ] );
        tChildren[ 0 ]->insert_basis( 47, tBasis[  58 ] );
        tChildren[ 0 ]->insert_basis( 48, tBasis[  42 ] );
        tChildren[ 0 ]->insert_basis( 49, tBasis[  41 ] );
        tChildren[ 0 ]->insert_basis( 50, tBasis[  66 ] );
        tChildren[ 0 ]->insert_basis( 51, tBasis[  67 ] );
        tChildren[ 0 ]->insert_basis( 52, tBasis[  81 ] );
        tChildren[ 0 ]->insert_basis( 53, tBasis[  82 ] );
        tChildren[ 0 ]->insert_basis( 54, tBasis[  87 ] );
        tChildren[ 0 ]->insert_basis( 55, tBasis[  86 ] );
        tChildren[ 0 ]->insert_basis( 56, tBasis[  31 ] );
        tChildren[ 0 ]->insert_basis( 57, tBasis[  32 ] );
        tChildren[ 0 ]->insert_basis( 58, tBasis[  37 ] );
        tChildren[ 0 ]->insert_basis( 59, tBasis[  36 ] );
        tChildren[ 0 ]->insert_basis( 60, tBasis[  56 ] );
        tChildren[ 0 ]->insert_basis( 61, tBasis[  57 ] );
        tChildren[ 0 ]->insert_basis( 62, tBasis[  62 ] );
        tChildren[ 0 ]->insert_basis( 63, tBasis[  61 ] );

        // assign basis to child 2
        tChildren[ 1 ]->insert_basis(  0, tBasis[   1 ] );
        tChildren[ 1 ]->insert_basis(  1, tBasis[   4 ] );
        tChildren[ 1 ]->insert_basis(  2, tBasis[  19 ] );
        tChildren[ 1 ]->insert_basis(  3, tBasis[  16 ] );
        tChildren[ 1 ]->insert_basis(  4, tBasis[  76 ] );
        tChildren[ 1 ]->insert_basis(  5, tBasis[  79 ] );
        tChildren[ 1 ]->insert_basis(  6, tBasis[  94 ] );
        tChildren[ 1 ]->insert_basis(  7, tBasis[  91 ] );
        tChildren[ 1 ]->insert_basis(  8, tBasis[   2 ] );
        tChildren[ 1 ]->insert_basis(  9, tBasis[   3 ] );
        tChildren[ 1 ]->insert_basis( 10, tBasis[   6 ] );
        tChildren[ 1 ]->insert_basis( 11, tBasis[  11 ] );
        tChildren[ 1 ]->insert_basis( 12, tBasis[  26 ] );
        tChildren[ 1 ]->insert_basis( 13, tBasis[  51 ] );
        tChildren[ 1 ]->insert_basis( 14, tBasis[   9 ] );
        tChildren[ 1 ]->insert_basis( 15, tBasis[  14 ] );
        tChildren[ 1 ]->insert_basis( 16, tBasis[  29 ] );
        tChildren[ 1 ]->insert_basis( 17, tBasis[  54 ] );
        tChildren[ 1 ]->insert_basis( 18, tBasis[  18 ] );
        tChildren[ 1 ]->insert_basis( 19, tBasis[  17 ] );
        tChildren[ 1 ]->insert_basis( 20, tBasis[  44 ] );
        tChildren[ 1 ]->insert_basis( 21, tBasis[  69 ] );
        tChildren[ 1 ]->insert_basis( 22, tBasis[  41 ] );
        tChildren[ 1 ]->insert_basis( 23, tBasis[  66 ] );
        tChildren[ 1 ]->insert_basis( 24, tBasis[  77 ] );
        tChildren[ 1 ]->insert_basis( 25, tBasis[  78 ] );
        tChildren[ 1 ]->insert_basis( 26, tBasis[  81 ] );
        tChildren[ 1 ]->insert_basis( 27, tBasis[  86 ] );
        tChildren[ 1 ]->insert_basis( 28, tBasis[  84 ] );
        tChildren[ 1 ]->insert_basis( 29, tBasis[  89 ] );
        tChildren[ 1 ]->insert_basis( 30, tBasis[  93 ] );
        tChildren[ 1 ]->insert_basis( 31, tBasis[  92 ] );
        tChildren[ 1 ]->insert_basis( 32, tBasis[   7 ] );
        tChildren[ 1 ]->insert_basis( 33, tBasis[  12 ] );
        tChildren[ 1 ]->insert_basis( 34, tBasis[  13 ] );
        tChildren[ 1 ]->insert_basis( 35, tBasis[   8 ] );
        tChildren[ 1 ]->insert_basis( 36, tBasis[  27 ] );
        tChildren[ 1 ]->insert_basis( 37, tBasis[  28 ] );
        tChildren[ 1 ]->insert_basis( 38, tBasis[  53 ] );
        tChildren[ 1 ]->insert_basis( 39, tBasis[  52 ] );
        tChildren[ 1 ]->insert_basis( 40, tBasis[  31 ] );
        tChildren[ 1 ]->insert_basis( 41, tBasis[  56 ] );
        tChildren[ 1 ]->insert_basis( 42, tBasis[  61 ] );
        tChildren[ 1 ]->insert_basis( 43, tBasis[  36 ] );
        tChildren[ 1 ]->insert_basis( 44, tBasis[  34 ] );
        tChildren[ 1 ]->insert_basis( 45, tBasis[  39 ] );
        tChildren[ 1 ]->insert_basis( 46, tBasis[  64 ] );
        tChildren[ 1 ]->insert_basis( 47, tBasis[  59 ] );
        tChildren[ 1 ]->insert_basis( 48, tBasis[  43 ] );
        tChildren[ 1 ]->insert_basis( 49, tBasis[  42 ] );
        tChildren[ 1 ]->insert_basis( 50, tBasis[  67 ] );
        tChildren[ 1 ]->insert_basis( 51, tBasis[  68 ] );
        tChildren[ 1 ]->insert_basis( 52, tBasis[  82 ] );
        tChildren[ 1 ]->insert_basis( 53, tBasis[  83 ] );
        tChildren[ 1 ]->insert_basis( 54, tBasis[  88 ] );
        tChildren[ 1 ]->insert_basis( 55, tBasis[  87 ] );
        tChildren[ 1 ]->insert_basis( 56, tBasis[  32 ] );
        tChildren[ 1 ]->insert_basis( 57, tBasis[  33 ] );
        tChildren[ 1 ]->insert_basis( 58, tBasis[  38 ] );
        tChildren[ 1 ]->insert_basis( 59, tBasis[  37 ] );
        tChildren[ 1 ]->insert_basis( 60, tBasis[  57 ] );
        tChildren[ 1 ]->insert_basis( 61, tBasis[  58 ] );
        tChildren[ 1 ]->insert_basis( 62, tBasis[  63 ] );
        tChildren[ 1 ]->insert_basis( 63, tBasis[  62 ] );

        // assign basis to child 3
        tChildren[ 2 ]->insert_basis(  0, tBasis[   5 ] );
        tChildren[ 2 ]->insert_basis(  1, tBasis[   8 ] );
        tChildren[ 2 ]->insert_basis(  2, tBasis[  23 ] );
        tChildren[ 2 ]->insert_basis(  3, tBasis[  20 ] );
        tChildren[ 2 ]->insert_basis(  4, tBasis[  80 ] );
        tChildren[ 2 ]->insert_basis(  5, tBasis[  83 ] );
        tChildren[ 2 ]->insert_basis(  6, tBasis[  98 ] );
        tChildren[ 2 ]->insert_basis(  7, tBasis[  95 ] );
        tChildren[ 2 ]->insert_basis(  8, tBasis[   6 ] );
        tChildren[ 2 ]->insert_basis(  9, tBasis[   7 ] );
        tChildren[ 2 ]->insert_basis( 10, tBasis[  10 ] );
        tChildren[ 2 ]->insert_basis( 11, tBasis[  15 ] );
        tChildren[ 2 ]->insert_basis( 12, tBasis[  30 ] );
        tChildren[ 2 ]->insert_basis( 13, tBasis[  55 ] );
        tChildren[ 2 ]->insert_basis( 14, tBasis[  13 ] );
        tChildren[ 2 ]->insert_basis( 15, tBasis[  18 ] );
        tChildren[ 2 ]->insert_basis( 16, tBasis[  33 ] );
        tChildren[ 2 ]->insert_basis( 17, tBasis[  58 ] );
        tChildren[ 2 ]->insert_basis( 18, tBasis[  22 ] );
        tChildren[ 2 ]->insert_basis( 19, tBasis[  21 ] );
        tChildren[ 2 ]->insert_basis( 20, tBasis[  48 ] );
        tChildren[ 2 ]->insert_basis( 21, tBasis[  73 ] );
        tChildren[ 2 ]->insert_basis( 22, tBasis[  45 ] );
        tChildren[ 2 ]->insert_basis( 23, tBasis[  70 ] );
        tChildren[ 2 ]->insert_basis( 24, tBasis[  81 ] );
        tChildren[ 2 ]->insert_basis( 25, tBasis[  82 ] );
        tChildren[ 2 ]->insert_basis( 26, tBasis[  85 ] );
        tChildren[ 2 ]->insert_basis( 27, tBasis[  90 ] );
        tChildren[ 2 ]->insert_basis( 28, tBasis[  88 ] );
        tChildren[ 2 ]->insert_basis( 29, tBasis[  93 ] );
        tChildren[ 2 ]->insert_basis( 30, tBasis[  97 ] );
        tChildren[ 2 ]->insert_basis( 31, tBasis[  96 ] );
        tChildren[ 2 ]->insert_basis( 32, tBasis[  11 ] );
        tChildren[ 2 ]->insert_basis( 33, tBasis[  16 ] );
        tChildren[ 2 ]->insert_basis( 34, tBasis[  17 ] );
        tChildren[ 2 ]->insert_basis( 35, tBasis[  12 ] );
        tChildren[ 2 ]->insert_basis( 36, tBasis[  31 ] );
        tChildren[ 2 ]->insert_basis( 37, tBasis[  32 ] );
        tChildren[ 2 ]->insert_basis( 38, tBasis[  57 ] );
        tChildren[ 2 ]->insert_basis( 39, tBasis[  56 ] );
        tChildren[ 2 ]->insert_basis( 40, tBasis[  35 ] );
        tChildren[ 2 ]->insert_basis( 41, tBasis[  60 ] );
        tChildren[ 2 ]->insert_basis( 42, tBasis[  65 ] );
        tChildren[ 2 ]->insert_basis( 43, tBasis[  40 ] );
        tChildren[ 2 ]->insert_basis( 44, tBasis[  38 ] );
        tChildren[ 2 ]->insert_basis( 45, tBasis[  43 ] );
        tChildren[ 2 ]->insert_basis( 46, tBasis[  68 ] );
        tChildren[ 2 ]->insert_basis( 47, tBasis[  63 ] );
        tChildren[ 2 ]->insert_basis( 48, tBasis[  47 ] );
        tChildren[ 2 ]->insert_basis( 49, tBasis[  46 ] );
        tChildren[ 2 ]->insert_basis( 50, tBasis[  71 ] );
        tChildren[ 2 ]->insert_basis( 51, tBasis[  72 ] );
        tChildren[ 2 ]->insert_basis( 52, tBasis[  86 ] );
        tChildren[ 2 ]->insert_basis( 53, tBasis[  87 ] );
        tChildren[ 2 ]->insert_basis( 54, tBasis[  92 ] );
        tChildren[ 2 ]->insert_basis( 55, tBasis[  91 ] );
        tChildren[ 2 ]->insert_basis( 56, tBasis[  36 ] );
        tChildren[ 2 ]->insert_basis( 57, tBasis[  37 ] );
        tChildren[ 2 ]->insert_basis( 58, tBasis[  42 ] );
        tChildren[ 2 ]->insert_basis( 59, tBasis[  41 ] );
        tChildren[ 2 ]->insert_basis( 60, tBasis[  61 ] );
        tChildren[ 2 ]->insert_basis( 61, tBasis[  62 ] );
        tChildren[ 2 ]->insert_basis( 62, tBasis[  67 ] );
        tChildren[ 2 ]->insert_basis( 63, tBasis[  66 ] );

        // assign basis to child 4
        tChildren[ 3 ]->insert_basis(  0, tBasis[   6 ] );
        tChildren[ 3 ]->insert_basis(  1, tBasis[   9 ] );
        tChildren[ 3 ]->insert_basis(  2, tBasis[  24 ] );
        tChildren[ 3 ]->insert_basis(  3, tBasis[  21 ] );
        tChildren[ 3 ]->insert_basis(  4, tBasis[  81 ] );
        tChildren[ 3 ]->insert_basis(  5, tBasis[  84 ] );
        tChildren[ 3 ]->insert_basis(  6, tBasis[  99 ] );
        tChildren[ 3 ]->insert_basis(  7, tBasis[  96 ] );
        tChildren[ 3 ]->insert_basis(  8, tBasis[   7 ] );
        tChildren[ 3 ]->insert_basis(  9, tBasis[   8 ] );
        tChildren[ 3 ]->insert_basis( 10, tBasis[  11 ] );
        tChildren[ 3 ]->insert_basis( 11, tBasis[  16 ] );
        tChildren[ 3 ]->insert_basis( 12, tBasis[  31 ] );
        tChildren[ 3 ]->insert_basis( 13, tBasis[  56 ] );
        tChildren[ 3 ]->insert_basis( 14, tBasis[  14 ] );
        tChildren[ 3 ]->insert_basis( 15, tBasis[  19 ] );
        tChildren[ 3 ]->insert_basis( 16, tBasis[  34 ] );
        tChildren[ 3 ]->insert_basis( 17, tBasis[  59 ] );
        tChildren[ 3 ]->insert_basis( 18, tBasis[  23 ] );
        tChildren[ 3 ]->insert_basis( 19, tBasis[  22 ] );
        tChildren[ 3 ]->insert_basis( 20, tBasis[  49 ] );
        tChildren[ 3 ]->insert_basis( 21, tBasis[  74 ] );
        tChildren[ 3 ]->insert_basis( 22, tBasis[  46 ] );
        tChildren[ 3 ]->insert_basis( 23, tBasis[  71 ] );
        tChildren[ 3 ]->insert_basis( 24, tBasis[  82 ] );
        tChildren[ 3 ]->insert_basis( 25, tBasis[  83 ] );
        tChildren[ 3 ]->insert_basis( 26, tBasis[  86 ] );
        tChildren[ 3 ]->insert_basis( 27, tBasis[  91 ] );
        tChildren[ 3 ]->insert_basis( 28, tBasis[  89 ] );
        tChildren[ 3 ]->insert_basis( 29, tBasis[  94 ] );
        tChildren[ 3 ]->insert_basis( 30, tBasis[  98 ] );
        tChildren[ 3 ]->insert_basis( 31, tBasis[  97 ] );
        tChildren[ 3 ]->insert_basis( 32, tBasis[  12 ] );
        tChildren[ 3 ]->insert_basis( 33, tBasis[  17 ] );
        tChildren[ 3 ]->insert_basis( 34, tBasis[  18 ] );
        tChildren[ 3 ]->insert_basis( 35, tBasis[  13 ] );
        tChildren[ 3 ]->insert_basis( 36, tBasis[  32 ] );
        tChildren[ 3 ]->insert_basis( 37, tBasis[  33 ] );
        tChildren[ 3 ]->insert_basis( 38, tBasis[  58 ] );
        tChildren[ 3 ]->insert_basis( 39, tBasis[  57 ] );
        tChildren[ 3 ]->insert_basis( 40, tBasis[  36 ] );
        tChildren[ 3 ]->insert_basis( 41, tBasis[  61 ] );
        tChildren[ 3 ]->insert_basis( 42, tBasis[  66 ] );
        tChildren[ 3 ]->insert_basis( 43, tBasis[  41 ] );
        tChildren[ 3 ]->insert_basis( 44, tBasis[  39 ] );
        tChildren[ 3 ]->insert_basis( 45, tBasis[  44 ] );
        tChildren[ 3 ]->insert_basis( 46, tBasis[  69 ] );
        tChildren[ 3 ]->insert_basis( 47, tBasis[  64 ] );
        tChildren[ 3 ]->insert_basis( 48, tBasis[  48 ] );
        tChildren[ 3 ]->insert_basis( 49, tBasis[  47 ] );
        tChildren[ 3 ]->insert_basis( 50, tBasis[  72 ] );
        tChildren[ 3 ]->insert_basis( 51, tBasis[  73 ] );
        tChildren[ 3 ]->insert_basis( 52, tBasis[  87 ] );
        tChildren[ 3 ]->insert_basis( 53, tBasis[  88 ] );
        tChildren[ 3 ]->insert_basis( 54, tBasis[  93 ] );
        tChildren[ 3 ]->insert_basis( 55, tBasis[  92 ] );
        tChildren[ 3 ]->insert_basis( 56, tBasis[  37 ] );
        tChildren[ 3 ]->insert_basis( 57, tBasis[  38 ] );
        tChildren[ 3 ]->insert_basis( 58, tBasis[  43 ] );
        tChildren[ 3 ]->insert_basis( 59, tBasis[  42 ] );
        tChildren[ 3 ]->insert_basis( 60, tBasis[  62 ] );
        tChildren[ 3 ]->insert_basis( 61, tBasis[  63 ] );
        tChildren[ 3 ]->insert_basis( 62, tBasis[  68 ] );
        tChildren[ 3 ]->insert_basis( 63, tBasis[  67 ] );

        // assign basis to child 5
        tChildren[ 4 ]->insert_basis(  0, tBasis[  25 ] );
        tChildren[ 4 ]->insert_basis(  1, tBasis[  28 ] );
        tChildren[ 4 ]->insert_basis(  2, tBasis[  43 ] );
        tChildren[ 4 ]->insert_basis(  3, tBasis[  40 ] );
        tChildren[ 4 ]->insert_basis(  4, tBasis[ 100 ] );
        tChildren[ 4 ]->insert_basis(  5, tBasis[ 103 ] );
        tChildren[ 4 ]->insert_basis(  6, tBasis[ 118 ] );
        tChildren[ 4 ]->insert_basis(  7, tBasis[ 115 ] );
        tChildren[ 4 ]->insert_basis(  8, tBasis[  26 ] );
        tChildren[ 4 ]->insert_basis(  9, tBasis[  27 ] );
        tChildren[ 4 ]->insert_basis( 10, tBasis[  30 ] );
        tChildren[ 4 ]->insert_basis( 11, tBasis[  35 ] );
        tChildren[ 4 ]->insert_basis( 12, tBasis[  50 ] );
        tChildren[ 4 ]->insert_basis( 13, tBasis[  75 ] );
        tChildren[ 4 ]->insert_basis( 14, tBasis[  33 ] );
        tChildren[ 4 ]->insert_basis( 15, tBasis[  38 ] );
        tChildren[ 4 ]->insert_basis( 16, tBasis[  53 ] );
        tChildren[ 4 ]->insert_basis( 17, tBasis[  78 ] );
        tChildren[ 4 ]->insert_basis( 18, tBasis[  42 ] );
        tChildren[ 4 ]->insert_basis( 19, tBasis[  41 ] );
        tChildren[ 4 ]->insert_basis( 20, tBasis[  68 ] );
        tChildren[ 4 ]->insert_basis( 21, tBasis[  93 ] );
        tChildren[ 4 ]->insert_basis( 22, tBasis[  65 ] );
        tChildren[ 4 ]->insert_basis( 23, tBasis[  90 ] );
        tChildren[ 4 ]->insert_basis( 24, tBasis[ 101 ] );
        tChildren[ 4 ]->insert_basis( 25, tBasis[ 102 ] );
        tChildren[ 4 ]->insert_basis( 26, tBasis[ 105 ] );
        tChildren[ 4 ]->insert_basis( 27, tBasis[ 110 ] );
        tChildren[ 4 ]->insert_basis( 28, tBasis[ 108 ] );
        tChildren[ 4 ]->insert_basis( 29, tBasis[ 113 ] );
        tChildren[ 4 ]->insert_basis( 30, tBasis[ 117 ] );
        tChildren[ 4 ]->insert_basis( 31, tBasis[ 116 ] );
        tChildren[ 4 ]->insert_basis( 32, tBasis[  31 ] );
        tChildren[ 4 ]->insert_basis( 33, tBasis[  36 ] );
        tChildren[ 4 ]->insert_basis( 34, tBasis[  37 ] );
        tChildren[ 4 ]->insert_basis( 35, tBasis[  32 ] );
        tChildren[ 4 ]->insert_basis( 36, tBasis[  51 ] );
        tChildren[ 4 ]->insert_basis( 37, tBasis[  52 ] );
        tChildren[ 4 ]->insert_basis( 38, tBasis[  77 ] );
        tChildren[ 4 ]->insert_basis( 39, tBasis[  76 ] );
        tChildren[ 4 ]->insert_basis( 40, tBasis[  55 ] );
        tChildren[ 4 ]->insert_basis( 41, tBasis[  80 ] );
        tChildren[ 4 ]->insert_basis( 42, tBasis[  85 ] );
        tChildren[ 4 ]->insert_basis( 43, tBasis[  60 ] );
        tChildren[ 4 ]->insert_basis( 44, tBasis[  58 ] );
        tChildren[ 4 ]->insert_basis( 45, tBasis[  63 ] );
        tChildren[ 4 ]->insert_basis( 46, tBasis[  88 ] );
        tChildren[ 4 ]->insert_basis( 47, tBasis[  83 ] );
        tChildren[ 4 ]->insert_basis( 48, tBasis[  67 ] );
        tChildren[ 4 ]->insert_basis( 49, tBasis[  66 ] );
        tChildren[ 4 ]->insert_basis( 50, tBasis[  91 ] );
        tChildren[ 4 ]->insert_basis( 51, tBasis[  92 ] );
        tChildren[ 4 ]->insert_basis( 52, tBasis[ 106 ] );
        tChildren[ 4 ]->insert_basis( 53, tBasis[ 107 ] );
        tChildren[ 4 ]->insert_basis( 54, tBasis[ 112 ] );
        tChildren[ 4 ]->insert_basis( 55, tBasis[ 111 ] );
        tChildren[ 4 ]->insert_basis( 56, tBasis[  56 ] );
        tChildren[ 4 ]->insert_basis( 57, tBasis[  57 ] );
        tChildren[ 4 ]->insert_basis( 58, tBasis[  62 ] );
        tChildren[ 4 ]->insert_basis( 59, tBasis[  61 ] );
        tChildren[ 4 ]->insert_basis( 60, tBasis[  81 ] );
        tChildren[ 4 ]->insert_basis( 61, tBasis[  82 ] );
        tChildren[ 4 ]->insert_basis( 62, tBasis[  87 ] );
        tChildren[ 4 ]->insert_basis( 63, tBasis[  86 ] );

        // assign basis to child 6
        tChildren[ 5 ]->insert_basis(  0, tBasis[  26 ] );
        tChildren[ 5 ]->insert_basis(  1, tBasis[  29 ] );
        tChildren[ 5 ]->insert_basis(  2, tBasis[  44 ] );
        tChildren[ 5 ]->insert_basis(  3, tBasis[  41 ] );
        tChildren[ 5 ]->insert_basis(  4, tBasis[ 101 ] );
        tChildren[ 5 ]->insert_basis(  5, tBasis[ 104 ] );
        tChildren[ 5 ]->insert_basis(  6, tBasis[ 119 ] );
        tChildren[ 5 ]->insert_basis(  7, tBasis[ 116 ] );
        tChildren[ 5 ]->insert_basis(  8, tBasis[  27 ] );
        tChildren[ 5 ]->insert_basis(  9, tBasis[  28 ] );
        tChildren[ 5 ]->insert_basis( 10, tBasis[  31 ] );
        tChildren[ 5 ]->insert_basis( 11, tBasis[  36 ] );
        tChildren[ 5 ]->insert_basis( 12, tBasis[  51 ] );
        tChildren[ 5 ]->insert_basis( 13, tBasis[  76 ] );
        tChildren[ 5 ]->insert_basis( 14, tBasis[  34 ] );
        tChildren[ 5 ]->insert_basis( 15, tBasis[  39 ] );
        tChildren[ 5 ]->insert_basis( 16, tBasis[  54 ] );
        tChildren[ 5 ]->insert_basis( 17, tBasis[  79 ] );
        tChildren[ 5 ]->insert_basis( 18, tBasis[  43 ] );
        tChildren[ 5 ]->insert_basis( 19, tBasis[  42 ] );
        tChildren[ 5 ]->insert_basis( 20, tBasis[  69 ] );
        tChildren[ 5 ]->insert_basis( 21, tBasis[  94 ] );
        tChildren[ 5 ]->insert_basis( 22, tBasis[  66 ] );
        tChildren[ 5 ]->insert_basis( 23, tBasis[  91 ] );
        tChildren[ 5 ]->insert_basis( 24, tBasis[ 102 ] );
        tChildren[ 5 ]->insert_basis( 25, tBasis[ 103 ] );
        tChildren[ 5 ]->insert_basis( 26, tBasis[ 106 ] );
        tChildren[ 5 ]->insert_basis( 27, tBasis[ 111 ] );
        tChildren[ 5 ]->insert_basis( 28, tBasis[ 109 ] );
        tChildren[ 5 ]->insert_basis( 29, tBasis[ 114 ] );
        tChildren[ 5 ]->insert_basis( 30, tBasis[ 118 ] );
        tChildren[ 5 ]->insert_basis( 31, tBasis[ 117 ] );
        tChildren[ 5 ]->insert_basis( 32, tBasis[  32 ] );
        tChildren[ 5 ]->insert_basis( 33, tBasis[  37 ] );
        tChildren[ 5 ]->insert_basis( 34, tBasis[  38 ] );
        tChildren[ 5 ]->insert_basis( 35, tBasis[  33 ] );
        tChildren[ 5 ]->insert_basis( 36, tBasis[  52 ] );
        tChildren[ 5 ]->insert_basis( 37, tBasis[  53 ] );
        tChildren[ 5 ]->insert_basis( 38, tBasis[  78 ] );
        tChildren[ 5 ]->insert_basis( 39, tBasis[  77 ] );
        tChildren[ 5 ]->insert_basis( 40, tBasis[  56 ] );
        tChildren[ 5 ]->insert_basis( 41, tBasis[  81 ] );
        tChildren[ 5 ]->insert_basis( 42, tBasis[  86 ] );
        tChildren[ 5 ]->insert_basis( 43, tBasis[  61 ] );
        tChildren[ 5 ]->insert_basis( 44, tBasis[  59 ] );
        tChildren[ 5 ]->insert_basis( 45, tBasis[  64 ] );
        tChildren[ 5 ]->insert_basis( 46, tBasis[  89 ] );
        tChildren[ 5 ]->insert_basis( 47, tBasis[  84 ] );
        tChildren[ 5 ]->insert_basis( 48, tBasis[  68 ] );
        tChildren[ 5 ]->insert_basis( 49, tBasis[  67 ] );
        tChildren[ 5 ]->insert_basis( 50, tBasis[  92 ] );
        tChildren[ 5 ]->insert_basis( 51, tBasis[  93 ] );
        tChildren[ 5 ]->insert_basis( 52, tBasis[ 107 ] );
        tChildren[ 5 ]->insert_basis( 53, tBasis[ 108 ] );
        tChildren[ 5 ]->insert_basis( 54, tBasis[ 113 ] );
        tChildren[ 5 ]->insert_basis( 55, tBasis[ 112 ] );
        tChildren[ 5 ]->insert_basis( 56, tBasis[  57 ] );
        tChildren[ 5 ]->insert_basis( 57, tBasis[  58 ] );
        tChildren[ 5 ]->insert_basis( 58, tBasis[  63 ] );
        tChildren[ 5 ]->insert_basis( 59, tBasis[  62 ] );
        tChildren[ 5 ]->insert_basis( 60, tBasis[  82 ] );
        tChildren[ 5 ]->insert_basis( 61, tBasis[  83 ] );
        tChildren[ 5 ]->insert_basis( 62, tBasis[  88 ] );
        tChildren[ 5 ]->insert_basis( 63, tBasis[  87 ] );

        // assign basis to child 7
        tChildren[ 6 ]->insert_basis(  0, tBasis[  30 ] );
        tChildren[ 6 ]->insert_basis(  1, tBasis[  33 ] );
        tChildren[ 6 ]->insert_basis(  2, tBasis[  48 ] );
        tChildren[ 6 ]->insert_basis(  3, tBasis[  45 ] );
        tChildren[ 6 ]->insert_basis(  4, tBasis[ 105 ] );
        tChildren[ 6 ]->insert_basis(  5, tBasis[ 108 ] );
        tChildren[ 6 ]->insert_basis(  6, tBasis[ 123 ] );
        tChildren[ 6 ]->insert_basis(  7, tBasis[ 120 ] );
        tChildren[ 6 ]->insert_basis(  8, tBasis[  31 ] );
        tChildren[ 6 ]->insert_basis(  9, tBasis[  32 ] );
        tChildren[ 6 ]->insert_basis( 10, tBasis[  35 ] );
        tChildren[ 6 ]->insert_basis( 11, tBasis[  40 ] );
        tChildren[ 6 ]->insert_basis( 12, tBasis[  55 ] );
        tChildren[ 6 ]->insert_basis( 13, tBasis[  80 ] );
        tChildren[ 6 ]->insert_basis( 14, tBasis[  38 ] );
        tChildren[ 6 ]->insert_basis( 15, tBasis[  43 ] );
        tChildren[ 6 ]->insert_basis( 16, tBasis[  58 ] );
        tChildren[ 6 ]->insert_basis( 17, tBasis[  83 ] );
        tChildren[ 6 ]->insert_basis( 18, tBasis[  47 ] );
        tChildren[ 6 ]->insert_basis( 19, tBasis[  46 ] );
        tChildren[ 6 ]->insert_basis( 20, tBasis[  73 ] );
        tChildren[ 6 ]->insert_basis( 21, tBasis[  98 ] );
        tChildren[ 6 ]->insert_basis( 22, tBasis[  70 ] );
        tChildren[ 6 ]->insert_basis( 23, tBasis[  95 ] );
        tChildren[ 6 ]->insert_basis( 24, tBasis[ 106 ] );
        tChildren[ 6 ]->insert_basis( 25, tBasis[ 107 ] );
        tChildren[ 6 ]->insert_basis( 26, tBasis[ 110 ] );
        tChildren[ 6 ]->insert_basis( 27, tBasis[ 115 ] );
        tChildren[ 6 ]->insert_basis( 28, tBasis[ 113 ] );
        tChildren[ 6 ]->insert_basis( 29, tBasis[ 118 ] );
        tChildren[ 6 ]->insert_basis( 30, tBasis[ 122 ] );
        tChildren[ 6 ]->insert_basis( 31, tBasis[ 121 ] );
        tChildren[ 6 ]->insert_basis( 32, tBasis[  36 ] );
        tChildren[ 6 ]->insert_basis( 33, tBasis[  41 ] );
        tChildren[ 6 ]->insert_basis( 34, tBasis[  42 ] );
        tChildren[ 6 ]->insert_basis( 35, tBasis[  37 ] );
        tChildren[ 6 ]->insert_basis( 36, tBasis[  56 ] );
        tChildren[ 6 ]->insert_basis( 37, tBasis[  57 ] );
        tChildren[ 6 ]->insert_basis( 38, tBasis[  82 ] );
        tChildren[ 6 ]->insert_basis( 39, tBasis[  81 ] );
        tChildren[ 6 ]->insert_basis( 40, tBasis[  60 ] );
        tChildren[ 6 ]->insert_basis( 41, tBasis[  85 ] );
        tChildren[ 6 ]->insert_basis( 42, tBasis[  90 ] );
        tChildren[ 6 ]->insert_basis( 43, tBasis[  65 ] );
        tChildren[ 6 ]->insert_basis( 44, tBasis[  63 ] );
        tChildren[ 6 ]->insert_basis( 45, tBasis[  68 ] );
        tChildren[ 6 ]->insert_basis( 46, tBasis[  93 ] );
        tChildren[ 6 ]->insert_basis( 47, tBasis[  88 ] );
        tChildren[ 6 ]->insert_basis( 48, tBasis[  72 ] );
        tChildren[ 6 ]->insert_basis( 49, tBasis[  71 ] );
        tChildren[ 6 ]->insert_basis( 50, tBasis[  96 ] );
        tChildren[ 6 ]->insert_basis( 51, tBasis[  97 ] );
        tChildren[ 6 ]->insert_basis( 52, tBasis[ 111 ] );
        tChildren[ 6 ]->insert_basis( 53, tBasis[ 112 ] );
        tChildren[ 6 ]->insert_basis( 54, tBasis[ 117 ] );
        tChildren[ 6 ]->insert_basis( 55, tBasis[ 116 ] );
        tChildren[ 6 ]->insert_basis( 56, tBasis[  61 ] );
        tChildren[ 6 ]->insert_basis( 57, tBasis[  62 ] );
        tChildren[ 6 ]->insert_basis( 58, tBasis[  67 ] );
        tChildren[ 6 ]->insert_basis( 59, tBasis[  66 ] );
        tChildren[ 6 ]->insert_basis( 60, tBasis[  86 ] );
        tChildren[ 6 ]->insert_basis( 61, tBasis[  87 ] );
        tChildren[ 6 ]->insert_basis( 62, tBasis[  92 ] );
        tChildren[ 6 ]->insert_basis( 63, tBasis[  91 ] );

        // assign basis to child 8
        tChildren[ 7 ]->insert_basis(  0, tBasis[  31 ] );
        tChildren[ 7 ]->insert_basis(  1, tBasis[  34 ] );
        tChildren[ 7 ]->insert_basis(  2, tBasis[  49 ] );
        tChildren[ 7 ]->insert_basis(  3, tBasis[  46 ] );
        tChildren[ 7 ]->insert_basis(  4, tBasis[ 106 ] );
        tChildren[ 7 ]->insert_basis(  5, tBasis[ 109 ] );
        tChildren[ 7 ]->insert_basis(  6, tBasis[ 124 ] );
        tChildren[ 7 ]->insert_basis(  7, tBasis[ 121 ] );
        tChildren[ 7 ]->insert_basis(  8, tBasis[  32 ] );
        tChildren[ 7 ]->insert_basis(  9, tBasis[  33 ] );
        tChildren[ 7 ]->insert_basis( 10, tBasis[  36 ] );
        tChildren[ 7 ]->insert_basis( 11, tBasis[  41 ] );
        tChildren[ 7 ]->insert_basis( 12, tBasis[  56 ] );
        tChildren[ 7 ]->insert_basis( 13, tBasis[  81 ] );
        tChildren[ 7 ]->insert_basis( 14, tBasis[  39 ] );
        tChildren[ 7 ]->insert_basis( 15, tBasis[  44 ] );
        tChildren[ 7 ]->insert_basis( 16, tBasis[  59 ] );
        tChildren[ 7 ]->insert_basis( 17, tBasis[  84 ] );
        tChildren[ 7 ]->insert_basis( 18, tBasis[  48 ] );
        tChildren[ 7 ]->insert_basis( 19, tBasis[  47 ] );
        tChildren[ 7 ]->insert_basis( 20, tBasis[  74 ] );
        tChildren[ 7 ]->insert_basis( 21, tBasis[  99 ] );
        tChildren[ 7 ]->insert_basis( 22, tBasis[  71 ] );
        tChildren[ 7 ]->insert_basis( 23, tBasis[  96 ] );
        tChildren[ 7 ]->insert_basis( 24, tBasis[ 107 ] );
        tChildren[ 7 ]->insert_basis( 25, tBasis[ 108 ] );
        tChildren[ 7 ]->insert_basis( 26, tBasis[ 111 ] );
        tChildren[ 7 ]->insert_basis( 27, tBasis[ 116 ] );
        tChildren[ 7 ]->insert_basis( 28, tBasis[ 114 ] );
        tChildren[ 7 ]->insert_basis( 29, tBasis[ 119 ] );
        tChildren[ 7 ]->insert_basis( 30, tBasis[ 123 ] );
        tChildren[ 7 ]->insert_basis( 31, tBasis[ 122 ] );
        tChildren[ 7 ]->insert_basis( 32, tBasis[  37 ] );
        tChildren[ 7 ]->insert_basis( 33, tBasis[  42 ] );
        tChildren[ 7 ]->insert_basis( 34, tBasis[  43 ] );
        tChildren[ 7 ]->insert_basis( 35, tBasis[  38 ] );
        tChildren[ 7 ]->insert_basis( 36, tBasis[  57 ] );
        tChildren[ 7 ]->insert_basis( 37, tBasis[  58 ] );
        tChildren[ 7 ]->insert_basis( 38, tBasis[  83 ] );
        tChildren[ 7 ]->insert_basis( 39, tBasis[  82 ] );
        tChildren[ 7 ]->insert_basis( 40, tBasis[  61 ] );
        tChildren[ 7 ]->insert_basis( 41, tBasis[  86 ] );
        tChildren[ 7 ]->insert_basis( 42, tBasis[  91 ] );
        tChildren[ 7 ]->insert_basis( 43, tBasis[  66 ] );
        tChildren[ 7 ]->insert_basis( 44, tBasis[  64 ] );
        tChildren[ 7 ]->insert_basis( 45, tBasis[  69 ] );
        tChildren[ 7 ]->insert_basis( 46, tBasis[  94 ] );
        tChildren[ 7 ]->insert_basis( 47, tBasis[  89 ] );
        tChildren[ 7 ]->insert_basis( 48, tBasis[  73 ] );
        tChildren[ 7 ]->insert_basis( 49, tBasis[  72 ] );
        tChildren[ 7 ]->insert_basis( 50, tBasis[  97 ] );
        tChildren[ 7 ]->insert_basis( 51, tBasis[  98 ] );
        tChildren[ 7 ]->insert_basis( 52, tBasis[ 112 ] );
        tChildren[ 7 ]->insert_basis( 53, tBasis[ 113 ] );
        tChildren[ 7 ]->insert_basis( 54, tBasis[ 118 ] );
        tChildren[ 7 ]->insert_basis( 55, tBasis[ 117 ] );
        tChildren[ 7 ]->insert_basis( 56, tBasis[  62 ] );
        tChildren[ 7 ]->insert_basis( 57, tBasis[  63 ] );
        tChildren[ 7 ]->insert_basis( 58, tBasis[  68 ] );
        tChildren[ 7 ]->insert_basis( 59, tBasis[  67 ] );
        tChildren[ 7 ]->insert_basis( 60, tBasis[  87 ] );
        tChildren[ 7 ]->insert_basis( 61, tBasis[  88 ] );
        tChildren[ 7 ]->insert_basis( 62, tBasis[  93 ] );
        tChildren[ 7 ]->insert_basis( 63, tBasis[  92 ] );

        // set basis flag of element
        mChildrenBasisFlag = true;
        
        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_HEX64_HPP_ */

