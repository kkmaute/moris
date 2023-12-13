/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_BSpline_Element_Quad9.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD9_HPP_
#define SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD9_HPP_

#include "cl_HMR_BSpline_Element.hpp"
#include "fn_HMR_get_basis_neighbors_2d.hpp"

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
    BSpline_Element< 2, 2, 0 >::get_geometry_type() const
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
    BSpline_Element< 2, 2, 0 >::get_ijk_of_basis_function(
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
                break;
            }
            case( 1 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case( 2 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case( 3 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case( 4 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 0 ;
                break;
            }
            case( 5 ) :
            {
                aIJK[ 0 ] = 2 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case( 6 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 2 ;
                break;
            }
            case( 7 ) :
            {
                aIJK[ 0 ] = 0 ;
                aIJK[ 1 ] = 1 ;
                break;
            }
            case( 8 ) :
            {
                aIJK[ 0 ] = 1 ;
                aIJK[ 1 ] = 1 ;
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
    BSpline_Element< 2, 2, 0 >::link_basis_functions_with_neighbors(
          moris::Cell< Element* > & aAllElementsOnProc )
    {
         // initialize frame of basis around basis from this element
         Basis_Function* tBasis[ 16 ] = { nullptr };

         // get pointer to neighbor  0
         Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

         // test if neighbor  0 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   1 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   2 ] = tNeighbor->get_basis_function(   4 );
             tBasis[   3 ] = tNeighbor->get_basis_function(   1 );
         }

         // get pointer to neighbor  1
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

         // test if neighbor  1 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   6 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   8 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  10 ] = tNeighbor->get_basis_function(   2 );
         }

         // get pointer to neighbor  2
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

         // test if neighbor  2 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[  12 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  13 ] = tNeighbor->get_basis_function(   6 );
             tBasis[  14 ] = tNeighbor->get_basis_function(   2 );
         }

         // get pointer to neighbor  3
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

         // test if neighbor  3 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   5 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   7 ] = tNeighbor->get_basis_function(   7 );
             tBasis[   9 ] = tNeighbor->get_basis_function(   3 );
         }

         // get pointer to neighbor  4
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

         // test if neighbor  4 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   0 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   1 ] = tNeighbor->get_basis_function(   4 );
             tBasis[   2 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   5 ] = tNeighbor->get_basis_function(   7 );
             tBasis[   7 ] = tNeighbor->get_basis_function(   3 );
         }

         // get pointer to neighbor  5
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

         // test if neighbor  5 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   2 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   3 ] = tNeighbor->get_basis_function(   4 );
             tBasis[   4 ] = tNeighbor->get_basis_function(   1 );
             tBasis[   6 ] = tNeighbor->get_basis_function(   5 );
             tBasis[   8 ] = tNeighbor->get_basis_function(   2 );
         }

         // get pointer to neighbor  6
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

         // test if neighbor  6 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   8 ] = tNeighbor->get_basis_function(   1 );
             tBasis[  10 ] = tNeighbor->get_basis_function(   5 );
             tBasis[  13 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  14 ] = tNeighbor->get_basis_function(   6 );
             tBasis[  15 ] = tNeighbor->get_basis_function(   2 );
         }

         // get pointer to neighbor  7
         tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

         // test if neighbor  7 exists
         if ( tNeighbor != nullptr )
         {
             // copy pointers into frame
             tBasis[   7 ] = tNeighbor->get_basis_function(   0 );
             tBasis[   9 ] = tNeighbor->get_basis_function(   7 );
             tBasis[  11 ] = tNeighbor->get_basis_function(   3 );
             tBasis[  12 ] = tNeighbor->get_basis_function(   6 );
             tBasis[  13 ] = tNeighbor->get_basis_function(   2 );
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
                 mBasis[   0 ]->insert_neighbor(  2, mBasis[   7 ] );
                 mBasis[   0 ]->insert_neighbor(  3, tBasis[   5 ] );
                 mBasis[   0 ]->insert_neighbor(  4, tBasis[   0 ] );
                 mBasis[   0 ]->insert_neighbor(  5, tBasis[   2 ] );
                 mBasis[   0 ]->insert_neighbor(  6, mBasis[   8 ] );
                 mBasis[   0 ]->insert_neighbor(  7, tBasis[   7 ] );

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
                 mBasis[   1 ]->insert_neighbor(  0, tBasis[   3 ] );
                 mBasis[   1 ]->insert_neighbor(  1, tBasis[   6 ] );
                 mBasis[   1 ]->insert_neighbor(  2, mBasis[   5 ] );
                 mBasis[   1 ]->insert_neighbor(  3, mBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor(  4, tBasis[   2 ] );
                 mBasis[   1 ]->insert_neighbor(  5, tBasis[   4 ] );
                 mBasis[   1 ]->insert_neighbor(  6, tBasis[   8 ] );
                 mBasis[   1 ]->insert_neighbor(  7, mBasis[   8 ] );

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
                 mBasis[   2 ]->insert_neighbor(  0, mBasis[   5 ] );
                 mBasis[   2 ]->insert_neighbor(  1, tBasis[  10 ] );
                 mBasis[   2 ]->insert_neighbor(  2, tBasis[  14 ] );
                 mBasis[   2 ]->insert_neighbor(  3, mBasis[   6 ] );
                 mBasis[   2 ]->insert_neighbor(  4, mBasis[   8 ] );
                 mBasis[   2 ]->insert_neighbor(  5, tBasis[   8 ] );
                 mBasis[   2 ]->insert_neighbor(  6, tBasis[  15 ] );
                 mBasis[   2 ]->insert_neighbor(  7, tBasis[  13 ] );

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
                 mBasis[   3 ]->insert_neighbor(  0, mBasis[   7 ] );
                 mBasis[   3 ]->insert_neighbor(  1, mBasis[   6 ] );
                 mBasis[   3 ]->insert_neighbor(  2, tBasis[  12 ] );
                 mBasis[   3 ]->insert_neighbor(  3, tBasis[   9 ] );
                 mBasis[   3 ]->insert_neighbor(  4, tBasis[   7 ] );
                 mBasis[   3 ]->insert_neighbor(  5, mBasis[   8 ] );
                 mBasis[   3 ]->insert_neighbor(  6, tBasis[  13 ] );
                 mBasis[   3 ]->insert_neighbor(  7, tBasis[  11 ] );

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
                 mBasis[   4 ]->insert_neighbor(  1, mBasis[   1 ] );
                 mBasis[   4 ]->insert_neighbor(  2, mBasis[   8 ] );
                 mBasis[   4 ]->insert_neighbor(  3, mBasis[   0 ] );
                 mBasis[   4 ]->insert_neighbor(  4, tBasis[   1 ] );
                 mBasis[   4 ]->insert_neighbor(  5, tBasis[   3 ] );
                 mBasis[   4 ]->insert_neighbor(  6, mBasis[   5 ] );
                 mBasis[   4 ]->insert_neighbor(  7, mBasis[   7 ] );

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
                 mBasis[   5 ]->insert_neighbor(  0, mBasis[   1 ] );
                 mBasis[   5 ]->insert_neighbor(  1, tBasis[   8 ] );
                 mBasis[   5 ]->insert_neighbor(  2, mBasis[   2 ] );
                 mBasis[   5 ]->insert_neighbor(  3, mBasis[   8 ] );
                 mBasis[   5 ]->insert_neighbor(  4, mBasis[   4 ] );
                 mBasis[   5 ]->insert_neighbor(  5, tBasis[   6 ] );
                 mBasis[   5 ]->insert_neighbor(  6, tBasis[  10 ] );
                 mBasis[   5 ]->insert_neighbor(  7, mBasis[   6 ] );

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
                 mBasis[   6 ]->insert_neighbor(  0, mBasis[   8 ] );
                 mBasis[   6 ]->insert_neighbor(  1, mBasis[   2 ] );
                 mBasis[   6 ]->insert_neighbor(  2, tBasis[  13 ] );
                 mBasis[   6 ]->insert_neighbor(  3, mBasis[   3 ] );
                 mBasis[   6 ]->insert_neighbor(  4, mBasis[   7 ] );
                 mBasis[   6 ]->insert_neighbor(  5, mBasis[   5 ] );
                 mBasis[   6 ]->insert_neighbor(  6, tBasis[  14 ] );
                 mBasis[   6 ]->insert_neighbor(  7, tBasis[  12 ] );

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
                 mBasis[   7 ]->insert_neighbor(  0, mBasis[   0 ] );
                 mBasis[   7 ]->insert_neighbor(  1, mBasis[   8 ] );
                 mBasis[   7 ]->insert_neighbor(  2, mBasis[   3 ] );
                 mBasis[   7 ]->insert_neighbor(  3, tBasis[   7 ] );
                 mBasis[   7 ]->insert_neighbor(  4, tBasis[   5 ] );
                 mBasis[   7 ]->insert_neighbor(  5, mBasis[   4 ] );
                 mBasis[   7 ]->insert_neighbor(  6, mBasis[   6 ] );
                 mBasis[   7 ]->insert_neighbor(  7, tBasis[   9 ] );

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
                 mBasis[   8 ]->insert_neighbor(  0, mBasis[   4 ] );
                 mBasis[   8 ]->insert_neighbor(  1, mBasis[   5 ] );
                 mBasis[   8 ]->insert_neighbor(  2, mBasis[   6 ] );
                 mBasis[   8 ]->insert_neighbor(  3, mBasis[   7 ] );
                 mBasis[   8 ]->insert_neighbor(  4, mBasis[   0 ] );
                 mBasis[   8 ]->insert_neighbor(  5, mBasis[   1 ] );
                 mBasis[   8 ]->insert_neighbor(  6, mBasis[   2 ] );
                 mBasis[   8 ]->insert_neighbor(  7, mBasis[   3 ] );

                 // flag this basis
                 mBasis[   8 ]->flag();
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
    luint BSpline_Element< 2, 2, 0 >::refine_basis( uint aBasisNumber )
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
                // create temporary container for children
                tBasis->init_children_container();

                // pointer to basis neighbor
                Basis_Function* tNeighbor;

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
                        tBasis->insert_child(  0, tNeighbor->get_child( 10 ) );
                        tBasis->insert_child(  1, tNeighbor->get_child( 11 ) );
                        tBasis->insert_child(  4, tNeighbor->get_child( 14 ) );
                        tBasis->insert_child(  5, tNeighbor->get_child( 15 ) );
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
                        tBasis->insert_child(  2, tNeighbor->get_child(  8 ) );
                        tBasis->insert_child(  3, tNeighbor->get_child(  9 ) );
                        tBasis->insert_child(  6, tNeighbor->get_child( 12 ) );
                        tBasis->insert_child(  7, tNeighbor->get_child( 13 ) );
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
                        tBasis->insert_child( 10, tNeighbor->get_child(  0 ) );
                        tBasis->insert_child( 11, tNeighbor->get_child(  1 ) );
                        tBasis->insert_child( 14, tNeighbor->get_child(  4 ) );
                        tBasis->insert_child( 15, tNeighbor->get_child(  5 ) );
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
                        tBasis->insert_child(  8, tNeighbor->get_child(  2 ) );
                        tBasis->insert_child(  9, tNeighbor->get_child(  3 ) );
                        tBasis->insert_child( 12, tNeighbor->get_child(  6 ) );
                        tBasis->insert_child( 13, tNeighbor->get_child(  7 ) );
                    }
                }

                // level of child basis
                uint tLevel = tBasis->get_level() + 1;

                // position of basis
                const luint* tParentIJ  = tBasis->get_ijk();

                // minumum i-position
                luint tIMin = 2*tParentIJ[ 0 ];

                // minumum j-position
                luint tJMin = 2*tParentIJ[ 1 ];

                // maximum i-position
                luint tIMax = tIMin + 4;

                // maximum j-position
                luint tJMax = tJMin + 4;

                // initialize counter
                uint tChildIndex = 0;

                // loop over all positions
                for( luint j=tJMin; j<tJMax; ++j )
                {
                    for( luint i=tIMin; i<tIMax; ++i )
                    {
                        // test if child does not exist
                        if( tBasis->get_child( tChildIndex ) == nullptr )
                        {
                            // calculate i-j position of child
                            luint tIJ[ 2 ] = { i, j };

                            // create child
                            tBasis->insert_child( tChildIndex,
                                 new BSpline< 2, 2, 0 >( tIJ, tLevel, gNoProcOwner ) );

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
    luint BSpline_Element< 2, 2, 0 >::refine( moris::Cell< Element* > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // refine basis if they have not been refined already
        for( uint k=0; k<9; ++k )
        {
            tBasisCounter += this->refine_basis( k );
        }

        // initialize temporary basis pattern
        Basis_Function* tBasis[ 16 ] = { nullptr };

        // populate basis pattern
        if ( mBasis[   0 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   0 ]->get_child(  10 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   0 ]->get_child(  11 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   0 ]->get_child(  14 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   0 ]->get_child(  15 );
            }
        }

        if ( mBasis[   1 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   1 ]->get_child(   8 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   1 ]->get_child(   9 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   1 ]->get_child(  12 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   1 ]->get_child(  13 );
            }
        }

        if ( mBasis[   2 ] != nullptr )
        {
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   2 ]->get_child(   0 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   2 ]->get_child(   1 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   2 ]->get_child(   4 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   2 ]->get_child(   5 );
            }
        }

        if ( mBasis[   3 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   3 ]->get_child(   2 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   3 ]->get_child(   3 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   3 ]->get_child(   6 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   3 ]->get_child(   7 );
            }
        }

        if ( mBasis[   4 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   4 ]->get_child(   8 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   4 ]->get_child(   9 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   4 ]->get_child(  10 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   4 ]->get_child(  11 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   4 ]->get_child(  12 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   4 ]->get_child(  13 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   4 ]->get_child(  14 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   4 ]->get_child(  15 );
            }
        }

        if ( mBasis[   5 ] != nullptr )
        {
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   5 ]->get_child(   0 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   5 ]->get_child(   1 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   5 ]->get_child(   4 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   5 ]->get_child(   5 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   5 ]->get_child(   8 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   5 ]->get_child(   9 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   5 ]->get_child(  12 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   5 ]->get_child(  13 );
            }
        }

        if ( mBasis[   6 ] != nullptr )
        {
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   6 ]->get_child(   0 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   6 ]->get_child(   1 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   6 ]->get_child(   2 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   6 ]->get_child(   3 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   6 ]->get_child(   4 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   6 ]->get_child(   5 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   6 ]->get_child(   6 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   6 ]->get_child(   7 );
            }
        }

        if ( mBasis[   7 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   7 ]->get_child(   2 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   7 ]->get_child(   3 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   7 ]->get_child(   6 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   7 ]->get_child(   7 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   7 ]->get_child(  10 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   7 ]->get_child(  11 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   7 ]->get_child(  14 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   7 ]->get_child(  15 );
            }
        }

        if ( mBasis[   8 ] != nullptr )
        {
            if ( tBasis[   0 ] == nullptr )
            {
                tBasis[   0 ] = mBasis[   8 ]->get_child(   0 );
            }
            if ( tBasis[   1 ] == nullptr )
            {
                tBasis[   1 ] = mBasis[   8 ]->get_child(   1 );
            }
            if ( tBasis[   2 ] == nullptr )
            {
                tBasis[   2 ] = mBasis[   8 ]->get_child(   2 );
            }
            if ( tBasis[   3 ] == nullptr )
            {
                tBasis[   3 ] = mBasis[   8 ]->get_child(   3 );
            }
            if ( tBasis[   4 ] == nullptr )
            {
                tBasis[   4 ] = mBasis[   8 ]->get_child(   4 );
            }
            if ( tBasis[   5 ] == nullptr )
            {
                tBasis[   5 ] = mBasis[   8 ]->get_child(   5 );
            }
            if ( tBasis[   6 ] == nullptr )
            {
                tBasis[   6 ] = mBasis[   8 ]->get_child(   6 );
            }
            if ( tBasis[   7 ] == nullptr )
            {
                tBasis[   7 ] = mBasis[   8 ]->get_child(   7 );
            }
            if ( tBasis[   8 ] == nullptr )
            {
                tBasis[   8 ] = mBasis[   8 ]->get_child(   8 );
            }
            if ( tBasis[   9 ] == nullptr )
            {
                tBasis[   9 ] = mBasis[   8 ]->get_child(   9 );
            }
            if ( tBasis[  10 ] == nullptr )
            {
                tBasis[  10 ] = mBasis[   8 ]->get_child(  10 );
            }
            if ( tBasis[  11 ] == nullptr )
            {
                tBasis[  11 ] = mBasis[   8 ]->get_child(  11 );
            }
            if ( tBasis[  12 ] == nullptr )
            {
                tBasis[  12 ] = mBasis[   8 ]->get_child(  12 );
            }
            if ( tBasis[  13 ] == nullptr )
            {
                tBasis[  13 ] = mBasis[   8 ]->get_child(  13 );
            }
            if ( tBasis[  14 ] == nullptr )
            {
                tBasis[  14 ] = mBasis[   8 ]->get_child(  14 );
            }
            if ( tBasis[  15 ] == nullptr )
            {
                tBasis[  15 ] = mBasis[   8 ]->get_child(  15 );
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
        tChildren[ 0 ]->insert_basis_function(  0, tBasis[   0 ] );
        tChildren[ 0 ]->insert_basis_function(  1, tBasis[   2 ] );
        tChildren[ 0 ]->insert_basis_function(  2, tBasis[  10 ] );
        tChildren[ 0 ]->insert_basis_function(  3, tBasis[   8 ] );
        tChildren[ 0 ]->insert_basis_function(  4, tBasis[   1 ] );
        tChildren[ 0 ]->insert_basis_function(  5, tBasis[   6 ] );
        tChildren[ 0 ]->insert_basis_function(  6, tBasis[   9 ] );
        tChildren[ 0 ]->insert_basis_function(  7, tBasis[   4 ] );
        tChildren[ 0 ]->insert_basis_function(  8, tBasis[   5 ] );

        // assign basis to child 2
        tChildren[ 1 ]->insert_basis_function(  0, tBasis[   1 ] );
        tChildren[ 1 ]->insert_basis_function(  1, tBasis[   3 ] );
        tChildren[ 1 ]->insert_basis_function(  2, tBasis[  11 ] );
        tChildren[ 1 ]->insert_basis_function(  3, tBasis[   9 ] );
        tChildren[ 1 ]->insert_basis_function(  4, tBasis[   2 ] );
        tChildren[ 1 ]->insert_basis_function(  5, tBasis[   7 ] );
        tChildren[ 1 ]->insert_basis_function(  6, tBasis[  10 ] );
        tChildren[ 1 ]->insert_basis_function(  7, tBasis[   5 ] );
        tChildren[ 1 ]->insert_basis_function(  8, tBasis[   6 ] );

        // assign basis to child 3
        tChildren[ 2 ]->insert_basis_function(  0, tBasis[   4 ] );
        tChildren[ 2 ]->insert_basis_function(  1, tBasis[   6 ] );
        tChildren[ 2 ]->insert_basis_function(  2, tBasis[  14 ] );
        tChildren[ 2 ]->insert_basis_function(  3, tBasis[  12 ] );
        tChildren[ 2 ]->insert_basis_function(  4, tBasis[   5 ] );
        tChildren[ 2 ]->insert_basis_function(  5, tBasis[  10 ] );
        tChildren[ 2 ]->insert_basis_function(  6, tBasis[  13 ] );
        tChildren[ 2 ]->insert_basis_function(  7, tBasis[   8 ] );
        tChildren[ 2 ]->insert_basis_function(  8, tBasis[   9 ] );

        // assign basis to child 4
        tChildren[ 3 ]->insert_basis_function(  0, tBasis[   5 ] );
        tChildren[ 3 ]->insert_basis_function(  1, tBasis[   7 ] );
        tChildren[ 3 ]->insert_basis_function(  2, tBasis[  15 ] );
        tChildren[ 3 ]->insert_basis_function(  3, tBasis[  13 ] );
        tChildren[ 3 ]->insert_basis_function(  4, tBasis[   6 ] );
        tChildren[ 3 ]->insert_basis_function(  5, tBasis[  11 ] );
        tChildren[ 3 ]->insert_basis_function(  6, tBasis[  14 ] );
        tChildren[ 3 ]->insert_basis_function(  7, tBasis[   9 ] );
        tChildren[ 3 ]->insert_basis_function(  8, tBasis[  10 ] );

        // set basis flag of element
        mChildrenBasisFlag = true;
        
        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_BSPLINE_ELEMENT_QUAD9_HPP_ */

