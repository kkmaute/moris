/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Hex64.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX64_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX64_HPP_

#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_MTK_Cell_Info_Hex64.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------
    template<>
    inline
    void
    Lagrange_Element< 3, 64 >::set_cell_info()
    {
        std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Hex64 >();

        this->set_mtk_cell_info( tCellInfo );
    }

// ----------------------------------------------------------------------------

    /**
     * string needed for gmsh output
     *
     * @return std::string
     *
     */
    template<>
    inline
    std::string
    Lagrange_Element< 3, 64 >::get_gmsh_string()
    {
        // gmsh type - number of tags - physical tag - geometry tag
        std::string aString = "92 2 0 1";

        // loop over all nodes
        for( uint k=0; k<64; ++k )
        {
            // add node index to string
            aString += " " + std::to_string(
                this->get_basis( k )->get_memory_index() + 1 );
        }

        // return the string that goes into the gmsh file
        return aString;
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
    inline
    void
    Lagrange_Element< 3, 64 >::get_ijk_of_basis(
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
        aIJK[ 0 ] += 3*tElIJK[ 0 ];
        aIJK[ 1 ] += 3*tElIJK[ 1 ];
        aIJK[ 2 ] += 3*tElIJK[ 2 ];
    }

    // ----------------------------------------------------------------------------

    /**
     * Creates all bases on the coarsest level.
     *
     * @param aAllElementsOnProc Cell containing all elements including the aura
     * @return Number of created bases
     */
    template<>
    inline
    luint Lagrange_Element< 3, 64 >::create_basis_on_level_zero(
          moris::Cell< Element * > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;

        // initialize container for nodes
        this->init_basis_container();

        // get pointer to neighbor 4
        Element* tNeighbor
            = this->get_neighbor( aAllElementsOnProc, 4 );

        // test if neighbor 4 exists
        if ( tNeighbor != nullptr )
        {
            // copy nodes from this neighbor
            mNodes[  0 ] = tNeighbor->get_basis(  4 );
            mNodes[  1 ] = tNeighbor->get_basis(  5 );
            mNodes[  2 ] = tNeighbor->get_basis(  6 );
            mNodes[  3 ] = tNeighbor->get_basis(  7 );
            mNodes[  8 ] = tNeighbor->get_basis( 24 );
            mNodes[  9 ] = tNeighbor->get_basis( 25 );
            mNodes[ 10 ] = tNeighbor->get_basis( 26 );
            mNodes[ 11 ] = tNeighbor->get_basis( 27 );
            mNodes[ 14 ] = tNeighbor->get_basis( 28 );
            mNodes[ 15 ] = tNeighbor->get_basis( 29 );
            mNodes[ 18 ] = tNeighbor->get_basis( 30 );
            mNodes[ 19 ] = tNeighbor->get_basis( 31 );
            mNodes[ 32 ] = tNeighbor->get_basis( 52 );
            mNodes[ 33 ] = tNeighbor->get_basis( 55 );
            mNodes[ 34 ] = tNeighbor->get_basis( 54 );
            mNodes[ 35 ] = tNeighbor->get_basis( 53 );
        }

        // get pointer to neighbor 0
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

        // test if neighbor 0 exists
        if ( tNeighbor != nullptr )
        {
            // copy nodes from this neighbor
            mNodes[  0 ] = tNeighbor->get_basis(  3 );
            mNodes[  1 ] = tNeighbor->get_basis(  2 );
            mNodes[  4 ] = tNeighbor->get_basis(  7 );
            mNodes[  5 ] = tNeighbor->get_basis(  6 );
            mNodes[  8 ] = tNeighbor->get_basis( 19 );
            mNodes[  9 ] = tNeighbor->get_basis( 18 );
            mNodes[ 12 ] = tNeighbor->get_basis( 22 );
            mNodes[ 13 ] = tNeighbor->get_basis( 23 );
            mNodes[ 16 ] = tNeighbor->get_basis( 20 );
            mNodes[ 17 ] = tNeighbor->get_basis( 21 );
            mNodes[ 24 ] = tNeighbor->get_basis( 31 );
            mNodes[ 25 ] = tNeighbor->get_basis( 30 );
            mNodes[ 36 ] = tNeighbor->get_basis( 49 );
            mNodes[ 37 ] = tNeighbor->get_basis( 48 );
            mNodes[ 38 ] = tNeighbor->get_basis( 51 );
            mNodes[ 39 ] = tNeighbor->get_basis( 50 );
        }

        // get pointer to neighbor 3
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

        // test if neighbor 3 exists
        if ( tNeighbor != nullptr )
        {
            // copy nodes from this neighbor
            mNodes[  0 ] = tNeighbor->get_basis(  1 );
            mNodes[  3 ] = tNeighbor->get_basis(  2 );
            mNodes[  4 ] = tNeighbor->get_basis(  5 );
            mNodes[  7 ] = tNeighbor->get_basis(  6 );
            mNodes[ 10 ] = tNeighbor->get_basis( 14 );
            mNodes[ 11 ] = tNeighbor->get_basis( 15 );
            mNodes[ 12 ] = tNeighbor->get_basis( 16 );
            mNodes[ 13 ] = tNeighbor->get_basis( 17 );
            mNodes[ 22 ] = tNeighbor->get_basis( 20 );
            mNodes[ 23 ] = tNeighbor->get_basis( 21 );
            mNodes[ 26 ] = tNeighbor->get_basis( 28 );
            mNodes[ 27 ] = tNeighbor->get_basis( 29 );
            mNodes[ 40 ] = tNeighbor->get_basis( 44 );
            mNodes[ 41 ] = tNeighbor->get_basis( 47 );
            mNodes[ 42 ] = tNeighbor->get_basis( 46 );
            mNodes[ 43 ] = tNeighbor->get_basis( 45 );
        }

        // loop over all nodes
        for( uint k=0; k<64; ++k )
        {
            // test if node exists
            if( mNodes[ k ] == nullptr )
            {
                // create node
                this->create_basis( k );

                // increment node counter
                tBasisCounter++;
            }
        }

        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------

    /**
     * Creates bases for children of refined elements.
     *
     * @param aAllElementsOnProc Cell containing all elements including the aura
     * @return Number of created bases
     */
    template<>
    inline
    luint Lagrange_Element< 3, 64 >::create_basis_for_children(
        moris::Cell< Element * > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // create temporary array containing all nodes
        Basis* tNodes[ 343 ] = { nullptr };

        // copy my own nodes into this array
        tNodes[   0 ] = mNodes[   0 ];
        tNodes[   2 ] = mNodes[   8 ];
        tNodes[   4 ] = mNodes[   9 ];
        tNodes[   6 ] = mNodes[   1 ];
        tNodes[  14 ] = mNodes[  10 ];
        tNodes[  16 ] = mNodes[  32 ];
        tNodes[  18 ] = mNodes[  35 ];
        tNodes[  20 ] = mNodes[  14 ];
        tNodes[  28 ] = mNodes[  11 ];
        tNodes[  30 ] = mNodes[  33 ];
        tNodes[  32 ] = mNodes[  34 ];
        tNodes[  34 ] = mNodes[  15 ];
        tNodes[  42 ] = mNodes[   3 ];
        tNodes[  44 ] = mNodes[  19 ];
        tNodes[  46 ] = mNodes[  18 ];
        tNodes[  48 ] = mNodes[   2 ];
        tNodes[  98 ] = mNodes[  12 ];
        tNodes[ 100 ] = mNodes[  36 ];
        tNodes[ 102 ] = mNodes[  37 ];
        tNodes[ 104 ] = mNodes[  16 ];
        tNodes[ 112 ] = mNodes[  40 ];
        tNodes[ 114 ] = mNodes[  56 ];
        tNodes[ 116 ] = mNodes[  57 ];
        tNodes[ 118 ] = mNodes[  44 ];
        tNodes[ 126 ] = mNodes[  43 ];
        tNodes[ 128 ] = mNodes[  59 ];
        tNodes[ 130 ] = mNodes[  58 ];
        tNodes[ 132 ] = mNodes[  45 ];
        tNodes[ 140 ] = mNodes[  22 ];
        tNodes[ 142 ] = mNodes[  49 ];
        tNodes[ 144 ] = mNodes[  48 ];
        tNodes[ 146 ] = mNodes[  20 ];
        tNodes[ 196 ] = mNodes[  13 ];
        tNodes[ 198 ] = mNodes[  39 ];
        tNodes[ 200 ] = mNodes[  38 ];
        tNodes[ 202 ] = mNodes[  17 ];
        tNodes[ 210 ] = mNodes[  41 ];
        tNodes[ 212 ] = mNodes[  60 ];
        tNodes[ 214 ] = mNodes[  61 ];
        tNodes[ 216 ] = mNodes[  47 ];
        tNodes[ 224 ] = mNodes[  42 ];
        tNodes[ 226 ] = mNodes[  63 ];
        tNodes[ 228 ] = mNodes[  62 ];
        tNodes[ 230 ] = mNodes[  46 ];
        tNodes[ 238 ] = mNodes[  23 ];
        tNodes[ 240 ] = mNodes[  50 ];
        tNodes[ 242 ] = mNodes[  51 ];
        tNodes[ 244 ] = mNodes[  21 ];
        tNodes[ 294 ] = mNodes[   4 ];
        tNodes[ 296 ] = mNodes[  24 ];
        tNodes[ 298 ] = mNodes[  25 ];
        tNodes[ 300 ] = mNodes[   5 ];
        tNodes[ 308 ] = mNodes[  26 ];
        tNodes[ 310 ] = mNodes[  52 ];
        tNodes[ 312 ] = mNodes[  53 ];
        tNodes[ 314 ] = mNodes[  28 ];
        tNodes[ 322 ] = mNodes[  27 ];
        tNodes[ 324 ] = mNodes[  55 ];
        tNodes[ 326 ] = mNodes[  54 ];
        tNodes[ 328 ] = mNodes[  29 ];
        tNodes[ 336 ] = mNodes[   7 ];
        tNodes[ 338 ] = mNodes[  31 ];
        tNodes[ 340 ] = mNodes[  30 ];
        tNodes[ 342 ] = mNodes[   6 ];

        // get pointer to neighbor
        Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

        // test if neighbor 0 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 0 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 2
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[   1 ] = tChild->get_basis(  19 );
                tNodes[   3 ] = tChild->get_basis(   2 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[   5 ] = tChild->get_basis(  18 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  49 ] = tChild->get_basis(  22 );
                tNodes[  50 ] = tChild->get_basis(  49 );
                tNodes[  51 ] = tChild->get_basis(  48 );
                tNodes[  52 ] = tChild->get_basis(  20 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  53 ] = tChild->get_basis(  49 );
                tNodes[  54 ] = tChild->get_basis(  48 );
                tNodes[  55 ] = tChild->get_basis(  20 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  99 ] = tChild->get_basis(  50 );
                tNodes[ 101 ] = tChild->get_basis(  21 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 103 ] = tChild->get_basis(  51 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 147 ] = tChild->get_basis(   7 );
                tNodes[ 148 ] = tChild->get_basis(  31 );
                tNodes[ 149 ] = tChild->get_basis(  30 );
                tNodes[ 150 ] = tChild->get_basis(   6 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 151 ] = tChild->get_basis(  31 );
                tNodes[ 152 ] = tChild->get_basis(  30 );
                tNodes[ 153 ] = tChild->get_basis(   6 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 197 ] = tChild->get_basis(  49 );
                tNodes[ 199 ] = tChild->get_basis(  20 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 201 ] = tChild->get_basis(  48 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 245 ] = tChild->get_basis(  23 );
                tNodes[ 246 ] = tChild->get_basis(  50 );
                tNodes[ 247 ] = tChild->get_basis(  51 );
                tNodes[ 248 ] = tChild->get_basis(  21 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 249 ] = tChild->get_basis(  50 );
                tNodes[ 250 ] = tChild->get_basis(  51 );
                tNodes[ 251 ] = tChild->get_basis(  21 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 295 ] = tChild->get_basis(  31 );
                tNodes[ 297 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 299 ] = tChild->get_basis(  30 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

        // test if neighbor 1 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 1 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[  13 ] = tChild->get_basis(  10 );
                tNodes[  27 ] = tChild->get_basis(   3 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  41 ] = tChild->get_basis(  11 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[  55 ] = tChild->get_basis(  12 );
                tNodes[  62 ] = tChild->get_basis(  40 );
                tNodes[  69 ] = tChild->get_basis(  43 );
                tNodes[  76 ] = tChild->get_basis(  22 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  83 ] = tChild->get_basis(  40 );
                tNodes[  90 ] = tChild->get_basis(  43 );
                tNodes[  97 ] = tChild->get_basis(  22 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 111 ] = tChild->get_basis(  41 );
                tNodes[ 125 ] = tChild->get_basis(  23 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 139 ] = tChild->get_basis(  42 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 153 ] = tChild->get_basis(   4 );
                tNodes[ 160 ] = tChild->get_basis(  26 );
                tNodes[ 167 ] = tChild->get_basis(  27 );
                tNodes[ 174 ] = tChild->get_basis(   7 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 181 ] = tChild->get_basis(  26 );
                tNodes[ 188 ] = tChild->get_basis(  27 );
                tNodes[ 195 ] = tChild->get_basis(   7 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 209 ] = tChild->get_basis(  40 );
                tNodes[ 223 ] = tChild->get_basis(  22 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 237 ] = tChild->get_basis(  43 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 251 ] = tChild->get_basis(  13 );
                tNodes[ 258 ] = tChild->get_basis(  41 );
                tNodes[ 265 ] = tChild->get_basis(  42 );
                tNodes[ 272 ] = tChild->get_basis(  23 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 279 ] = tChild->get_basis(  41 );
                tNodes[ 286 ] = tChild->get_basis(  42 );
                tNodes[ 293 ] = tChild->get_basis(  23 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 307 ] = tChild->get_basis(  26 );
                tNodes[ 321 ] = tChild->get_basis(   7 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 335 ] = tChild->get_basis(  27 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

        // test if neighbor 2 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 2 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[  43 ] = tChild->get_basis(   8 );
                tNodes[  45 ] = tChild->get_basis(   1 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  47 ] = tChild->get_basis(   9 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[  91 ] = tChild->get_basis(  12 );
                tNodes[  92 ] = tChild->get_basis(  36 );
                tNodes[  93 ] = tChild->get_basis(  37 );
                tNodes[  94 ] = tChild->get_basis(  16 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  95 ] = tChild->get_basis(  36 );
                tNodes[  96 ] = tChild->get_basis(  37 );
                tNodes[  97 ] = tChild->get_basis(  16 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 141 ] = tChild->get_basis(  39 );
                tNodes[ 143 ] = tChild->get_basis(  17 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 145 ] = tChild->get_basis(  38 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 189 ] = tChild->get_basis(   4 );
                tNodes[ 190 ] = tChild->get_basis(  24 );
                tNodes[ 191 ] = tChild->get_basis(  25 );
                tNodes[ 192 ] = tChild->get_basis(   5 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 193 ] = tChild->get_basis(  24 );
                tNodes[ 194 ] = tChild->get_basis(  25 );
                tNodes[ 195 ] = tChild->get_basis(   5 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 239 ] = tChild->get_basis(  36 );
                tNodes[ 241 ] = tChild->get_basis(  16 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 243 ] = tChild->get_basis(  37 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 287 ] = tChild->get_basis(  13 );
                tNodes[ 288 ] = tChild->get_basis(  39 );
                tNodes[ 289 ] = tChild->get_basis(  38 );
                tNodes[ 290 ] = tChild->get_basis(  17 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 291 ] = tChild->get_basis(  39 );
                tNodes[ 292 ] = tChild->get_basis(  38 );
                tNodes[ 293 ] = tChild->get_basis(  17 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 337 ] = tChild->get_basis(  24 );
                tNodes[ 339 ] = tChild->get_basis(   5 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 341 ] = tChild->get_basis(  25 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

        // test if neighbor 3 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 3 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 1
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[   7 ] = tChild->get_basis(  14 );
                tNodes[  21 ] = tChild->get_basis(   2 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  35 ] = tChild->get_basis(  15 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  49 ] = tChild->get_basis(  16 );
                tNodes[  56 ] = tChild->get_basis(  44 );
                tNodes[  63 ] = tChild->get_basis(  45 );
                tNodes[  70 ] = tChild->get_basis(  20 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  77 ] = tChild->get_basis(  44 );
                tNodes[  84 ] = tChild->get_basis(  45 );
                tNodes[  91 ] = tChild->get_basis(  20 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 105 ] = tChild->get_basis(  47 );
                tNodes[ 119 ] = tChild->get_basis(  21 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 133 ] = tChild->get_basis(  46 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 147 ] = tChild->get_basis(   5 );
                tNodes[ 154 ] = tChild->get_basis(  28 );
                tNodes[ 161 ] = tChild->get_basis(  29 );
                tNodes[ 168 ] = tChild->get_basis(   6 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 175 ] = tChild->get_basis(  28 );
                tNodes[ 182 ] = tChild->get_basis(  29 );
                tNodes[ 189 ] = tChild->get_basis(   6 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 203 ] = tChild->get_basis(  44 );
                tNodes[ 217 ] = tChild->get_basis(  20 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 231 ] = tChild->get_basis(  45 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 245 ] = tChild->get_basis(  17 );
                tNodes[ 252 ] = tChild->get_basis(  47 );
                tNodes[ 259 ] = tChild->get_basis(  46 );
                tNodes[ 266 ] = tChild->get_basis(  21 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 273 ] = tChild->get_basis(  47 );
                tNodes[ 280 ] = tChild->get_basis(  46 );
                tNodes[ 287 ] = tChild->get_basis(  21 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 301 ] = tChild->get_basis(  28 );
                tNodes[ 315 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 329 ] = tChild->get_basis(  29 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

        // test if neighbor 4 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 4 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 4
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[   1 ] = tChild->get_basis(  24 );
                tNodes[   3 ] = tChild->get_basis(   5 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[   5 ] = tChild->get_basis(  25 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[   7 ] = tChild->get_basis(  26 );
                tNodes[   8 ] = tChild->get_basis(  52 );
                tNodes[   9 ] = tChild->get_basis(  53 );
                tNodes[  10 ] = tChild->get_basis(  28 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[  11 ] = tChild->get_basis(  52 );
                tNodes[  12 ] = tChild->get_basis(  53 );
                tNodes[  13 ] = tChild->get_basis(  28 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  15 ] = tChild->get_basis(  55 );
                tNodes[  17 ] = tChild->get_basis(  29 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[  19 ] = tChild->get_basis(  54 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  21 ] = tChild->get_basis(   7 );
                tNodes[  22 ] = tChild->get_basis(  31 );
                tNodes[  23 ] = tChild->get_basis(  30 );
                tNodes[  24 ] = tChild->get_basis(   6 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[  25 ] = tChild->get_basis(  31 );
                tNodes[  26 ] = tChild->get_basis(  30 );
                tNodes[  27 ] = tChild->get_basis(   6 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[  29 ] = tChild->get_basis(  52 );
                tNodes[  31 ] = tChild->get_basis(  28 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[  33 ] = tChild->get_basis(  53 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[  35 ] = tChild->get_basis(  27 );
                tNodes[  36 ] = tChild->get_basis(  55 );
                tNodes[  37 ] = tChild->get_basis(  54 );
                tNodes[  38 ] = tChild->get_basis(  29 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[  39 ] = tChild->get_basis(  55 );
                tNodes[  40 ] = tChild->get_basis(  54 );
                tNodes[  41 ] = tChild->get_basis(  29 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[  43 ] = tChild->get_basis(  31 );
                tNodes[  45 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[  47 ] = tChild->get_basis(  30 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

        // test if neighbor 5 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on face 5 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 295 ] = tChild->get_basis(   8 );
                tNodes[ 297 ] = tChild->get_basis(   1 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 299 ] = tChild->get_basis(   9 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 301 ] = tChild->get_basis(  10 );
                tNodes[ 302 ] = tChild->get_basis(  32 );
                tNodes[ 303 ] = tChild->get_basis(  35 );
                tNodes[ 304 ] = tChild->get_basis(  14 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 305 ] = tChild->get_basis(  32 );
                tNodes[ 306 ] = tChild->get_basis(  35 );
                tNodes[ 307 ] = tChild->get_basis(  14 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 309 ] = tChild->get_basis(  33 );
                tNodes[ 311 ] = tChild->get_basis(  15 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 313 ] = tChild->get_basis(  34 );

                // get pointer to child 0
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 315 ] = tChild->get_basis(   3 );
                tNodes[ 316 ] = tChild->get_basis(  19 );
                tNodes[ 317 ] = tChild->get_basis(  18 );
                tNodes[ 318 ] = tChild->get_basis(   2 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 319 ] = tChild->get_basis(  19 );
                tNodes[ 320 ] = tChild->get_basis(  18 );
                tNodes[ 321 ] = tChild->get_basis(   2 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 323 ] = tChild->get_basis(  32 );
                tNodes[ 325 ] = tChild->get_basis(  14 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 327 ] = tChild->get_basis(  35 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 329 ] = tChild->get_basis(  11 );
                tNodes[ 330 ] = tChild->get_basis(  33 );
                tNodes[ 331 ] = tChild->get_basis(  34 );
                tNodes[ 332 ] = tChild->get_basis(  15 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 333 ] = tChild->get_basis(  33 );
                tNodes[ 334 ] = tChild->get_basis(  34 );
                tNodes[ 335 ] = tChild->get_basis(  15 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 337 ] = tChild->get_basis(  19 );
                tNodes[ 339 ] = tChild->get_basis(   2 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 341 ] = tChild->get_basis(  18 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

        // test if neighbor 6 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 0 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 6
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[   1 ] = tChild->get_basis(  31 );
                tNodes[   3 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[   5 ] = tChild->get_basis(  30 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

        // test if neighbor 7 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 1 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 4
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  13 ] = tChild->get_basis(  26 );
                tNodes[  27 ] = tChild->get_basis(   7 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[  41 ] = tChild->get_basis(  27 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 8 );

        // test if neighbor 8 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 2 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 4
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  43 ] = tChild->get_basis(  24 );
                tNodes[  45 ] = tChild->get_basis(   5 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[  47 ] = tChild->get_basis(  25 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 9 );

        // test if neighbor 9 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 3 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 5
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[   7 ] = tChild->get_basis(  28 );
                tNodes[  21 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[  35 ] = tChild->get_basis(  29 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 10 );

        // test if neighbor 10 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 4 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 3
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  49 ] = tChild->get_basis(  20 );
                tNodes[ 147 ] = tChild->get_basis(   6 );

                // get pointer to child 7
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 7 )->get_memory_index() );

                tNodes[ 245 ] = tChild->get_basis(  21 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 11 );

        // test if neighbor 11 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 5 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 2
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  55 ] = tChild->get_basis(  22 );
                tNodes[ 153 ] = tChild->get_basis(   7 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[ 251 ] = tChild->get_basis(  23 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 12 );

        // test if neighbor 12 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 6 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[  97 ] = tChild->get_basis(  12 );
                tNodes[ 195 ] = tChild->get_basis(   4 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[ 293 ] = tChild->get_basis(  13 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 13 );

        // test if neighbor 13 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 7 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 1
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  91 ] = tChild->get_basis(  16 );
                tNodes[ 189 ] = tChild->get_basis(   5 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[ 287 ] = tChild->get_basis(  17 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 14 );

        // test if neighbor 14 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 8 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 2
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 295 ] = tChild->get_basis(  19 );
                tNodes[ 297 ] = tChild->get_basis(   2 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 299 ] = tChild->get_basis(  18 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 15 );

        // test if neighbor 15 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 9 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 307 ] = tChild->get_basis(  10 );
                tNodes[ 321 ] = tChild->get_basis(   3 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[ 335 ] = tChild->get_basis(  11 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 16 );

        // test if neighbor 16 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 10 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 0
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 0 )->get_memory_index() );

                tNodes[ 337 ] = tChild->get_basis(   8 );
                tNodes[ 339 ] = tChild->get_basis(   1 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 341 ] = tChild->get_basis(   9 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 17 );

        // test if neighbor 17 exists
        if ( tNeighbor != nullptr )
        {
            // test if nodes on edge 11 exist
            if ( tNeighbor->children_have_basis() )
            {
                // get pointer to background element
                Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                // get pointer to child 1
                Element* tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[ 301 ] = tChild->get_basis(  14 );
                tNodes[ 315 ] = tChild->get_basis(   2 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[ 329 ] = tChild->get_basis(  15 );
            }
        }

        // level of node
        auto tLevel = mElement->get_level() + 1;

        // owner of element
        auto tOwner = mElement->get_owner();

        // get position of element
        const luint * tElIJK = mElement->get_ijk();

        // anchor point of nodes
        luint tAnchor[ 3 ];
        tAnchor[ 0 ] = 6 * tElIJK[ 0 ];
        tAnchor[ 1 ] = 6 * tElIJK[ 1 ];
        tAnchor[ 2 ] = 6 * tElIJK[ 2 ];

        // array containing node position;
        luint tIJK[ 3 ] = { 0, 0, 0 };

        // test if node 1 exists
        if ( tNodes[ 1 ] == nullptr )
        {
             // calculate position of node 1
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 1
             tNodes[ 1 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 3 exists
        if ( tNodes[ 3 ] == nullptr )
        {
             // calculate position of node 3
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 3
             tNodes[ 3 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 5 exists
        if ( tNodes[ 5 ] == nullptr )
        {
             // calculate position of node 5
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 5
             tNodes[ 5 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 7 exists
        if ( tNodes[ 7 ] == nullptr )
        {
             // calculate position of node 7
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 7
             tNodes[ 7 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 8 exists
        if ( tNodes[ 8 ] == nullptr )
        {
             // calculate position of node 8
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 8
             tNodes[ 8 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 9 exists
        if ( tNodes[ 9 ] == nullptr )
        {
             // calculate position of node 9
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 9
             tNodes[ 9 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 10 exists
        if ( tNodes[ 10 ] == nullptr )
        {
             // calculate position of node 10
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 10
             tNodes[ 10 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 11 exists
        if ( tNodes[ 11 ] == nullptr )
        {
             // calculate position of node 11
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 11
             tNodes[ 11 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 12 exists
        if ( tNodes[ 12 ] == nullptr )
        {
             // calculate position of node 12
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 12
             tNodes[ 12 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 13 exists
        if ( tNodes[ 13 ] == nullptr )
        {
             // calculate position of node 13
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 13
             tNodes[ 13 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 15 exists
        if ( tNodes[ 15 ] == nullptr )
        {
             // calculate position of node 15
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 15
             tNodes[ 15 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 17 exists
        if ( tNodes[ 17 ] == nullptr )
        {
             // calculate position of node 17
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 17
             tNodes[ 17 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 19 exists
        if ( tNodes[ 19 ] == nullptr )
        {
             // calculate position of node 19
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 19
             tNodes[ 19 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 21 exists
        if ( tNodes[ 21 ] == nullptr )
        {
             // calculate position of node 21
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 21
             tNodes[ 21 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 22 exists
        if ( tNodes[ 22 ] == nullptr )
        {
             // calculate position of node 22
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 22
             tNodes[ 22 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 23 exists
        if ( tNodes[ 23 ] == nullptr )
        {
             // calculate position of node 23
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 23
             tNodes[ 23 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 24 exists
        if ( tNodes[ 24 ] == nullptr )
        {
             // calculate position of node 24
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 24
             tNodes[ 24 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 25 exists
        if ( tNodes[ 25 ] == nullptr )
        {
             // calculate position of node 25
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 25
             tNodes[ 25 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 26 exists
        if ( tNodes[ 26 ] == nullptr )
        {
             // calculate position of node 26
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 26
             tNodes[ 26 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 27 exists
        if ( tNodes[ 27 ] == nullptr )
        {
             // calculate position of node 27
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 27
             tNodes[ 27 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 29 exists
        if ( tNodes[ 29 ] == nullptr )
        {
             // calculate position of node 29
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 29
             tNodes[ 29 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 31 exists
        if ( tNodes[ 31 ] == nullptr )
        {
             // calculate position of node 31
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 31
             tNodes[ 31 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 33 exists
        if ( tNodes[ 33 ] == nullptr )
        {
             // calculate position of node 33
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 33
             tNodes[ 33 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 35 exists
        if ( tNodes[ 35 ] == nullptr )
        {
             // calculate position of node 35
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 35
             tNodes[ 35 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 36 exists
        if ( tNodes[ 36 ] == nullptr )
        {
             // calculate position of node 36
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 36
             tNodes[ 36 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 37 exists
        if ( tNodes[ 37 ] == nullptr )
        {
             // calculate position of node 37
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 37
             tNodes[ 37 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 38 exists
        if ( tNodes[ 38 ] == nullptr )
        {
             // calculate position of node 38
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 38
             tNodes[ 38 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 39 exists
        if ( tNodes[ 39 ] == nullptr )
        {
             // calculate position of node 39
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 39
             tNodes[ 39 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 40 exists
        if ( tNodes[ 40 ] == nullptr )
        {
             // calculate position of node 40
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 40
             tNodes[ 40 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 41 exists
        if ( tNodes[ 41 ] == nullptr )
        {
             // calculate position of node 41
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 41
             tNodes[ 41 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 43 exists
        if ( tNodes[ 43 ] == nullptr )
        {
             // calculate position of node 43
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 43
             tNodes[ 43 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 45 exists
        if ( tNodes[ 45 ] == nullptr )
        {
             // calculate position of node 45
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 45
             tNodes[ 45 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 47 exists
        if ( tNodes[ 47 ] == nullptr )
        {
             // calculate position of node 47
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 47
             tNodes[ 47 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 49 exists
        if ( tNodes[ 49 ] == nullptr )
        {
             // calculate position of node 49
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 49
             tNodes[ 49 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 50 exists
        if ( tNodes[ 50 ] == nullptr )
        {
             // calculate position of node 50
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 50
             tNodes[ 50 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 51 exists
        if ( tNodes[ 51 ] == nullptr )
        {
             // calculate position of node 51
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 51
             tNodes[ 51 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 52 exists
        if ( tNodes[ 52 ] == nullptr )
        {
             // calculate position of node 52
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 52
             tNodes[ 52 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 53 exists
        if ( tNodes[ 53 ] == nullptr )
        {
             // calculate position of node 53
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 53
             tNodes[ 53 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 54 exists
        if ( tNodes[ 54 ] == nullptr )
        {
             // calculate position of node 54
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 54
             tNodes[ 54 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 55 exists
        if ( tNodes[ 55 ] == nullptr )
        {
             // calculate position of node 55
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 55
             tNodes[ 55 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 56 exists
        if ( tNodes[ 56 ] == nullptr )
        {
             // calculate position of node 56
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 56
             tNodes[ 56 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 57
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 57
         tNodes[ 57 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 58
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 58
         tNodes[ 58 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 59
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 59
         tNodes[ 59 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 60
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 60
         tNodes[ 60 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 61
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 61
         tNodes[ 61 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 62 exists
        if ( tNodes[ 62 ] == nullptr )
        {
             // calculate position of node 62
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 62
             tNodes[ 62 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 63 exists
        if ( tNodes[ 63 ] == nullptr )
        {
             // calculate position of node 63
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 63
             tNodes[ 63 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 64
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 64
         tNodes[ 64 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 65
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 65
         tNodes[ 65 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 66
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 66
         tNodes[ 66 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 67
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 67
         tNodes[ 67 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 68
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 68
         tNodes[ 68 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 69 exists
        if ( tNodes[ 69 ] == nullptr )
        {
             // calculate position of node 69
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 69
             tNodes[ 69 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 70 exists
        if ( tNodes[ 70 ] == nullptr )
        {
             // calculate position of node 70
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 70
             tNodes[ 70 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 71
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 71
         tNodes[ 71 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 72
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 72
         tNodes[ 72 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 73
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 73
         tNodes[ 73 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 74
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 74
         tNodes[ 74 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 75
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 75
         tNodes[ 75 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 76 exists
        if ( tNodes[ 76 ] == nullptr )
        {
             // calculate position of node 76
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 76
             tNodes[ 76 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 77 exists
        if ( tNodes[ 77 ] == nullptr )
        {
             // calculate position of node 77
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 77
             tNodes[ 77 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 78
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 78
         tNodes[ 78 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 79
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 79
         tNodes[ 79 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 80
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 80
         tNodes[ 80 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 81
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 81
         tNodes[ 81 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 82
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 82
         tNodes[ 82 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 83 exists
        if ( tNodes[ 83 ] == nullptr )
        {
             // calculate position of node 83
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 83
             tNodes[ 83 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 84 exists
        if ( tNodes[ 84 ] == nullptr )
        {
             // calculate position of node 84
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 84
             tNodes[ 84 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 85
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 85
         tNodes[ 85 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 86
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 86
         tNodes[ 86 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 87
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 87
         tNodes[ 87 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 88
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 88
         tNodes[ 88 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 89
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 89
         tNodes[ 89 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 90 exists
        if ( tNodes[ 90 ] == nullptr )
        {
             // calculate position of node 90
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 90
             tNodes[ 90 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 91 exists
        if ( tNodes[ 91 ] == nullptr )
        {
             // calculate position of node 91
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 91
             tNodes[ 91 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 92 exists
        if ( tNodes[ 92 ] == nullptr )
        {
             // calculate position of node 92
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 92
             tNodes[ 92 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 93 exists
        if ( tNodes[ 93 ] == nullptr )
        {
             // calculate position of node 93
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 93
             tNodes[ 93 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 94 exists
        if ( tNodes[ 94 ] == nullptr )
        {
             // calculate position of node 94
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 94
             tNodes[ 94 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 95 exists
        if ( tNodes[ 95 ] == nullptr )
        {
             // calculate position of node 95
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 95
             tNodes[ 95 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 96 exists
        if ( tNodes[ 96 ] == nullptr )
        {
             // calculate position of node 96
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 96
             tNodes[ 96 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 97 exists
        if ( tNodes[ 97 ] == nullptr )
        {
             // calculate position of node 97
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 97
             tNodes[ 97 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 99 exists
        if ( tNodes[ 99 ] == nullptr )
        {
             // calculate position of node 99
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 99
             tNodes[ 99 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 101 exists
        if ( tNodes[ 101 ] == nullptr )
        {
             // calculate position of node 101
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 101
             tNodes[ 101 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 103 exists
        if ( tNodes[ 103 ] == nullptr )
        {
             // calculate position of node 103
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 103
             tNodes[ 103 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 105 exists
        if ( tNodes[ 105 ] == nullptr )
        {
             // calculate position of node 105
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 105
             tNodes[ 105 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 106
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 106
         tNodes[ 106 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 107
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 107
         tNodes[ 107 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 108
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 108
         tNodes[ 108 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 109
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 109
         tNodes[ 109 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 110
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 110
         tNodes[ 110 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 111 exists
        if ( tNodes[ 111 ] == nullptr )
        {
             // calculate position of node 111
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 111
             tNodes[ 111 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 113
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 113
         tNodes[ 113 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 115
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 115
         tNodes[ 115 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 117
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 117
         tNodes[ 117 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 119 exists
        if ( tNodes[ 119 ] == nullptr )
        {
             // calculate position of node 119
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 119
             tNodes[ 119 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 120
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 120
         tNodes[ 120 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 121
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 121
         tNodes[ 121 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 122
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 122
         tNodes[ 122 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 123
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 123
         tNodes[ 123 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 124
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 124
         tNodes[ 124 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 125 exists
        if ( tNodes[ 125 ] == nullptr )
        {
             // calculate position of node 125
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 125
             tNodes[ 125 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 127
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 127
         tNodes[ 127 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 129
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 129
         tNodes[ 129 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 131
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 131
         tNodes[ 131 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 133 exists
        if ( tNodes[ 133 ] == nullptr )
        {
             // calculate position of node 133
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 133
             tNodes[ 133 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 134
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 134
         tNodes[ 134 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 135
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 135
         tNodes[ 135 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 136
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 136
         tNodes[ 136 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 137
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 137
         tNodes[ 137 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 138
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 2;

         // create node 138
         tNodes[ 138 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 139 exists
        if ( tNodes[ 139 ] == nullptr )
        {
             // calculate position of node 139
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 139
             tNodes[ 139 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 141 exists
        if ( tNodes[ 141 ] == nullptr )
        {
             // calculate position of node 141
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 141
             tNodes[ 141 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 143 exists
        if ( tNodes[ 143 ] == nullptr )
        {
             // calculate position of node 143
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 143
             tNodes[ 143 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 145 exists
        if ( tNodes[ 145 ] == nullptr )
        {
             // calculate position of node 145
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 145
             tNodes[ 145 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 147 exists
        if ( tNodes[ 147 ] == nullptr )
        {
             // calculate position of node 147
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 147
             tNodes[ 147 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 148 exists
        if ( tNodes[ 148 ] == nullptr )
        {
             // calculate position of node 148
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 148
             tNodes[ 148 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 149 exists
        if ( tNodes[ 149 ] == nullptr )
        {
             // calculate position of node 149
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 149
             tNodes[ 149 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 150 exists
        if ( tNodes[ 150 ] == nullptr )
        {
             // calculate position of node 150
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 150
             tNodes[ 150 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 151 exists
        if ( tNodes[ 151 ] == nullptr )
        {
             // calculate position of node 151
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 151
             tNodes[ 151 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 152 exists
        if ( tNodes[ 152 ] == nullptr )
        {
             // calculate position of node 152
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 152
             tNodes[ 152 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 153 exists
        if ( tNodes[ 153 ] == nullptr )
        {
             // calculate position of node 153
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 153
             tNodes[ 153 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 154 exists
        if ( tNodes[ 154 ] == nullptr )
        {
             // calculate position of node 154
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 154
             tNodes[ 154 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 155
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 155
         tNodes[ 155 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 156
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 156
         tNodes[ 156 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 157
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 157
         tNodes[ 157 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 158
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 158
         tNodes[ 158 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 159
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 159
         tNodes[ 159 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 160 exists
        if ( tNodes[ 160 ] == nullptr )
        {
             // calculate position of node 160
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 160
             tNodes[ 160 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 161 exists
        if ( tNodes[ 161 ] == nullptr )
        {
             // calculate position of node 161
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 161
             tNodes[ 161 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 162
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 162
         tNodes[ 162 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 163
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 163
         tNodes[ 163 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 164
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 164
         tNodes[ 164 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 165
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 165
         tNodes[ 165 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 166
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 166
         tNodes[ 166 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 167 exists
        if ( tNodes[ 167 ] == nullptr )
        {
             // calculate position of node 167
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 167
             tNodes[ 167 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 168 exists
        if ( tNodes[ 168 ] == nullptr )
        {
             // calculate position of node 168
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 168
             tNodes[ 168 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 169
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 169
         tNodes[ 169 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 170
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 170
         tNodes[ 170 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 171
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 171
         tNodes[ 171 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 172
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 172
         tNodes[ 172 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 173
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 173
         tNodes[ 173 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 174 exists
        if ( tNodes[ 174 ] == nullptr )
        {
             // calculate position of node 174
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 174
             tNodes[ 174 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 175 exists
        if ( tNodes[ 175 ] == nullptr )
        {
             // calculate position of node 175
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 175
             tNodes[ 175 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 176
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 176
         tNodes[ 176 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 177
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 177
         tNodes[ 177 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 178
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 178
         tNodes[ 178 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 179
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 179
         tNodes[ 179 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 180
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 180
         tNodes[ 180 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 181 exists
        if ( tNodes[ 181 ] == nullptr )
        {
             // calculate position of node 181
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 181
             tNodes[ 181 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 182 exists
        if ( tNodes[ 182 ] == nullptr )
        {
             // calculate position of node 182
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 182
             tNodes[ 182 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 183
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 183
         tNodes[ 183 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 184
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 184
         tNodes[ 184 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 185
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 185
         tNodes[ 185 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 186
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 186
         tNodes[ 186 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 187
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 3;

         // create node 187
         tNodes[ 187 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 188 exists
        if ( tNodes[ 188 ] == nullptr )
        {
             // calculate position of node 188
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 188
             tNodes[ 188 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 189 exists
        if ( tNodes[ 189 ] == nullptr )
        {
             // calculate position of node 189
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 189
             tNodes[ 189 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 190 exists
        if ( tNodes[ 190 ] == nullptr )
        {
             // calculate position of node 190
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 190
             tNodes[ 190 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 191 exists
        if ( tNodes[ 191 ] == nullptr )
        {
             // calculate position of node 191
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 191
             tNodes[ 191 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 192 exists
        if ( tNodes[ 192 ] == nullptr )
        {
             // calculate position of node 192
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 192
             tNodes[ 192 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 193 exists
        if ( tNodes[ 193 ] == nullptr )
        {
             // calculate position of node 193
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 193
             tNodes[ 193 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 194 exists
        if ( tNodes[ 194 ] == nullptr )
        {
             // calculate position of node 194
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 194
             tNodes[ 194 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 195 exists
        if ( tNodes[ 195 ] == nullptr )
        {
             // calculate position of node 195
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 195
             tNodes[ 195 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 197 exists
        if ( tNodes[ 197 ] == nullptr )
        {
             // calculate position of node 197
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 197
             tNodes[ 197 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 199 exists
        if ( tNodes[ 199 ] == nullptr )
        {
             // calculate position of node 199
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 199
             tNodes[ 199 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 201 exists
        if ( tNodes[ 201 ] == nullptr )
        {
             // calculate position of node 201
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 201
             tNodes[ 201 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 203 exists
        if ( tNodes[ 203 ] == nullptr )
        {
             // calculate position of node 203
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 203
             tNodes[ 203 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 204
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 204
         tNodes[ 204 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 205
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 205
         tNodes[ 205 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 206
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 206
         tNodes[ 206 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 207
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 207
         tNodes[ 207 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 208
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 208
         tNodes[ 208 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 209 exists
        if ( tNodes[ 209 ] == nullptr )
        {
             // calculate position of node 209
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 209
             tNodes[ 209 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 211
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 211
         tNodes[ 211 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 213
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 213
         tNodes[ 213 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 215
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 215
         tNodes[ 215 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 217 exists
        if ( tNodes[ 217 ] == nullptr )
        {
             // calculate position of node 217
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 217
             tNodes[ 217 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 218
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 218
         tNodes[ 218 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 219
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 219
         tNodes[ 219 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 220
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 220
         tNodes[ 220 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 221
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 221
         tNodes[ 221 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 222
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 222
         tNodes[ 222 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 223 exists
        if ( tNodes[ 223 ] == nullptr )
        {
             // calculate position of node 223
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 223
             tNodes[ 223 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 225
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 225
         tNodes[ 225 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 227
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 227
         tNodes[ 227 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 229
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 229
         tNodes[ 229 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 231 exists
        if ( tNodes[ 231 ] == nullptr )
        {
             // calculate position of node 231
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 231
             tNodes[ 231 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 232
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 232
         tNodes[ 232 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 233
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 233
         tNodes[ 233 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 234
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 234
         tNodes[ 234 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 235
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 235
         tNodes[ 235 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 236
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 4;

         // create node 236
         tNodes[ 236 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 237 exists
        if ( tNodes[ 237 ] == nullptr )
        {
             // calculate position of node 237
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 237
             tNodes[ 237 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 239 exists
        if ( tNodes[ 239 ] == nullptr )
        {
             // calculate position of node 239
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 239
             tNodes[ 239 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 241 exists
        if ( tNodes[ 241 ] == nullptr )
        {
             // calculate position of node 241
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 241
             tNodes[ 241 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 243 exists
        if ( tNodes[ 243 ] == nullptr )
        {
             // calculate position of node 243
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 4;

             // create node 243
             tNodes[ 243 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 245 exists
        if ( tNodes[ 245 ] == nullptr )
        {
             // calculate position of node 245
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 245
             tNodes[ 245 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 246 exists
        if ( tNodes[ 246 ] == nullptr )
        {
             // calculate position of node 246
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 246
             tNodes[ 246 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 247 exists
        if ( tNodes[ 247 ] == nullptr )
        {
             // calculate position of node 247
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 247
             tNodes[ 247 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 248 exists
        if ( tNodes[ 248 ] == nullptr )
        {
             // calculate position of node 248
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 248
             tNodes[ 248 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 249 exists
        if ( tNodes[ 249 ] == nullptr )
        {
             // calculate position of node 249
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 249
             tNodes[ 249 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 250 exists
        if ( tNodes[ 250 ] == nullptr )
        {
             // calculate position of node 250
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 250
             tNodes[ 250 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 251 exists
        if ( tNodes[ 251 ] == nullptr )
        {
             // calculate position of node 251
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 251
             tNodes[ 251 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 252 exists
        if ( tNodes[ 252 ] == nullptr )
        {
             // calculate position of node 252
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 252
             tNodes[ 252 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 253
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 253
         tNodes[ 253 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 254
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 254
         tNodes[ 254 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 255
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 255
         tNodes[ 255 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 256
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 256
         tNodes[ 256 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 257
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 257
         tNodes[ 257 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 258 exists
        if ( tNodes[ 258 ] == nullptr )
        {
             // calculate position of node 258
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 258
             tNodes[ 258 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 259 exists
        if ( tNodes[ 259 ] == nullptr )
        {
             // calculate position of node 259
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 259
             tNodes[ 259 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 260
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 260
         tNodes[ 260 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 261
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 261
         tNodes[ 261 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 262
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 262
         tNodes[ 262 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 263
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 263
         tNodes[ 263 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 264
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 2;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 264
         tNodes[ 264 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 265 exists
        if ( tNodes[ 265 ] == nullptr )
        {
             // calculate position of node 265
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 265
             tNodes[ 265 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 266 exists
        if ( tNodes[ 266 ] == nullptr )
        {
             // calculate position of node 266
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 266
             tNodes[ 266 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 267
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 267
         tNodes[ 267 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 268
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 268
         tNodes[ 268 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 269
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 269
         tNodes[ 269 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 270
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 270
         tNodes[ 270 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 271
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 3;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 271
         tNodes[ 271 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 272 exists
        if ( tNodes[ 272 ] == nullptr )
        {
             // calculate position of node 272
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 272
             tNodes[ 272 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 273 exists
        if ( tNodes[ 273 ] == nullptr )
        {
             // calculate position of node 273
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 273
             tNodes[ 273 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 274
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 274
         tNodes[ 274 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 275
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 275
         tNodes[ 275 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 276
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 276
         tNodes[ 276 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 277
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 277
         tNodes[ 277 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 278
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 4;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 278
         tNodes[ 278 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 279 exists
        if ( tNodes[ 279 ] == nullptr )
        {
             // calculate position of node 279
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 279
             tNodes[ 279 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 280 exists
        if ( tNodes[ 280 ] == nullptr )
        {
             // calculate position of node 280
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 280
             tNodes[ 280 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 281
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 281
         tNodes[ 281 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 282
         tIJK[ 0 ] = tAnchor[ 0 ] + 2;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 282
         tNodes[ 282 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 283
         tIJK[ 0 ] = tAnchor[ 0 ] + 3;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 283
         tNodes[ 283 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 284
         tIJK[ 0 ] = tAnchor[ 0 ] + 4;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 284
         tNodes[ 284 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 285
         tIJK[ 0 ] = tAnchor[ 0 ] + 5;
         tIJK[ 1 ] = tAnchor[ 1 ] + 5;
         tIJK[ 2 ] = tAnchor[ 2 ] + 5;

         // create node 285
         tNodes[ 285 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 286 exists
        if ( tNodes[ 286 ] == nullptr )
        {
             // calculate position of node 286
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 286
             tNodes[ 286 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 287 exists
        if ( tNodes[ 287 ] == nullptr )
        {
             // calculate position of node 287
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 287
             tNodes[ 287 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 288 exists
        if ( tNodes[ 288 ] == nullptr )
        {
             // calculate position of node 288
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 288
             tNodes[ 288 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 289 exists
        if ( tNodes[ 289 ] == nullptr )
        {
             // calculate position of node 289
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 289
             tNodes[ 289 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 290 exists
        if ( tNodes[ 290 ] == nullptr )
        {
             // calculate position of node 290
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 290
             tNodes[ 290 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 291 exists
        if ( tNodes[ 291 ] == nullptr )
        {
             // calculate position of node 291
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 291
             tNodes[ 291 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 292 exists
        if ( tNodes[ 292 ] == nullptr )
        {
             // calculate position of node 292
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 292
             tNodes[ 292 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 293 exists
        if ( tNodes[ 293 ] == nullptr )
        {
             // calculate position of node 293
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 5;

             // create node 293
             tNodes[ 293 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 295 exists
        if ( tNodes[ 295 ] == nullptr )
        {
             // calculate position of node 295
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 295
             tNodes[ 295 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 297 exists
        if ( tNodes[ 297 ] == nullptr )
        {
             // calculate position of node 297
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 297
             tNodes[ 297 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 299 exists
        if ( tNodes[ 299 ] == nullptr )
        {
             // calculate position of node 299
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 299
             tNodes[ 299 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 301 exists
        if ( tNodes[ 301 ] == nullptr )
        {
             // calculate position of node 301
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 301
             tNodes[ 301 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 302 exists
        if ( tNodes[ 302 ] == nullptr )
        {
             // calculate position of node 302
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 302
             tNodes[ 302 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 303 exists
        if ( tNodes[ 303 ] == nullptr )
        {
             // calculate position of node 303
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 303
             tNodes[ 303 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 304 exists
        if ( tNodes[ 304 ] == nullptr )
        {
             // calculate position of node 304
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 304
             tNodes[ 304 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 305 exists
        if ( tNodes[ 305 ] == nullptr )
        {
             // calculate position of node 305
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 305
             tNodes[ 305 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 306 exists
        if ( tNodes[ 306 ] == nullptr )
        {
             // calculate position of node 306
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 306
             tNodes[ 306 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 307 exists
        if ( tNodes[ 307 ] == nullptr )
        {
             // calculate position of node 307
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 307
             tNodes[ 307 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 309 exists
        if ( tNodes[ 309 ] == nullptr )
        {
             // calculate position of node 309
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 309
             tNodes[ 309 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 311 exists
        if ( tNodes[ 311 ] == nullptr )
        {
             // calculate position of node 311
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 311
             tNodes[ 311 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 313 exists
        if ( tNodes[ 313 ] == nullptr )
        {
             // calculate position of node 313
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 313
             tNodes[ 313 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 315 exists
        if ( tNodes[ 315 ] == nullptr )
        {
             // calculate position of node 315
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 315
             tNodes[ 315 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 316 exists
        if ( tNodes[ 316 ] == nullptr )
        {
             // calculate position of node 316
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 316
             tNodes[ 316 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 317 exists
        if ( tNodes[ 317 ] == nullptr )
        {
             // calculate position of node 317
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 317
             tNodes[ 317 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 318 exists
        if ( tNodes[ 318 ] == nullptr )
        {
             // calculate position of node 318
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 318
             tNodes[ 318 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 319 exists
        if ( tNodes[ 319 ] == nullptr )
        {
             // calculate position of node 319
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 319
             tNodes[ 319 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 320 exists
        if ( tNodes[ 320 ] == nullptr )
        {
             // calculate position of node 320
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 320
             tNodes[ 320 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 321 exists
        if ( tNodes[ 321 ] == nullptr )
        {
             // calculate position of node 321
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 321
             tNodes[ 321 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 323 exists
        if ( tNodes[ 323 ] == nullptr )
        {
             // calculate position of node 323
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 323
             tNodes[ 323 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 325 exists
        if ( tNodes[ 325 ] == nullptr )
        {
             // calculate position of node 325
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 325
             tNodes[ 325 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 327 exists
        if ( tNodes[ 327 ] == nullptr )
        {
             // calculate position of node 327
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 4;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 327
             tNodes[ 327 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 329 exists
        if ( tNodes[ 329 ] == nullptr )
        {
             // calculate position of node 329
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 329
             tNodes[ 329 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 330 exists
        if ( tNodes[ 330 ] == nullptr )
        {
             // calculate position of node 330
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 330
             tNodes[ 330 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 331 exists
        if ( tNodes[ 331 ] == nullptr )
        {
             // calculate position of node 331
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 331
             tNodes[ 331 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 332 exists
        if ( tNodes[ 332 ] == nullptr )
        {
             // calculate position of node 332
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 332
             tNodes[ 332 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 333 exists
        if ( tNodes[ 333 ] == nullptr )
        {
             // calculate position of node 333
             tIJK[ 0 ] = tAnchor[ 0 ] + 4;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 333
             tNodes[ 333 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 334 exists
        if ( tNodes[ 334 ] == nullptr )
        {
             // calculate position of node 334
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 334
             tNodes[ 334 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 335 exists
        if ( tNodes[ 335 ] == nullptr )
        {
             // calculate position of node 335
             tIJK[ 0 ] = tAnchor[ 0 ] + 6;
             tIJK[ 1 ] = tAnchor[ 1 ] + 5;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 335
             tNodes[ 335 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 337 exists
        if ( tNodes[ 337 ] == nullptr )
        {
             // calculate position of node 337
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 337
             tNodes[ 337 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 339 exists
        if ( tNodes[ 339 ] == nullptr )
        {
             // calculate position of node 339
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 339
             tNodes[ 339 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 341 exists
        if ( tNodes[ 341 ] == nullptr )
        {
             // calculate position of node 341
             tIJK[ 0 ] = tAnchor[ 0 ] + 5;
             tIJK[ 1 ] = tAnchor[ 1 ] + 6;
             tIJK[ 2 ] = tAnchor[ 2 ] + 6;

             // create node 341
             tNodes[ 341 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // pointer to child
         Element* tChild;

         // get pointer to child 0
         tChild = aAllElementsOnProc(
             mElement->get_child( 0 )->get_memory_index() );

         // init basis container for child 0
         tChild->init_basis_container();

         // link child 0 to nodes
         tChild->insert_basis(   0, tNodes[   0 ] );
         tChild->insert_basis(   1, tNodes[   3 ] );
         tChild->insert_basis(   2, tNodes[  24 ] );
         tChild->insert_basis(   3, tNodes[  21 ] );
         tChild->insert_basis(   4, tNodes[ 147 ] );
         tChild->insert_basis(   5, tNodes[ 150 ] );
         tChild->insert_basis(   6, tNodes[ 171 ] );
         tChild->insert_basis(   7, tNodes[ 168 ] );
         tChild->insert_basis(   8, tNodes[   1 ] );
         tChild->insert_basis(   9, tNodes[   2 ] );
         tChild->insert_basis(  10, tNodes[   7 ] );
         tChild->insert_basis(  11, tNodes[  14 ] );
         tChild->insert_basis(  12, tNodes[  49 ] );
         tChild->insert_basis(  13, tNodes[  98 ] );
         tChild->insert_basis(  14, tNodes[  10 ] );
         tChild->insert_basis(  15, tNodes[  17 ] );
         tChild->insert_basis(  16, tNodes[  52 ] );
         tChild->insert_basis(  17, tNodes[ 101 ] );
         tChild->insert_basis(  18, tNodes[  23 ] );
         tChild->insert_basis(  19, tNodes[  22 ] );
         tChild->insert_basis(  20, tNodes[  73 ] );
         tChild->insert_basis(  21, tNodes[ 122 ] );
         tChild->insert_basis(  22, tNodes[  70 ] );
         tChild->insert_basis(  23, tNodes[ 119 ] );
         tChild->insert_basis(  24, tNodes[ 148 ] );
         tChild->insert_basis(  25, tNodes[ 149 ] );
         tChild->insert_basis(  26, tNodes[ 154 ] );
         tChild->insert_basis(  27, tNodes[ 161 ] );
         tChild->insert_basis(  28, tNodes[ 157 ] );
         tChild->insert_basis(  29, tNodes[ 164 ] );
         tChild->insert_basis(  30, tNodes[ 170 ] );
         tChild->insert_basis(  31, tNodes[ 169 ] );
         tChild->insert_basis(  32, tNodes[   8 ] );
         tChild->insert_basis(  33, tNodes[  15 ] );
         tChild->insert_basis(  34, tNodes[  16 ] );
         tChild->insert_basis(  35, tNodes[   9 ] );
         tChild->insert_basis(  36, tNodes[  50 ] );
         tChild->insert_basis(  37, tNodes[  51 ] );
         tChild->insert_basis(  38, tNodes[ 100 ] );
         tChild->insert_basis(  39, tNodes[  99 ] );
         tChild->insert_basis(  40, tNodes[  56 ] );
         tChild->insert_basis(  41, tNodes[ 105 ] );
         tChild->insert_basis(  42, tNodes[ 112 ] );
         tChild->insert_basis(  43, tNodes[  63 ] );
         tChild->insert_basis(  44, tNodes[  59 ] );
         tChild->insert_basis(  45, tNodes[  66 ] );
         tChild->insert_basis(  46, tNodes[ 115 ] );
         tChild->insert_basis(  47, tNodes[ 108 ] );
         tChild->insert_basis(  48, tNodes[  72 ] );
         tChild->insert_basis(  49, tNodes[  71 ] );
         tChild->insert_basis(  50, tNodes[ 120 ] );
         tChild->insert_basis(  51, tNodes[ 121 ] );
         tChild->insert_basis(  52, tNodes[ 155 ] );
         tChild->insert_basis(  53, tNodes[ 156 ] );
         tChild->insert_basis(  54, tNodes[ 163 ] );
         tChild->insert_basis(  55, tNodes[ 162 ] );
         tChild->insert_basis(  56, tNodes[  57 ] );
         tChild->insert_basis(  57, tNodes[  58 ] );
         tChild->insert_basis(  58, tNodes[  65 ] );
         tChild->insert_basis(  59, tNodes[  64 ] );
         tChild->insert_basis(  60, tNodes[ 106 ] );
         tChild->insert_basis(  61, tNodes[ 107 ] );
         tChild->insert_basis(  62, tNodes[ 114 ] );
         tChild->insert_basis(  63, tNodes[ 113 ] );

         // get pointer to child 1
         tChild = aAllElementsOnProc(
             mElement->get_child( 1 )->get_memory_index() );

         // init basis container for child 1
         tChild->init_basis_container();

         // link child 1 to nodes
         tChild->insert_basis(   0, tNodes[   3 ] );
         tChild->insert_basis(   1, tNodes[   6 ] );
         tChild->insert_basis(   2, tNodes[  27 ] );
         tChild->insert_basis(   3, tNodes[  24 ] );
         tChild->insert_basis(   4, tNodes[ 150 ] );
         tChild->insert_basis(   5, tNodes[ 153 ] );
         tChild->insert_basis(   6, tNodes[ 174 ] );
         tChild->insert_basis(   7, tNodes[ 171 ] );
         tChild->insert_basis(   8, tNodes[   4 ] );
         tChild->insert_basis(   9, tNodes[   5 ] );
         tChild->insert_basis(  10, tNodes[  10 ] );
         tChild->insert_basis(  11, tNodes[  17 ] );
         tChild->insert_basis(  12, tNodes[  52 ] );
         tChild->insert_basis(  13, tNodes[ 101 ] );
         tChild->insert_basis(  14, tNodes[  13 ] );
         tChild->insert_basis(  15, tNodes[  20 ] );
         tChild->insert_basis(  16, tNodes[  55 ] );
         tChild->insert_basis(  17, tNodes[ 104 ] );
         tChild->insert_basis(  18, tNodes[  26 ] );
         tChild->insert_basis(  19, tNodes[  25 ] );
         tChild->insert_basis(  20, tNodes[  76 ] );
         tChild->insert_basis(  21, tNodes[ 125 ] );
         tChild->insert_basis(  22, tNodes[  73 ] );
         tChild->insert_basis(  23, tNodes[ 122 ] );
         tChild->insert_basis(  24, tNodes[ 151 ] );
         tChild->insert_basis(  25, tNodes[ 152 ] );
         tChild->insert_basis(  26, tNodes[ 157 ] );
         tChild->insert_basis(  27, tNodes[ 164 ] );
         tChild->insert_basis(  28, tNodes[ 160 ] );
         tChild->insert_basis(  29, tNodes[ 167 ] );
         tChild->insert_basis(  30, tNodes[ 173 ] );
         tChild->insert_basis(  31, tNodes[ 172 ] );
         tChild->insert_basis(  32, tNodes[  11 ] );
         tChild->insert_basis(  33, tNodes[  18 ] );
         tChild->insert_basis(  34, tNodes[  19 ] );
         tChild->insert_basis(  35, tNodes[  12 ] );
         tChild->insert_basis(  36, tNodes[  53 ] );
         tChild->insert_basis(  37, tNodes[  54 ] );
         tChild->insert_basis(  38, tNodes[ 103 ] );
         tChild->insert_basis(  39, tNodes[ 102 ] );
         tChild->insert_basis(  40, tNodes[  59 ] );
         tChild->insert_basis(  41, tNodes[ 108 ] );
         tChild->insert_basis(  42, tNodes[ 115 ] );
         tChild->insert_basis(  43, tNodes[  66 ] );
         tChild->insert_basis(  44, tNodes[  62 ] );
         tChild->insert_basis(  45, tNodes[  69 ] );
         tChild->insert_basis(  46, tNodes[ 118 ] );
         tChild->insert_basis(  47, tNodes[ 111 ] );
         tChild->insert_basis(  48, tNodes[  75 ] );
         tChild->insert_basis(  49, tNodes[  74 ] );
         tChild->insert_basis(  50, tNodes[ 123 ] );
         tChild->insert_basis(  51, tNodes[ 124 ] );
         tChild->insert_basis(  52, tNodes[ 158 ] );
         tChild->insert_basis(  53, tNodes[ 159 ] );
         tChild->insert_basis(  54, tNodes[ 166 ] );
         tChild->insert_basis(  55, tNodes[ 165 ] );
         tChild->insert_basis(  56, tNodes[  60 ] );
         tChild->insert_basis(  57, tNodes[  61 ] );
         tChild->insert_basis(  58, tNodes[  68 ] );
         tChild->insert_basis(  59, tNodes[  67 ] );
         tChild->insert_basis(  60, tNodes[ 109 ] );
         tChild->insert_basis(  61, tNodes[ 110 ] );
         tChild->insert_basis(  62, tNodes[ 117 ] );
         tChild->insert_basis(  63, tNodes[ 116 ] );

         // get pointer to child 2
         tChild = aAllElementsOnProc(
             mElement->get_child( 2 )->get_memory_index() );

         // init basis container for child 2
         tChild->init_basis_container();

         // link child 2 to nodes
         tChild->insert_basis(   0, tNodes[  21 ] );
         tChild->insert_basis(   1, tNodes[  24 ] );
         tChild->insert_basis(   2, tNodes[  45 ] );
         tChild->insert_basis(   3, tNodes[  42 ] );
         tChild->insert_basis(   4, tNodes[ 168 ] );
         tChild->insert_basis(   5, tNodes[ 171 ] );
         tChild->insert_basis(   6, tNodes[ 192 ] );
         tChild->insert_basis(   7, tNodes[ 189 ] );
         tChild->insert_basis(   8, tNodes[  22 ] );
         tChild->insert_basis(   9, tNodes[  23 ] );
         tChild->insert_basis(  10, tNodes[  28 ] );
         tChild->insert_basis(  11, tNodes[  35 ] );
         tChild->insert_basis(  12, tNodes[  70 ] );
         tChild->insert_basis(  13, tNodes[ 119 ] );
         tChild->insert_basis(  14, tNodes[  31 ] );
         tChild->insert_basis(  15, tNodes[  38 ] );
         tChild->insert_basis(  16, tNodes[  73 ] );
         tChild->insert_basis(  17, tNodes[ 122 ] );
         tChild->insert_basis(  18, tNodes[  44 ] );
         tChild->insert_basis(  19, tNodes[  43 ] );
         tChild->insert_basis(  20, tNodes[  94 ] );
         tChild->insert_basis(  21, tNodes[ 143 ] );
         tChild->insert_basis(  22, tNodes[  91 ] );
         tChild->insert_basis(  23, tNodes[ 140 ] );
         tChild->insert_basis(  24, tNodes[ 169 ] );
         tChild->insert_basis(  25, tNodes[ 170 ] );
         tChild->insert_basis(  26, tNodes[ 175 ] );
         tChild->insert_basis(  27, tNodes[ 182 ] );
         tChild->insert_basis(  28, tNodes[ 178 ] );
         tChild->insert_basis(  29, tNodes[ 185 ] );
         tChild->insert_basis(  30, tNodes[ 191 ] );
         tChild->insert_basis(  31, tNodes[ 190 ] );
         tChild->insert_basis(  32, tNodes[  29 ] );
         tChild->insert_basis(  33, tNodes[  36 ] );
         tChild->insert_basis(  34, tNodes[  37 ] );
         tChild->insert_basis(  35, tNodes[  30 ] );
         tChild->insert_basis(  36, tNodes[  71 ] );
         tChild->insert_basis(  37, tNodes[  72 ] );
         tChild->insert_basis(  38, tNodes[ 121 ] );
         tChild->insert_basis(  39, tNodes[ 120 ] );
         tChild->insert_basis(  40, tNodes[  77 ] );
         tChild->insert_basis(  41, tNodes[ 126 ] );
         tChild->insert_basis(  42, tNodes[ 133 ] );
         tChild->insert_basis(  43, tNodes[  84 ] );
         tChild->insert_basis(  44, tNodes[  80 ] );
         tChild->insert_basis(  45, tNodes[  87 ] );
         tChild->insert_basis(  46, tNodes[ 136 ] );
         tChild->insert_basis(  47, tNodes[ 129 ] );
         tChild->insert_basis(  48, tNodes[  93 ] );
         tChild->insert_basis(  49, tNodes[  92 ] );
         tChild->insert_basis(  50, tNodes[ 141 ] );
         tChild->insert_basis(  51, tNodes[ 142 ] );
         tChild->insert_basis(  52, tNodes[ 176 ] );
         tChild->insert_basis(  53, tNodes[ 177 ] );
         tChild->insert_basis(  54, tNodes[ 184 ] );
         tChild->insert_basis(  55, tNodes[ 183 ] );
         tChild->insert_basis(  56, tNodes[  78 ] );
         tChild->insert_basis(  57, tNodes[  79 ] );
         tChild->insert_basis(  58, tNodes[  86 ] );
         tChild->insert_basis(  59, tNodes[  85 ] );
         tChild->insert_basis(  60, tNodes[ 127 ] );
         tChild->insert_basis(  61, tNodes[ 128 ] );
         tChild->insert_basis(  62, tNodes[ 135 ] );
         tChild->insert_basis(  63, tNodes[ 134 ] );

         // get pointer to child 3
         tChild = aAllElementsOnProc(
             mElement->get_child( 3 )->get_memory_index() );

         // init basis container for child 3
         tChild->init_basis_container();

         // link child 3 to nodes
         tChild->insert_basis(   0, tNodes[  24 ] );
         tChild->insert_basis(   1, tNodes[  27 ] );
         tChild->insert_basis(   2, tNodes[  48 ] );
         tChild->insert_basis(   3, tNodes[  45 ] );
         tChild->insert_basis(   4, tNodes[ 171 ] );
         tChild->insert_basis(   5, tNodes[ 174 ] );
         tChild->insert_basis(   6, tNodes[ 195 ] );
         tChild->insert_basis(   7, tNodes[ 192 ] );
         tChild->insert_basis(   8, tNodes[  25 ] );
         tChild->insert_basis(   9, tNodes[  26 ] );
         tChild->insert_basis(  10, tNodes[  31 ] );
         tChild->insert_basis(  11, tNodes[  38 ] );
         tChild->insert_basis(  12, tNodes[  73 ] );
         tChild->insert_basis(  13, tNodes[ 122 ] );
         tChild->insert_basis(  14, tNodes[  34 ] );
         tChild->insert_basis(  15, tNodes[  41 ] );
         tChild->insert_basis(  16, tNodes[  76 ] );
         tChild->insert_basis(  17, tNodes[ 125 ] );
         tChild->insert_basis(  18, tNodes[  47 ] );
         tChild->insert_basis(  19, tNodes[  46 ] );
         tChild->insert_basis(  20, tNodes[  97 ] );
         tChild->insert_basis(  21, tNodes[ 146 ] );
         tChild->insert_basis(  22, tNodes[  94 ] );
         tChild->insert_basis(  23, tNodes[ 143 ] );
         tChild->insert_basis(  24, tNodes[ 172 ] );
         tChild->insert_basis(  25, tNodes[ 173 ] );
         tChild->insert_basis(  26, tNodes[ 178 ] );
         tChild->insert_basis(  27, tNodes[ 185 ] );
         tChild->insert_basis(  28, tNodes[ 181 ] );
         tChild->insert_basis(  29, tNodes[ 188 ] );
         tChild->insert_basis(  30, tNodes[ 194 ] );
         tChild->insert_basis(  31, tNodes[ 193 ] );
         tChild->insert_basis(  32, tNodes[  32 ] );
         tChild->insert_basis(  33, tNodes[  39 ] );
         tChild->insert_basis(  34, tNodes[  40 ] );
         tChild->insert_basis(  35, tNodes[  33 ] );
         tChild->insert_basis(  36, tNodes[  74 ] );
         tChild->insert_basis(  37, tNodes[  75 ] );
         tChild->insert_basis(  38, tNodes[ 124 ] );
         tChild->insert_basis(  39, tNodes[ 123 ] );
         tChild->insert_basis(  40, tNodes[  80 ] );
         tChild->insert_basis(  41, tNodes[ 129 ] );
         tChild->insert_basis(  42, tNodes[ 136 ] );
         tChild->insert_basis(  43, tNodes[  87 ] );
         tChild->insert_basis(  44, tNodes[  83 ] );
         tChild->insert_basis(  45, tNodes[  90 ] );
         tChild->insert_basis(  46, tNodes[ 139 ] );
         tChild->insert_basis(  47, tNodes[ 132 ] );
         tChild->insert_basis(  48, tNodes[  96 ] );
         tChild->insert_basis(  49, tNodes[  95 ] );
         tChild->insert_basis(  50, tNodes[ 144 ] );
         tChild->insert_basis(  51, tNodes[ 145 ] );
         tChild->insert_basis(  52, tNodes[ 179 ] );
         tChild->insert_basis(  53, tNodes[ 180 ] );
         tChild->insert_basis(  54, tNodes[ 187 ] );
         tChild->insert_basis(  55, tNodes[ 186 ] );
         tChild->insert_basis(  56, tNodes[  81 ] );
         tChild->insert_basis(  57, tNodes[  82 ] );
         tChild->insert_basis(  58, tNodes[  89 ] );
         tChild->insert_basis(  59, tNodes[  88 ] );
         tChild->insert_basis(  60, tNodes[ 130 ] );
         tChild->insert_basis(  61, tNodes[ 131 ] );
         tChild->insert_basis(  62, tNodes[ 138 ] );
         tChild->insert_basis(  63, tNodes[ 137 ] );

         // get pointer to child 4
         tChild = aAllElementsOnProc(
             mElement->get_child( 4 )->get_memory_index() );

         // init basis container for child 4
         tChild->init_basis_container();

         // link child 4 to nodes
         tChild->insert_basis(   0, tNodes[ 147 ] );
         tChild->insert_basis(   1, tNodes[ 150 ] );
         tChild->insert_basis(   2, tNodes[ 171 ] );
         tChild->insert_basis(   3, tNodes[ 168 ] );
         tChild->insert_basis(   4, tNodes[ 294 ] );
         tChild->insert_basis(   5, tNodes[ 297 ] );
         tChild->insert_basis(   6, tNodes[ 318 ] );
         tChild->insert_basis(   7, tNodes[ 315 ] );
         tChild->insert_basis(   8, tNodes[ 148 ] );
         tChild->insert_basis(   9, tNodes[ 149 ] );
         tChild->insert_basis(  10, tNodes[ 154 ] );
         tChild->insert_basis(  11, tNodes[ 161 ] );
         tChild->insert_basis(  12, tNodes[ 196 ] );
         tChild->insert_basis(  13, tNodes[ 245 ] );
         tChild->insert_basis(  14, tNodes[ 157 ] );
         tChild->insert_basis(  15, tNodes[ 164 ] );
         tChild->insert_basis(  16, tNodes[ 199 ] );
         tChild->insert_basis(  17, tNodes[ 248 ] );
         tChild->insert_basis(  18, tNodes[ 170 ] );
         tChild->insert_basis(  19, tNodes[ 169 ] );
         tChild->insert_basis(  20, tNodes[ 220 ] );
         tChild->insert_basis(  21, tNodes[ 269 ] );
         tChild->insert_basis(  22, tNodes[ 217 ] );
         tChild->insert_basis(  23, tNodes[ 266 ] );
         tChild->insert_basis(  24, tNodes[ 295 ] );
         tChild->insert_basis(  25, tNodes[ 296 ] );
         tChild->insert_basis(  26, tNodes[ 301 ] );
         tChild->insert_basis(  27, tNodes[ 308 ] );
         tChild->insert_basis(  28, tNodes[ 304 ] );
         tChild->insert_basis(  29, tNodes[ 311 ] );
         tChild->insert_basis(  30, tNodes[ 317 ] );
         tChild->insert_basis(  31, tNodes[ 316 ] );
         tChild->insert_basis(  32, tNodes[ 155 ] );
         tChild->insert_basis(  33, tNodes[ 162 ] );
         tChild->insert_basis(  34, tNodes[ 163 ] );
         tChild->insert_basis(  35, tNodes[ 156 ] );
         tChild->insert_basis(  36, tNodes[ 197 ] );
         tChild->insert_basis(  37, tNodes[ 198 ] );
         tChild->insert_basis(  38, tNodes[ 247 ] );
         tChild->insert_basis(  39, tNodes[ 246 ] );
         tChild->insert_basis(  40, tNodes[ 203 ] );
         tChild->insert_basis(  41, tNodes[ 252 ] );
         tChild->insert_basis(  42, tNodes[ 259 ] );
         tChild->insert_basis(  43, tNodes[ 210 ] );
         tChild->insert_basis(  44, tNodes[ 206 ] );
         tChild->insert_basis(  45, tNodes[ 213 ] );
         tChild->insert_basis(  46, tNodes[ 262 ] );
         tChild->insert_basis(  47, tNodes[ 255 ] );
         tChild->insert_basis(  48, tNodes[ 219 ] );
         tChild->insert_basis(  49, tNodes[ 218 ] );
         tChild->insert_basis(  50, tNodes[ 267 ] );
         tChild->insert_basis(  51, tNodes[ 268 ] );
         tChild->insert_basis(  52, tNodes[ 302 ] );
         tChild->insert_basis(  53, tNodes[ 303 ] );
         tChild->insert_basis(  54, tNodes[ 310 ] );
         tChild->insert_basis(  55, tNodes[ 309 ] );
         tChild->insert_basis(  56, tNodes[ 204 ] );
         tChild->insert_basis(  57, tNodes[ 205 ] );
         tChild->insert_basis(  58, tNodes[ 212 ] );
         tChild->insert_basis(  59, tNodes[ 211 ] );
         tChild->insert_basis(  60, tNodes[ 253 ] );
         tChild->insert_basis(  61, tNodes[ 254 ] );
         tChild->insert_basis(  62, tNodes[ 261 ] );
         tChild->insert_basis(  63, tNodes[ 260 ] );

         // get pointer to child 5
         tChild = aAllElementsOnProc(
             mElement->get_child( 5 )->get_memory_index() );

         // init basis container for child 5
         tChild->init_basis_container();

         // link child 5 to nodes
         tChild->insert_basis(   0, tNodes[ 150 ] );
         tChild->insert_basis(   1, tNodes[ 153 ] );
         tChild->insert_basis(   2, tNodes[ 174 ] );
         tChild->insert_basis(   3, tNodes[ 171 ] );
         tChild->insert_basis(   4, tNodes[ 297 ] );
         tChild->insert_basis(   5, tNodes[ 300 ] );
         tChild->insert_basis(   6, tNodes[ 321 ] );
         tChild->insert_basis(   7, tNodes[ 318 ] );
         tChild->insert_basis(   8, tNodes[ 151 ] );
         tChild->insert_basis(   9, tNodes[ 152 ] );
         tChild->insert_basis(  10, tNodes[ 157 ] );
         tChild->insert_basis(  11, tNodes[ 164 ] );
         tChild->insert_basis(  12, tNodes[ 199 ] );
         tChild->insert_basis(  13, tNodes[ 248 ] );
         tChild->insert_basis(  14, tNodes[ 160 ] );
         tChild->insert_basis(  15, tNodes[ 167 ] );
         tChild->insert_basis(  16, tNodes[ 202 ] );
         tChild->insert_basis(  17, tNodes[ 251 ] );
         tChild->insert_basis(  18, tNodes[ 173 ] );
         tChild->insert_basis(  19, tNodes[ 172 ] );
         tChild->insert_basis(  20, tNodes[ 223 ] );
         tChild->insert_basis(  21, tNodes[ 272 ] );
         tChild->insert_basis(  22, tNodes[ 220 ] );
         tChild->insert_basis(  23, tNodes[ 269 ] );
         tChild->insert_basis(  24, tNodes[ 298 ] );
         tChild->insert_basis(  25, tNodes[ 299 ] );
         tChild->insert_basis(  26, tNodes[ 304 ] );
         tChild->insert_basis(  27, tNodes[ 311 ] );
         tChild->insert_basis(  28, tNodes[ 307 ] );
         tChild->insert_basis(  29, tNodes[ 314 ] );
         tChild->insert_basis(  30, tNodes[ 320 ] );
         tChild->insert_basis(  31, tNodes[ 319 ] );
         tChild->insert_basis(  32, tNodes[ 158 ] );
         tChild->insert_basis(  33, tNodes[ 165 ] );
         tChild->insert_basis(  34, tNodes[ 166 ] );
         tChild->insert_basis(  35, tNodes[ 159 ] );
         tChild->insert_basis(  36, tNodes[ 200 ] );
         tChild->insert_basis(  37, tNodes[ 201 ] );
         tChild->insert_basis(  38, tNodes[ 250 ] );
         tChild->insert_basis(  39, tNodes[ 249 ] );
         tChild->insert_basis(  40, tNodes[ 206 ] );
         tChild->insert_basis(  41, tNodes[ 255 ] );
         tChild->insert_basis(  42, tNodes[ 262 ] );
         tChild->insert_basis(  43, tNodes[ 213 ] );
         tChild->insert_basis(  44, tNodes[ 209 ] );
         tChild->insert_basis(  45, tNodes[ 216 ] );
         tChild->insert_basis(  46, tNodes[ 265 ] );
         tChild->insert_basis(  47, tNodes[ 258 ] );
         tChild->insert_basis(  48, tNodes[ 222 ] );
         tChild->insert_basis(  49, tNodes[ 221 ] );
         tChild->insert_basis(  50, tNodes[ 270 ] );
         tChild->insert_basis(  51, tNodes[ 271 ] );
         tChild->insert_basis(  52, tNodes[ 305 ] );
         tChild->insert_basis(  53, tNodes[ 306 ] );
         tChild->insert_basis(  54, tNodes[ 313 ] );
         tChild->insert_basis(  55, tNodes[ 312 ] );
         tChild->insert_basis(  56, tNodes[ 207 ] );
         tChild->insert_basis(  57, tNodes[ 208 ] );
         tChild->insert_basis(  58, tNodes[ 215 ] );
         tChild->insert_basis(  59, tNodes[ 214 ] );
         tChild->insert_basis(  60, tNodes[ 256 ] );
         tChild->insert_basis(  61, tNodes[ 257 ] );
         tChild->insert_basis(  62, tNodes[ 264 ] );
         tChild->insert_basis(  63, tNodes[ 263 ] );

         // get pointer to child 6
         tChild = aAllElementsOnProc(
             mElement->get_child( 6 )->get_memory_index() );

         // init basis container for child 6
         tChild->init_basis_container();

         // link child 6 to nodes
         tChild->insert_basis(   0, tNodes[ 168 ] );
         tChild->insert_basis(   1, tNodes[ 171 ] );
         tChild->insert_basis(   2, tNodes[ 192 ] );
         tChild->insert_basis(   3, tNodes[ 189 ] );
         tChild->insert_basis(   4, tNodes[ 315 ] );
         tChild->insert_basis(   5, tNodes[ 318 ] );
         tChild->insert_basis(   6, tNodes[ 339 ] );
         tChild->insert_basis(   7, tNodes[ 336 ] );
         tChild->insert_basis(   8, tNodes[ 169 ] );
         tChild->insert_basis(   9, tNodes[ 170 ] );
         tChild->insert_basis(  10, tNodes[ 175 ] );
         tChild->insert_basis(  11, tNodes[ 182 ] );
         tChild->insert_basis(  12, tNodes[ 217 ] );
         tChild->insert_basis(  13, tNodes[ 266 ] );
         tChild->insert_basis(  14, tNodes[ 178 ] );
         tChild->insert_basis(  15, tNodes[ 185 ] );
         tChild->insert_basis(  16, tNodes[ 220 ] );
         tChild->insert_basis(  17, tNodes[ 269 ] );
         tChild->insert_basis(  18, tNodes[ 191 ] );
         tChild->insert_basis(  19, tNodes[ 190 ] );
         tChild->insert_basis(  20, tNodes[ 241 ] );
         tChild->insert_basis(  21, tNodes[ 290 ] );
         tChild->insert_basis(  22, tNodes[ 238 ] );
         tChild->insert_basis(  23, tNodes[ 287 ] );
         tChild->insert_basis(  24, tNodes[ 316 ] );
         tChild->insert_basis(  25, tNodes[ 317 ] );
         tChild->insert_basis(  26, tNodes[ 322 ] );
         tChild->insert_basis(  27, tNodes[ 329 ] );
         tChild->insert_basis(  28, tNodes[ 325 ] );
         tChild->insert_basis(  29, tNodes[ 332 ] );
         tChild->insert_basis(  30, tNodes[ 338 ] );
         tChild->insert_basis(  31, tNodes[ 337 ] );
         tChild->insert_basis(  32, tNodes[ 176 ] );
         tChild->insert_basis(  33, tNodes[ 183 ] );
         tChild->insert_basis(  34, tNodes[ 184 ] );
         tChild->insert_basis(  35, tNodes[ 177 ] );
         tChild->insert_basis(  36, tNodes[ 218 ] );
         tChild->insert_basis(  37, tNodes[ 219 ] );
         tChild->insert_basis(  38, tNodes[ 268 ] );
         tChild->insert_basis(  39, tNodes[ 267 ] );
         tChild->insert_basis(  40, tNodes[ 224 ] );
         tChild->insert_basis(  41, tNodes[ 273 ] );
         tChild->insert_basis(  42, tNodes[ 280 ] );
         tChild->insert_basis(  43, tNodes[ 231 ] );
         tChild->insert_basis(  44, tNodes[ 227 ] );
         tChild->insert_basis(  45, tNodes[ 234 ] );
         tChild->insert_basis(  46, tNodes[ 283 ] );
         tChild->insert_basis(  47, tNodes[ 276 ] );
         tChild->insert_basis(  48, tNodes[ 240 ] );
         tChild->insert_basis(  49, tNodes[ 239 ] );
         tChild->insert_basis(  50, tNodes[ 288 ] );
         tChild->insert_basis(  51, tNodes[ 289 ] );
         tChild->insert_basis(  52, tNodes[ 323 ] );
         tChild->insert_basis(  53, tNodes[ 324 ] );
         tChild->insert_basis(  54, tNodes[ 331 ] );
         tChild->insert_basis(  55, tNodes[ 330 ] );
         tChild->insert_basis(  56, tNodes[ 225 ] );
         tChild->insert_basis(  57, tNodes[ 226 ] );
         tChild->insert_basis(  58, tNodes[ 233 ] );
         tChild->insert_basis(  59, tNodes[ 232 ] );
         tChild->insert_basis(  60, tNodes[ 274 ] );
         tChild->insert_basis(  61, tNodes[ 275 ] );
         tChild->insert_basis(  62, tNodes[ 282 ] );
         tChild->insert_basis(  63, tNodes[ 281 ] );

         // get pointer to child 7
         tChild = aAllElementsOnProc(
             mElement->get_child( 7 )->get_memory_index() );

         // init basis container for child 7
         tChild->init_basis_container();

         // link child 7 to nodes
         tChild->insert_basis(   0, tNodes[ 171 ] );
         tChild->insert_basis(   1, tNodes[ 174 ] );
         tChild->insert_basis(   2, tNodes[ 195 ] );
         tChild->insert_basis(   3, tNodes[ 192 ] );
         tChild->insert_basis(   4, tNodes[ 318 ] );
         tChild->insert_basis(   5, tNodes[ 321 ] );
         tChild->insert_basis(   6, tNodes[ 342 ] );
         tChild->insert_basis(   7, tNodes[ 339 ] );
         tChild->insert_basis(   8, tNodes[ 172 ] );
         tChild->insert_basis(   9, tNodes[ 173 ] );
         tChild->insert_basis(  10, tNodes[ 178 ] );
         tChild->insert_basis(  11, tNodes[ 185 ] );
         tChild->insert_basis(  12, tNodes[ 220 ] );
         tChild->insert_basis(  13, tNodes[ 269 ] );
         tChild->insert_basis(  14, tNodes[ 181 ] );
         tChild->insert_basis(  15, tNodes[ 188 ] );
         tChild->insert_basis(  16, tNodes[ 223 ] );
         tChild->insert_basis(  17, tNodes[ 272 ] );
         tChild->insert_basis(  18, tNodes[ 194 ] );
         tChild->insert_basis(  19, tNodes[ 193 ] );
         tChild->insert_basis(  20, tNodes[ 244 ] );
         tChild->insert_basis(  21, tNodes[ 293 ] );
         tChild->insert_basis(  22, tNodes[ 241 ] );
         tChild->insert_basis(  23, tNodes[ 290 ] );
         tChild->insert_basis(  24, tNodes[ 319 ] );
         tChild->insert_basis(  25, tNodes[ 320 ] );
         tChild->insert_basis(  26, tNodes[ 325 ] );
         tChild->insert_basis(  27, tNodes[ 332 ] );
         tChild->insert_basis(  28, tNodes[ 328 ] );
         tChild->insert_basis(  29, tNodes[ 335 ] );
         tChild->insert_basis(  30, tNodes[ 341 ] );
         tChild->insert_basis(  31, tNodes[ 340 ] );
         tChild->insert_basis(  32, tNodes[ 179 ] );
         tChild->insert_basis(  33, tNodes[ 186 ] );
         tChild->insert_basis(  34, tNodes[ 187 ] );
         tChild->insert_basis(  35, tNodes[ 180 ] );
         tChild->insert_basis(  36, tNodes[ 221 ] );
         tChild->insert_basis(  37, tNodes[ 222 ] );
         tChild->insert_basis(  38, tNodes[ 271 ] );
         tChild->insert_basis(  39, tNodes[ 270 ] );
         tChild->insert_basis(  40, tNodes[ 227 ] );
         tChild->insert_basis(  41, tNodes[ 276 ] );
         tChild->insert_basis(  42, tNodes[ 283 ] );
         tChild->insert_basis(  43, tNodes[ 234 ] );
         tChild->insert_basis(  44, tNodes[ 230 ] );
         tChild->insert_basis(  45, tNodes[ 237 ] );
         tChild->insert_basis(  46, tNodes[ 286 ] );
         tChild->insert_basis(  47, tNodes[ 279 ] );
         tChild->insert_basis(  48, tNodes[ 243 ] );
         tChild->insert_basis(  49, tNodes[ 242 ] );
         tChild->insert_basis(  50, tNodes[ 291 ] );
         tChild->insert_basis(  51, tNodes[ 292 ] );
         tChild->insert_basis(  52, tNodes[ 326 ] );
         tChild->insert_basis(  53, tNodes[ 327 ] );
         tChild->insert_basis(  54, tNodes[ 334 ] );
         tChild->insert_basis(  55, tNodes[ 333 ] );
         tChild->insert_basis(  56, tNodes[ 228 ] );
         tChild->insert_basis(  57, tNodes[ 229 ] );
         tChild->insert_basis(  58, tNodes[ 236 ] );
         tChild->insert_basis(  59, tNodes[ 235 ] );
         tChild->insert_basis(  60, tNodes[ 277 ] );
         tChild->insert_basis(  61, tNodes[ 278 ] );
         tChild->insert_basis(  62, tNodes[ 285 ] );
         tChild->insert_basis(  63, tNodes[ 284 ] );

        // set flag that this element has been processed
        this->set_children_basis_flag();

        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX64_HPP_ */

