/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Quad16.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD16_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD16_HPP_

#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_MTK_Cell_Info_Quad16.hpp"
namespace moris
{
    namespace hmr
    {

        // ----------------------------------------------------------------------------

        template<>
        inline
        void
        Lagrange_Element< 2, 16 >::set_cell_info()
        {
            std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Quad16 >();

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
        Lagrange_Element< 2, 16 >::get_gmsh_string()
        {
            // gmsh type - number of tags - physical tag - geometry tag
            std::string aString = "36 2 0 1";

            // loop over all nodes
            for( uint k=0; k<16; ++k )
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
        Lagrange_Element< 2, 16 >::get_ijk_of_basis(
                uint aBasisNumber,
                luint      * aIJK )
        {
            // get element local coordinate
            switch ( aBasisNumber )
            {
                case 0 :
                {
                    aIJK[ 0 ] = 0 ;
                    aIJK[ 1 ] = 0 ;
                    break;
                }
                case 1 :
                {
                    aIJK[ 0 ] = 3 ;
                    aIJK[ 1 ] = 0 ;
                    break;
                }
                case 2 :
                {
                    aIJK[ 0 ] = 3 ;
                    aIJK[ 1 ] = 3 ;
                    break;
                }
                case 3 :
                {
                    aIJK[ 0 ] = 0 ;
                    aIJK[ 1 ] = 3 ;
                    break;
                }
                case 4 :
                {
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 0 ;
                    break;
                }
                case 5 :
                {
                    aIJK[ 0 ] = 2 ;
                    aIJK[ 1 ] = 0 ;
                    break;
                }
                case 6 :
                {
                    aIJK[ 0 ] = 3 ;
                    aIJK[ 1 ] = 1 ;
                    break;
                }
                case 7 :
                {
                    aIJK[ 0 ] = 3 ;
                    aIJK[ 1 ] = 2 ;
                    break;
                }
                case 8 :
                {
                    aIJK[ 0 ] = 2 ;
                    aIJK[ 1 ] = 3 ;
                    break;
                }
                case 9 :
                {
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 3 ;
                    break;
                }
                case 10 :
                {
                    aIJK[ 0 ] = 0 ;
                    aIJK[ 1 ] = 2 ;
                    break;
                }
                case 11 :
                {
                    aIJK[ 0 ] = 0 ;
                    aIJK[ 1 ] = 1 ;
                    break;
                }
                case 12 :
                {
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 1 ;
                    break;
                }
                case 13 :
                {
                    aIJK[ 0 ] = 2 ;
                    aIJK[ 1 ] = 1 ;
                    break;
                }
                case 14 :
                {
                    aIJK[ 0 ] = 2 ;
                    aIJK[ 1 ] = 2 ;
                    break;
                }
                case 15 :
                {
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 2 ;
                    break;
                }
            }

            // get position of element on background mesh
            const luint * tElIJK = mElement->get_ijk();

            // add element offset
            aIJK[ 0 ] += 3*tElIJK[ 0 ];
            aIJK[ 1 ] += 3*tElIJK[ 1 ];
        }

        // ----------------------------------------------------------------------------

        /**
         * Creates all nodes on the coarsest level.
         * Called by Lagrange mesh create_nodes_on_level_zero().
         *
         * @param[inout] aAllElementsOnProc   cell containing all Lagrange
         *                                    elements including the aura
         * @param[inout] aBasisCounter         counter to keep track of
         *                                    how many nodes were generated
         * @return void
         */
        template<>
        inline
        void
        Lagrange_Element< 2, 16 >::create_basis_on_level_zero(
                moris::Cell< Element * > & aAllElementsOnProc,
                luint                           & aBasisCounter )
        {
            // initialize container for nodes
            this->init_basis_container();

            // get pointer to neighbor 0
            Element* tNeighbor
            = this->get_neighbor( aAllElementsOnProc, 0 );

            // test if neighbor 0 exists
            if ( tNeighbor != NULL )
            {
                // copy nodes from this neighbor
                mNodes[  0 ] = tNeighbor->get_basis(  3 );
                mNodes[  1 ] = tNeighbor->get_basis(  2 );
                mNodes[  4 ] = tNeighbor->get_basis(  9 );
                mNodes[  5 ] = tNeighbor->get_basis(  8 );
            }

            // get pointer to neighbor 3
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

            // test if neighbor 3 exists
            if ( tNeighbor != NULL )
            {
                // copy nodes from this neighbor
                mNodes[  0 ] = tNeighbor->get_basis(  1 );
                mNodes[  3 ] = tNeighbor->get_basis(  2 );
                mNodes[ 10 ] = tNeighbor->get_basis(  7 );
                mNodes[ 11 ] = tNeighbor->get_basis(  6 );
            }

            // loop over all nodes
            for( uint k=0; k<16; ++k )
            {
                // test if node exists
                if( mNodes[ k ] == NULL )
                {
                    // create node
                    this->create_basis( k );

                    // increment node counter
                    ++aBasisCounter;
                }
            }
        }

        // ----------------------------------------------------------------------------

        /**
         * Creates nodes for children of refined elements.
         * Called by Lagrange mesh create_nodes_on_higher_levels().
         *
         * @param[inout] aAllElementsOnProc   cell containing all Lagrange
         *                                    elements including the aura
         * @param[inout] aNodeCounter         counter to keep track of
         *                                    how many nodes were generated
         * @return void
         */
        template<>
        inline
        void
        Lagrange_Element< 2, 16 >::create_basis_for_children(
                moris::Cell< Element * > & aAllElementsOnProc,
                luint             & aBasisCounter )
        {
            // create temporary array containing all nodes
            Basis* tNodes[ 49 ] = { nullptr };

            // copy my own nodes into this array
            tNodes[   0 ] = mNodes[   0 ];
            tNodes[   2 ] = mNodes[   4 ];
            tNodes[   4 ] = mNodes[   5 ];
            tNodes[   6 ] = mNodes[   1 ];
            tNodes[  14 ] = mNodes[  11 ];
            tNodes[  16 ] = mNodes[  12 ];
            tNodes[  18 ] = mNodes[  13 ];
            tNodes[  20 ] = mNodes[   6 ];
            tNodes[  28 ] = mNodes[  10 ];
            tNodes[  30 ] = mNodes[  15 ];
            tNodes[  32 ] = mNodes[  14 ];
            tNodes[  34 ] = mNodes[   7 ];
            tNodes[  42 ] = mNodes[   3 ];
            tNodes[  44 ] = mNodes[   9 ];
            tNodes[  46 ] = mNodes[   8 ];
            tNodes[  48 ] = mNodes[   2 ];

            // get pointer to neighbor
            Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

            // test if neighbor 0 exists
            if ( tNeighbor != NULL )
            {
                // test if nodes on edge 0 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                    // get pointer to child 2
                    Element* tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[   1 ] = tChild->get_basis(   9 );
                    tNodes[   3 ] = tChild->get_basis(   2 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[   5 ] = tChild->get_basis(   8 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

            // test if neighbor 1 exists
            if ( tNeighbor != NULL )
            {
                // test if nodes on edge 1 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                    // get pointer to child 0
                    Element* tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  13 ] = tChild->get_basis(  11 );
                    tNodes[  27 ] = tChild->get_basis(   3 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  41 ] = tChild->get_basis(  10 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

            // test if neighbor 2 exists
            if ( tNeighbor != NULL )
            {
                // test if nodes on edge 2 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                    // get pointer to child 0
                    Element* tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  43 ] = tChild->get_basis(   4 );
                    tNodes[  45 ] = tChild->get_basis(   1 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  47 ] = tChild->get_basis(   5 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

            // test if neighbor 3 exists
            if ( tNeighbor != NULL )
            {
                // test if nodes on edge 3 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor
                    = tNeighbor->get_background_element();

                    // get pointer to child 1
                    Element* tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[   7 ] = tChild->get_basis(   6 );
                    tNodes[  21 ] = tChild->get_basis(   2 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                            tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  35 ] = tChild->get_basis(   7 );
                }
            }

            // level of node
            auto tLevel = mElement->get_level() + 1;

            // owner of element
            auto tOwner = mElement->get_owner();

            // get position of element
            const luint * tElIJ = mElement->get_ijk();

            // anchor point of nodes
            luint tAnchor[ 2 ];
            tAnchor[ 0 ] = 6 * tElIJ[ 0 ];
            tAnchor[ 1 ] = 6 * tElIJ[ 1 ];

            // array containing node position;
            luint tIJ[ 2 ] = { 0, 0 } ;

            // test if node 1 exists
            if ( tNodes[ 1 ] == NULL )
            {
                // calculate position of node 1
                tIJ[ 0 ] = tAnchor[ 0 ] + 1;
                tIJ[ 1 ] = tAnchor[ 1 ];

                // create node 1
                tNodes[ 1 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 3 exists
            if ( tNodes[ 3 ] == NULL )
            {
                // calculate position of node 3
                tIJ[ 0 ] = tAnchor[ 0 ] + 3;
                tIJ[ 1 ] = tAnchor[ 1 ];

                // create node 3
                tNodes[ 3 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 5 exists
            if ( tNodes[ 5 ] == NULL )
            {
                // calculate position of node 5
                tIJ[ 0 ] = tAnchor[ 0 ] + 5;
                tIJ[ 1 ] = tAnchor[ 1 ];

                // create node 5
                tNodes[ 5 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 7 exists
            if ( tNodes[ 7 ] == NULL )
            {
                // calculate position of node 7
                tIJ[ 0 ] = tAnchor[ 0 ];
                tIJ[ 1 ] = tAnchor[ 1 ] + 1;

                // create node 7
                tNodes[ 7 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // calculate position of node 8
            tIJ[ 0 ] = tAnchor[ 0 ] + 1;
            tIJ[ 1 ] = tAnchor[ 1 ] + 1;

            // create node 8
            tNodes[ 8 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 9
            tIJ[ 0 ] = tAnchor[ 0 ] + 2;
            tIJ[ 1 ] = tAnchor[ 1 ] + 1;

            // create node 9
            tNodes[ 9 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 10
            tIJ[ 0 ] = tAnchor[ 0 ] + 3;
            tIJ[ 1 ] = tAnchor[ 1 ] + 1;

            // create node 10
            tNodes[ 10 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 11
            tIJ[ 0 ] = tAnchor[ 0 ] + 4;
            tIJ[ 1 ] = tAnchor[ 1 ] + 1;

            // create node 11
            tNodes[ 11 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 12
            tIJ[ 0 ] = tAnchor[ 0 ] + 5;
            tIJ[ 1 ] = tAnchor[ 1 ] + 1;

            // create node 12
            tNodes[ 12 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // test if node 13 exists
            if ( tNodes[ 13 ] == NULL )
            {
                // calculate position of node 13
                tIJ[ 0 ] = tAnchor[ 0 ] + 6;
                tIJ[ 1 ] = tAnchor[ 1 ] + 1;

                // create node 13
                tNodes[ 13 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // calculate position of node 15
            tIJ[ 0 ] = tAnchor[ 0 ] + 1;
            tIJ[ 1 ] = tAnchor[ 1 ] + 2;

            // create node 15
            tNodes[ 15 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 17
            tIJ[ 0 ] = tAnchor[ 0 ] + 3;
            tIJ[ 1 ] = tAnchor[ 1 ] + 2;

            // create node 17
            tNodes[ 17 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 19
            tIJ[ 0 ] = tAnchor[ 0 ] + 5;
            tIJ[ 1 ] = tAnchor[ 1 ] + 2;

            // create node 19
            tNodes[ 19 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // test if node 21 exists
            if ( tNodes[ 21 ] == NULL )
            {
                // calculate position of node 21
                tIJ[ 0 ] = tAnchor[ 0 ];
                tIJ[ 1 ] = tAnchor[ 1 ] + 3;

                // create node 21
                tNodes[ 21 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // calculate position of node 22
            tIJ[ 0 ] = tAnchor[ 0 ] + 1;
            tIJ[ 1 ] = tAnchor[ 1 ] + 3;

            // create node 22
            tNodes[ 22 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 23
            tIJ[ 0 ] = tAnchor[ 0 ] + 2;
            tIJ[ 1 ] = tAnchor[ 1 ] + 3;

            // create node 23
            tNodes[ 23 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 24
            tIJ[ 0 ] = tAnchor[ 0 ] + 3;
            tIJ[ 1 ] = tAnchor[ 1 ] + 3;

            // create node 24
            tNodes[ 24 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 25
            tIJ[ 0 ] = tAnchor[ 0 ] + 4;
            tIJ[ 1 ] = tAnchor[ 1 ] + 3;

            // create node 25
            tNodes[ 25 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 26
            tIJ[ 0 ] = tAnchor[ 0 ] + 5;
            tIJ[ 1 ] = tAnchor[ 1 ] + 3;

            // create node 26
            tNodes[ 26 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // test if node 27 exists
            if ( tNodes[ 27 ] == NULL )
            {
                // calculate position of node 27
                tIJ[ 0 ] = tAnchor[ 0 ] + 6;
                tIJ[ 1 ] = tAnchor[ 1 ] + 3;

                // create node 27
                tNodes[ 27 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // calculate position of node 29
            tIJ[ 0 ] = tAnchor[ 0 ] + 1;
            tIJ[ 1 ] = tAnchor[ 1 ] + 4;

            // create node 29
            tNodes[ 29 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 31
            tIJ[ 0 ] = tAnchor[ 0 ] + 3;
            tIJ[ 1 ] = tAnchor[ 1 ] + 4;

            // create node 31
            tNodes[ 31 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 33
            tIJ[ 0 ] = tAnchor[ 0 ] + 5;
            tIJ[ 1 ] = tAnchor[ 1 ] + 4;

            // create node 33
            tNodes[ 33 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // test if node 35 exists
            if ( tNodes[ 35 ] == NULL )
            {
                // calculate position of node 35
                tIJ[ 0 ] = tAnchor[ 0 ];
                tIJ[ 1 ] = tAnchor[ 1 ] + 5;

                // create node 35
                tNodes[ 35 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // calculate position of node 36
            tIJ[ 0 ] = tAnchor[ 0 ] + 1;
            tIJ[ 1 ] = tAnchor[ 1 ] + 5;

            // create node 36
            tNodes[ 36 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 37
            tIJ[ 0 ] = tAnchor[ 0 ] + 2;
            tIJ[ 1 ] = tAnchor[ 1 ] + 5;

            // create node 37
            tNodes[ 37 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 38
            tIJ[ 0 ] = tAnchor[ 0 ] + 3;
            tIJ[ 1 ] = tAnchor[ 1 ] + 5;

            // create node 38
            tNodes[ 38 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 39
            tIJ[ 0 ] = tAnchor[ 0 ] + 4;
            tIJ[ 1 ] = tAnchor[ 1 ] + 5;

            // create node 39
            tNodes[ 39 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // calculate position of node 40
            tIJ[ 0 ] = tAnchor[ 0 ] + 5;
            tIJ[ 1 ] = tAnchor[ 1 ] + 5;

            // create node 40
            tNodes[ 40 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

            // increment node counter
            ++aBasisCounter;

            // test if node 41 exists
            if ( tNodes[ 41 ] == NULL )
            {
                // calculate position of node 41
                tIJ[ 0 ] = tAnchor[ 0 ] + 6;
                tIJ[ 1 ] = tAnchor[ 1 ] + 5;

                // create node 41
                tNodes[ 41 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 43 exists
            if ( tNodes[ 43 ] == NULL )
            {
                // calculate position of node 43
                tIJ[ 0 ] = tAnchor[ 0 ] + 1;
                tIJ[ 1 ] = tAnchor[ 1 ] + 6;

                // create node 43
                tNodes[ 43 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 45 exists
            if ( tNodes[ 45 ] == NULL )
            {
                // calculate position of node 45
                tIJ[ 0 ] = tAnchor[ 0 ] + 3;
                tIJ[ 1 ] = tAnchor[ 1 ] + 6;

                // create node 45
                tNodes[ 45 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
            }

            // test if node 47 exists
            if ( tNodes[ 47 ] == NULL )
            {
                // calculate position of node 47
                tIJ[ 0 ] = tAnchor[ 0 ] + 5;
                tIJ[ 1 ] = tAnchor[ 1 ] + 6;

                // create node 47
                tNodes[ 47 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                // increment node counter
                ++aBasisCounter;
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
            tChild->insert_basis(   4, tNodes[   1 ] );
            tChild->insert_basis(   5, tNodes[   2 ] );
            tChild->insert_basis(   6, tNodes[  10 ] );
            tChild->insert_basis(   7, tNodes[  17 ] );
            tChild->insert_basis(   8, tNodes[  23 ] );
            tChild->insert_basis(   9, tNodes[  22 ] );
            tChild->insert_basis(  10, tNodes[  14 ] );
            tChild->insert_basis(  11, tNodes[   7 ] );
            tChild->insert_basis(  12, tNodes[   8 ] );
            tChild->insert_basis(  13, tNodes[   9 ] );
            tChild->insert_basis(  14, tNodes[  16 ] );
            tChild->insert_basis(  15, tNodes[  15 ] );

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
            tChild->insert_basis(   4, tNodes[   4 ] );
            tChild->insert_basis(   5, tNodes[   5 ] );
            tChild->insert_basis(   6, tNodes[  13 ] );
            tChild->insert_basis(   7, tNodes[  20 ] );
            tChild->insert_basis(   8, tNodes[  26 ] );
            tChild->insert_basis(   9, tNodes[  25 ] );
            tChild->insert_basis(  10, tNodes[  17 ] );
            tChild->insert_basis(  11, tNodes[  10 ] );
            tChild->insert_basis(  12, tNodes[  11 ] );
            tChild->insert_basis(  13, tNodes[  12 ] );
            tChild->insert_basis(  14, tNodes[  19 ] );
            tChild->insert_basis(  15, tNodes[  18 ] );

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
            tChild->insert_basis(   4, tNodes[  22 ] );
            tChild->insert_basis(   5, tNodes[  23 ] );
            tChild->insert_basis(   6, tNodes[  31 ] );
            tChild->insert_basis(   7, tNodes[  38 ] );
            tChild->insert_basis(   8, tNodes[  44 ] );
            tChild->insert_basis(   9, tNodes[  43 ] );
            tChild->insert_basis(  10, tNodes[  35 ] );
            tChild->insert_basis(  11, tNodes[  28 ] );
            tChild->insert_basis(  12, tNodes[  29 ] );
            tChild->insert_basis(  13, tNodes[  30 ] );
            tChild->insert_basis(  14, tNodes[  37 ] );
            tChild->insert_basis(  15, tNodes[  36 ] );

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
            tChild->insert_basis(   4, tNodes[  25 ] );
            tChild->insert_basis(   5, tNodes[  26 ] );
            tChild->insert_basis(   6, tNodes[  34 ] );
            tChild->insert_basis(   7, tNodes[  41 ] );
            tChild->insert_basis(   8, tNodes[  47 ] );
            tChild->insert_basis(   9, tNodes[  46 ] );
            tChild->insert_basis(  10, tNodes[  38 ] );
            tChild->insert_basis(  11, tNodes[  31 ] );
            tChild->insert_basis(  12, tNodes[  32 ] );
            tChild->insert_basis(  13, tNodes[  33 ] );
            tChild->insert_basis(  14, tNodes[  40 ] );
            tChild->insert_basis(  15, tNodes[  39 ] );

            // set flag that this element has been processed
            this->set_children_basis_flag();
        }

        // ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD16_HPP_ */

