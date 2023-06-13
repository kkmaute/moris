/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Hex8.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX8_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX8_HPP_

#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_div.hpp"

namespace moris::hmr
{
// ----------------------------------------------------------------------------
    template<>
    inline
    void
    Lagrange_Element< 3, 8 >::set_cell_info()
    {
        std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Hex8 >();

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
    Lagrange_Element< 3, 8 >::get_gmsh_string()
    {
        // gmsh type - number of tags - physical tag - geometry tag
        std::string aString = "5 2 0 1";

        // loop over all nodes
        for( uint k=0; k<8; ++k )
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
     * VTK ID needed for VTK output
     *
     * @return uint
     */
    template<>
    inline
    uint
    Lagrange_Element< 3, 8 >::get_vtk_type()
    {
        return 12;
    }

// ----------------------------------------------------------------------------

    /**
     * node IDs needed for VTK output
     *
     * @param[out] moris::Matrix< DDLUMat >
     *
     * @return void
     *
     */
    template<>
    inline
    void
    Lagrange_Element< 3, 8 >::get_basis_indices_for_vtk(
        Matrix< DDLUMat > & aBasis )
    {
        // loop over all nodes
        for( uint k=0; k<8; ++k )
        {
            aBasis( k ) = mNodes[ k ]->get_memory_index();
        }
    }
// ----------------------------------------------------------------------------
    template<>
    inline
    void
    Lagrange_Element< 3, 8 >::get_ijk_of_basis(
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
    Lagrange_Element< 3, 8 >::create_basis_on_level_zero(
          moris::Cell< Element * > & aAllElementsOnProc,
          luint                           & aBasisCounter )
    {
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
         }

         // loop over all nodes
         for( uint k=0; k<8; ++k )
         {
             // test if node exists
             if( mNodes[ k ] == nullptr )
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
    Lagrange_Element< 3, 8 >::create_basis_for_children(
        moris::Cell< Element * > & aAllElementsOnProc,
        luint             & aBasisCounter )
    {
        // create temporary array containing all nodes
        Basis* tNodes[ 27 ] = { nullptr };

        // copy my own nodes into this array
        tNodes[   0 ] = mNodes[   0 ];
        tNodes[   2 ] = mNodes[   1 ];
        tNodes[   6 ] = mNodes[   3 ];
        tNodes[   8 ] = mNodes[   2 ];
        tNodes[  18 ] = mNodes[   4 ];
        tNodes[  20 ] = mNodes[   5 ];
        tNodes[  24 ] = mNodes[   7 ];
        tNodes[  26 ] = mNodes[   6 ];

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

                tNodes[   1 ] = tChild->get_basis(   2 );
                tNodes[   9 ] = tChild->get_basis(   7 );
                tNodes[  10 ] = tChild->get_basis(   6 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  11 ] = tChild->get_basis(   6 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[  19 ] = tChild->get_basis(   6 );
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

                tNodes[   5 ] = tChild->get_basis(   3 );
                tNodes[  11 ] = tChild->get_basis(   4 );
                tNodes[  14 ] = tChild->get_basis(   7 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  17 ] = tChild->get_basis(   7 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  23 ] = tChild->get_basis(   7 );
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

                tNodes[   7 ] = tChild->get_basis(   1 );
                tNodes[  15 ] = tChild->get_basis(   4 );
                tNodes[  16 ] = tChild->get_basis(   5 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  17 ] = tChild->get_basis(   5 );

                // get pointer to child 4
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 4 )->get_memory_index() );

                tNodes[  25 ] = tChild->get_basis(   5 );
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

                tNodes[   3 ] = tChild->get_basis(   2 );
                tNodes[   9 ] = tChild->get_basis(   5 );
                tNodes[  12 ] = tChild->get_basis(   6 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  15 ] = tChild->get_basis(   6 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[  21 ] = tChild->get_basis(   6 );
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

                tNodes[   1 ] = tChild->get_basis(   5 );
                tNodes[   3 ] = tChild->get_basis(   7 );
                tNodes[   4 ] = tChild->get_basis(   6 );

                // get pointer to child 5
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 5 )->get_memory_index() );

                tNodes[   5 ] = tChild->get_basis(   6 );

                // get pointer to child 6
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 6 )->get_memory_index() );

                tNodes[   7 ] = tChild->get_basis(   6 );
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

                tNodes[  19 ] = tChild->get_basis(   1 );
                tNodes[  21 ] = tChild->get_basis(   3 );
                tNodes[  22 ] = tChild->get_basis(   2 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  23 ] = tChild->get_basis(   2 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  25 ] = tChild->get_basis(   2 );
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

                tNodes[   1 ] = tChild->get_basis(   6 );
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

                tNodes[   5 ] = tChild->get_basis(   7 );
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

                tNodes[   7 ] = tChild->get_basis(   5 );
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

                tNodes[   3 ] = tChild->get_basis(   6 );
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

                tNodes[   9 ] = tChild->get_basis(   6 );
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

                tNodes[  11 ] = tChild->get_basis(   7 );
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

                tNodes[  17 ] = tChild->get_basis(   4 );
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

                tNodes[  15 ] = tChild->get_basis(   5 );
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

                tNodes[  19 ] = tChild->get_basis(   2 );
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

                tNodes[  23 ] = tChild->get_basis(   3 );
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

                tNodes[  25 ] = tChild->get_basis(   1 );
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

                tNodes[  21 ] = tChild->get_basis(   2 );
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
        tAnchor[ 0 ] = 2 * tElIJK[ 0 ];
        tAnchor[ 1 ] = 2 * tElIJK[ 1 ];
        tAnchor[ 2 ] = 2 * tElIJK[ 2 ];

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
             ++aBasisCounter;
         }

        // test if node 3 exists
        if ( tNodes[ 3 ] == nullptr )
        {
             // calculate position of node 3
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 3
             tNodes[ 3 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 4 exists
        if ( tNodes[ 4 ] == nullptr )
        {
             // calculate position of node 4
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 4
             tNodes[ 4 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 5 exists
        if ( tNodes[ 5 ] == nullptr )
        {
             // calculate position of node 5
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 5
             tNodes[ 5 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 7 exists
        if ( tNodes[ 7 ] == nullptr )
        {
             // calculate position of node 7
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ];

             // create node 7
             tNodes[ 7 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 9 exists
        if ( tNodes[ 9 ] == nullptr )
        {
             // calculate position of node 9
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 9
             tNodes[ 9 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 10 exists
        if ( tNodes[ 10 ] == nullptr )
        {
             // calculate position of node 10
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 10
             tNodes[ 10 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 11 exists
        if ( tNodes[ 11 ] == nullptr )
        {
             // calculate position of node 11
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 11
             tNodes[ 11 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 12 exists
        if ( tNodes[ 12 ] == nullptr )
        {
             // calculate position of node 12
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 12
             tNodes[ 12 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

         // calculate position of node 13
         tIJK[ 0 ] = tAnchor[ 0 ] + 1;
         tIJK[ 1 ] = tAnchor[ 1 ] + 1;
         tIJK[ 2 ] = tAnchor[ 2 ] + 1;

         // create node 13
         tNodes[ 13 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

         // increment node counter
         ++aBasisCounter;

        // test if node 14 exists
        if ( tNodes[ 14 ] == nullptr )
        {
             // calculate position of node 14
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 14
             tNodes[ 14 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 15 exists
        if ( tNodes[ 15 ] == nullptr )
        {
             // calculate position of node 15
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 15
             tNodes[ 15 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 16 exists
        if ( tNodes[ 16 ] == nullptr )
        {
             // calculate position of node 16
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 16
             tNodes[ 16 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 17 exists
        if ( tNodes[ 17 ] == nullptr )
        {
             // calculate position of node 17
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 17
             tNodes[ 17 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 19 exists
        if ( tNodes[ 19 ] == nullptr )
        {
             // calculate position of node 19
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ];
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 19
             tNodes[ 19 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 21 exists
        if ( tNodes[ 21 ] == nullptr )
        {
             // calculate position of node 21
             tIJK[ 0 ] = tAnchor[ 0 ];
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 21
             tNodes[ 21 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 22 exists
        if ( tNodes[ 22 ] == nullptr )
        {
             // calculate position of node 22
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 22
             tNodes[ 22 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 23 exists
        if ( tNodes[ 23 ] == nullptr )
        {
             // calculate position of node 23
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 23
             tNodes[ 23 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;
         }

        // test if node 25 exists
        if ( tNodes[ 25 ] == nullptr )
        {
             // calculate position of node 25
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 25
             tNodes[ 25 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

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
         tChild->insert_basis(   1, tNodes[   1 ] );
         tChild->insert_basis(   2, tNodes[   4 ] );
         tChild->insert_basis(   3, tNodes[   3 ] );
         tChild->insert_basis(   4, tNodes[   9 ] );
         tChild->insert_basis(   5, tNodes[  10 ] );
         tChild->insert_basis(   6, tNodes[  13 ] );
         tChild->insert_basis(   7, tNodes[  12 ] );

         // get pointer to child 1
         tChild = aAllElementsOnProc(
             mElement->get_child( 1 )->get_memory_index() );

         // init basis container for child 1
         tChild->init_basis_container();

         // link child 1 to nodes
         tChild->insert_basis(   0, tNodes[   1 ] );
         tChild->insert_basis(   1, tNodes[   2 ] );
         tChild->insert_basis(   2, tNodes[   5 ] );
         tChild->insert_basis(   3, tNodes[   4 ] );
         tChild->insert_basis(   4, tNodes[  10 ] );
         tChild->insert_basis(   5, tNodes[  11 ] );
         tChild->insert_basis(   6, tNodes[  14 ] );
         tChild->insert_basis(   7, tNodes[  13 ] );

         // get pointer to child 2
         tChild = aAllElementsOnProc(
             mElement->get_child( 2 )->get_memory_index() );

         // init basis container for child 2
         tChild->init_basis_container();

         // link child 2 to nodes
         tChild->insert_basis(   0, tNodes[   3 ] );
         tChild->insert_basis(   1, tNodes[   4 ] );
         tChild->insert_basis(   2, tNodes[   7 ] );
         tChild->insert_basis(   3, tNodes[   6 ] );
         tChild->insert_basis(   4, tNodes[  12 ] );
         tChild->insert_basis(   5, tNodes[  13 ] );
         tChild->insert_basis(   6, tNodes[  16 ] );
         tChild->insert_basis(   7, tNodes[  15 ] );

         // get pointer to child 3
         tChild = aAllElementsOnProc(
             mElement->get_child( 3 )->get_memory_index() );

         // init basis container for child 3
         tChild->init_basis_container();

         // link child 3 to nodes
         tChild->insert_basis(   0, tNodes[   4 ] );
         tChild->insert_basis(   1, tNodes[   5 ] );
         tChild->insert_basis(   2, tNodes[   8 ] );
         tChild->insert_basis(   3, tNodes[   7 ] );
         tChild->insert_basis(   4, tNodes[  13 ] );
         tChild->insert_basis(   5, tNodes[  14 ] );
         tChild->insert_basis(   6, tNodes[  17 ] );
         tChild->insert_basis(   7, tNodes[  16 ] );

         // get pointer to child 4
         tChild = aAllElementsOnProc(
             mElement->get_child( 4 )->get_memory_index() );

         // init basis container for child 4
         tChild->init_basis_container();

         // link child 4 to nodes
         tChild->insert_basis(   0, tNodes[   9 ] );
         tChild->insert_basis(   1, tNodes[  10 ] );
         tChild->insert_basis(   2, tNodes[  13 ] );
         tChild->insert_basis(   3, tNodes[  12 ] );
         tChild->insert_basis(   4, tNodes[  18 ] );
         tChild->insert_basis(   5, tNodes[  19 ] );
         tChild->insert_basis(   6, tNodes[  22 ] );
         tChild->insert_basis(   7, tNodes[  21 ] );

         // get pointer to child 5
         tChild = aAllElementsOnProc(
             mElement->get_child( 5 )->get_memory_index() );

         // init basis container for child 5
         tChild->init_basis_container();

         // link child 5 to nodes
         tChild->insert_basis(   0, tNodes[  10 ] );
         tChild->insert_basis(   1, tNodes[  11 ] );
         tChild->insert_basis(   2, tNodes[  14 ] );
         tChild->insert_basis(   3, tNodes[  13 ] );
         tChild->insert_basis(   4, tNodes[  19 ] );
         tChild->insert_basis(   5, tNodes[  20 ] );
         tChild->insert_basis(   6, tNodes[  23 ] );
         tChild->insert_basis(   7, tNodes[  22 ] );

         // get pointer to child 6
         tChild = aAllElementsOnProc(
             mElement->get_child( 6 )->get_memory_index() );

         // init basis container for child 6
         tChild->init_basis_container();

         // link child 6 to nodes
         tChild->insert_basis(   0, tNodes[  12 ] );
         tChild->insert_basis(   1, tNodes[  13 ] );
         tChild->insert_basis(   2, tNodes[  16 ] );
         tChild->insert_basis(   3, tNodes[  15 ] );
         tChild->insert_basis(   4, tNodes[  21 ] );
         tChild->insert_basis(   5, tNodes[  22 ] );
         tChild->insert_basis(   6, tNodes[  25 ] );
         tChild->insert_basis(   7, tNodes[  24 ] );

         // get pointer to child 7
         tChild = aAllElementsOnProc(
             mElement->get_child( 7 )->get_memory_index() );

         // init basis container for child 7
         tChild->init_basis_container();

         // link child 7 to nodes
         tChild->insert_basis(   0, tNodes[  13 ] );
         tChild->insert_basis(   1, tNodes[  14 ] );
         tChild->insert_basis(   2, tNodes[  17 ] );
         tChild->insert_basis(   3, tNodes[  16 ] );
         tChild->insert_basis(   4, tNodes[  22 ] );
         tChild->insert_basis(   5, tNodes[  23 ] );
         tChild->insert_basis(   6, tNodes[  26 ] );
         tChild->insert_basis(   7, tNodes[  25 ] );

        // set flag that this element has been processed
        this->set_children_basis_flag();
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX8_HPP_ */

