/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Quad9.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD9_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD9_HPP_

#include "cl_HMR_Lagrange_Element.hpp"

#include "cl_MTK_Cell_Info_Quad9.hpp"
namespace moris::hmr
{
// ----------------------------------------------------------------------------
    template<>
    inline
    void
    Lagrange_Element< 2, 9 >::set_cell_info()
    {
        std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Quad9 >();
        this->set_mtk_cell_info(tCellInfo);
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
    Lagrange_Element< 2, 9 >::get_gmsh_string()
    {
        // gmsh type - number of tags - physical tag - geometry tag
        std::string aString = "10 2 0 1";

        // loop over all nodes
        for( uint k=0; k<9; ++k )
        {
            // add node index to string
            aString += " " + std::to_string(
                this->get_basis_function( k )->get_memory_index() + 1 );
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
    Lagrange_Element< 2, 9 >::get_vtk_type()
    {
        return 28;
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
    Lagrange_Element< 2, 9 >::get_ijk_of_basis_function(
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
        aIJK[ 0 ] += 2*tElIJK[ 0 ];
        aIJK[ 1 ] += 2*tElIJK[ 1 ];
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
    luint Lagrange_Element< 2, 9 >::create_basis_on_level_zero(
          moris::Cell< Element * > & aAllElementsOnProc )
    {
        // Start basis counter
        luint tBasisCounter = 0;
        
        // initialize container for nodes
        this->init_basis_container();

        // get pointer to neighbor 0
        Element* tNeighbor
            = this->get_neighbor( aAllElementsOnProc, 0 );

        // test if neighbor 0 exists
        if ( tNeighbor != nullptr )
        {
            // copy nodes from this neighbor
            mNodes[  0 ] = tNeighbor->get_basis_function(  3 );
            mNodes[  1 ] = tNeighbor->get_basis_function(  2 );
            mNodes[  4 ] = tNeighbor->get_basis_function(  6 );
        }

        // get pointer to neighbor 3
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

        // test if neighbor 3 exists
        if ( tNeighbor != nullptr )
        {
            // copy nodes from this neighbor
            mNodes[  0 ] = tNeighbor->get_basis_function(  1 );
            mNodes[  3 ] = tNeighbor->get_basis_function(  2 );
            mNodes[  7 ] = tNeighbor->get_basis_function(  5 );
        }

        // loop over all nodes
        for( uint k=0; k<9; ++k )
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
    luint Lagrange_Element< 2, 9 >::create_basis_for_children(
        moris::Cell< Element * > & aAllElementsOnProc )
    {
        // Start basis counter 
        luint tBasisCounter = 0;
        
        // create temporary array containing all nodes
        Basis_Function* tNodes[ 25 ] = { nullptr };

        // copy my own nodes into this array
        tNodes[   0 ] = mNodes[   0 ];
        tNodes[   2 ] = mNodes[   4 ];
        tNodes[   4 ] = mNodes[   1 ];
        tNodes[  10 ] = mNodes[   7 ];
        tNodes[  12 ] = mNodes[   8 ];
        tNodes[  14 ] = mNodes[   5 ];
        tNodes[  20 ] = mNodes[   3 ];
        tNodes[  22 ] = mNodes[   6 ];
        tNodes[  24 ] = mNodes[   2 ];

        // get pointer to neighbor
        Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

        // test if neighbor 0 exists
        if ( tNeighbor != nullptr )
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

                tNodes[   1 ] = tChild->get_basis_function(   6 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[   3 ] = tChild->get_basis_function(   6 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

        // test if neighbor 1 exists
        if ( tNeighbor != nullptr )
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

                tNodes[   9 ] = tChild->get_basis_function(   7 );

                // get pointer to child 2
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 2 )->get_memory_index() );

                tNodes[  19 ] = tChild->get_basis_function(   7 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

        // test if neighbor 2 exists
        if ( tNeighbor != nullptr )
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

                tNodes[  21 ] = tChild->get_basis_function(   4 );

                // get pointer to child 1
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 1 )->get_memory_index() );

                tNodes[  23 ] = tChild->get_basis_function(   4 );
            }
        }

        // get pointer to neighbor
        tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

        // test if neighbor 3 exists
        if ( tNeighbor != nullptr )
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

                tNodes[   5 ] = tChild->get_basis_function(   5 );

                // get pointer to child 3
                tChild = aAllElementsOnProc(
                    tBackNeighbor->get_child( 3 )->get_memory_index() );

                tNodes[  15 ] = tChild->get_basis_function(   5 );
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
        tAnchor[ 0 ] = 4 * tElIJ[ 0 ];
        tAnchor[ 1 ] = 4 * tElIJ[ 1 ];

        // array containing node position;
        luint tIJ[ 2 ] = { 0, 0 } ;

        // test if node 1 exists
        if ( tNodes[ 1 ] == nullptr )
        {
             // calculate position of node 1
             tIJ[ 0 ] = tAnchor[ 0 ] + 1;
             tIJ[ 1 ] = tAnchor[ 1 ];

             // create node 1
             tNodes[ 1 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 3 exists
        if ( tNodes[ 3 ] == nullptr )
        {
             // calculate position of node 3
             tIJ[ 0 ] = tAnchor[ 0 ] + 3;
             tIJ[ 1 ] = tAnchor[ 1 ];

             // create node 3
             tNodes[ 3 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 5 exists
        if ( tNodes[ 5 ] == nullptr )
        {
             // calculate position of node 5
             tIJ[ 0 ] = tAnchor[ 0 ];
             tIJ[ 1 ] = tAnchor[ 1 ] + 1;

             // create node 5
             tNodes[ 5 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 6
         tIJ[ 0 ] = tAnchor[ 0 ] + 1;
         tIJ[ 1 ] = tAnchor[ 1 ] + 1;

         // create node 6
         tNodes[ 6 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 7
         tIJ[ 0 ] = tAnchor[ 0 ] + 2;
         tIJ[ 1 ] = tAnchor[ 1 ] + 1;

         // create node 7
         tNodes[ 7 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 8
         tIJ[ 0 ] = tAnchor[ 0 ] + 3;
         tIJ[ 1 ] = tAnchor[ 1 ] + 1;

         // create node 8
         tNodes[ 8 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 9 exists
        if ( tNodes[ 9 ] == nullptr )
        {
             // calculate position of node 9
             tIJ[ 0 ] = tAnchor[ 0 ] + 4;
             tIJ[ 1 ] = tAnchor[ 1 ] + 1;

             // create node 9
             tNodes[ 9 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 11
         tIJ[ 0 ] = tAnchor[ 0 ] + 1;
         tIJ[ 1 ] = tAnchor[ 1 ] + 2;

         // create node 11
         tNodes[ 11 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 13
         tIJ[ 0 ] = tAnchor[ 0 ] + 3;
         tIJ[ 1 ] = tAnchor[ 1 ] + 2;

         // create node 13
         tNodes[ 13 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 15 exists
        if ( tNodes[ 15 ] == nullptr )
        {
             // calculate position of node 15
             tIJ[ 0 ] = tAnchor[ 0 ];
             tIJ[ 1 ] = tAnchor[ 1 ] + 3;

             // create node 15
             tNodes[ 15 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

         // calculate position of node 16
         tIJ[ 0 ] = tAnchor[ 0 ] + 1;
         tIJ[ 1 ] = tAnchor[ 1 ] + 3;

         // create node 16
         tNodes[ 16 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 17
         tIJ[ 0 ] = tAnchor[ 0 ] + 2;
         tIJ[ 1 ] = tAnchor[ 1 ] + 3;

         // create node 17
         tNodes[ 17 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

         // calculate position of node 18
         tIJ[ 0 ] = tAnchor[ 0 ] + 3;
         tIJ[ 1 ] = tAnchor[ 1 ] + 3;

         // create node 18
         tNodes[ 18 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

         // increment node counter
         tBasisCounter++;

        // test if node 19 exists
        if ( tNodes[ 19 ] == nullptr )
        {
             // calculate position of node 19
             tIJ[ 0 ] = tAnchor[ 0 ] + 4;
             tIJ[ 1 ] = tAnchor[ 1 ] + 3;

             // create node 19
             tNodes[ 19 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 21 exists
        if ( tNodes[ 21 ] == nullptr )
        {
             // calculate position of node 21
             tIJ[ 0 ] = tAnchor[ 0 ] + 1;
             tIJ[ 1 ] = tAnchor[ 1 ] + 4;

             // create node 21
             tNodes[ 21 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             tBasisCounter++;
         }

        // test if node 23 exists
        if ( tNodes[ 23 ] == nullptr )
        {
             // calculate position of node 23
             tIJ[ 0 ] = tAnchor[ 0 ] + 3;
             tIJ[ 1 ] = tAnchor[ 1 ] + 4;

             // create node 23
             tNodes[ 23 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

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
         tChild->insert_basis_function(   0, tNodes[   0 ] );
         tChild->insert_basis_function(   1, tNodes[   2 ] );
         tChild->insert_basis_function(   2, tNodes[  12 ] );
         tChild->insert_basis_function(   3, tNodes[  10 ] );
         tChild->insert_basis_function(   4, tNodes[   1 ] );
         tChild->insert_basis_function(   5, tNodes[   7 ] );
         tChild->insert_basis_function(   6, tNodes[  11 ] );
         tChild->insert_basis_function(   7, tNodes[   5 ] );
         tChild->insert_basis_function(   8, tNodes[   6 ] );

         // get pointer to child 1
         tChild = aAllElementsOnProc(
             mElement->get_child( 1 )->get_memory_index() );

         // init basis container for child 1
         tChild->init_basis_container();

         // link child 1 to nodes
         tChild->insert_basis_function(   0, tNodes[   2 ] );
         tChild->insert_basis_function(   1, tNodes[   4 ] );
         tChild->insert_basis_function(   2, tNodes[  14 ] );
         tChild->insert_basis_function(   3, tNodes[  12 ] );
         tChild->insert_basis_function(   4, tNodes[   3 ] );
         tChild->insert_basis_function(   5, tNodes[   9 ] );
         tChild->insert_basis_function(   6, tNodes[  13 ] );
         tChild->insert_basis_function(   7, tNodes[   7 ] );
         tChild->insert_basis_function(   8, tNodes[   8 ] );

         // get pointer to child 2
         tChild = aAllElementsOnProc(
             mElement->get_child( 2 )->get_memory_index() );

         // init basis container for child 2
         tChild->init_basis_container();

         // link child 2 to nodes
         tChild->insert_basis_function(   0, tNodes[  10 ] );
         tChild->insert_basis_function(   1, tNodes[  12 ] );
         tChild->insert_basis_function(   2, tNodes[  22 ] );
         tChild->insert_basis_function(   3, tNodes[  20 ] );
         tChild->insert_basis_function(   4, tNodes[  11 ] );
         tChild->insert_basis_function(   5, tNodes[  17 ] );
         tChild->insert_basis_function(   6, tNodes[  21 ] );
         tChild->insert_basis_function(   7, tNodes[  15 ] );
         tChild->insert_basis_function(   8, tNodes[  16 ] );

         // get pointer to child 3
         tChild = aAllElementsOnProc(
             mElement->get_child( 3 )->get_memory_index() );

         // init basis container for child 3
         tChild->init_basis_container();

         // link child 3 to nodes
         tChild->insert_basis_function(   0, tNodes[  12 ] );
         tChild->insert_basis_function(   1, tNodes[  14 ] );
         tChild->insert_basis_function(   2, tNodes[  24 ] );
         tChild->insert_basis_function(   3, tNodes[  22 ] );
         tChild->insert_basis_function(   4, tNodes[  13 ] );
         tChild->insert_basis_function(   5, tNodes[  19 ] );
         tChild->insert_basis_function(   6, tNodes[  23 ] );
         tChild->insert_basis_function(   7, tNodes[  17 ] );
         tChild->insert_basis_function(   8, tNodes[  18 ] );

        // set flag that this element has been processed
        this->set_children_basis_flag();

        // Return basis counter
        return tBasisCounter;
    }

// ----------------------------------------------------------------------------
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD9_HPP_ */

