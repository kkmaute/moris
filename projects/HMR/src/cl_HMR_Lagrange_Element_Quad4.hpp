/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Quad4.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD4_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD4_HPP_

#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------
        template<>
        inline
        void
        Lagrange_Element< 2, 4 >::set_cell_info()
        {
            std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Quad4 >();
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
        Lagrange_Element< 2, 4 >::get_gmsh_string()
        {
            // gmsh type - number of tags - physical tag - geometry tag
            std::string aString = "3 2 0 1";

            // loop over all nodes
            for( uint k=0; k<4; ++k )
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
        Lagrange_Element< 2, 4 >::get_vtk_type()
        {
            return 9;
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
        Lagrange_Element< 2, 4 >::get_basis_indices_for_vtk(
            Matrix< DDLUMat > & aBasis )
        {
            // loop over all nodes
            for( uint k=0; k<4; ++k )
            {
                aBasis( k ) = mNodes[ k ]->get_memory_index();
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
        inline
        void
        Lagrange_Element< 2, 4 >::get_ijk_of_basis(
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
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 0 ;
                    break;
                }
                case( 2 ) :
                {
                    aIJK[ 0 ] = 1 ;
                    aIJK[ 1 ] = 1 ;
                    break;
                }
                case( 3 ) :
                {
                    aIJK[ 0 ] = 0 ;
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
        Lagrange_Element< 2, 4 >::create_basis_on_level_zero(
              moris::Cell< Element * > & aAllElementsOnProc,
              luint                           & aBasisCounter )
        {
             // initialize container for nodes
             this->init_basis_container();

             // get pointer to neighbor 0
             Element* tNeighbor
                 = this->get_neighbor( aAllElementsOnProc, 0 );

             // test if neighbor 0 exists
             if ( tNeighbor != nullptr )
             {
                 // copy nodes from this neighbor
                 mNodes[  0 ] = tNeighbor->get_basis(  3 );
                 mNodes[  1 ] = tNeighbor->get_basis(  2 );
             }

             // get pointer to neighbor 3
             tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

             // test if neighbor 3 exists
             if ( tNeighbor != nullptr )
             {
                 // copy nodes from this neighbor
                 mNodes[  0 ] = tNeighbor->get_basis(  1 );
                 mNodes[  3 ] = tNeighbor->get_basis(  2 );
             }

             // loop over all nodes
             for( uint k=0; k<4; ++k )
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
        void Lagrange_Element< 2, 4 >::create_basis_for_children(
            moris::Cell< Element * > & aAllElementsOnProc,
            luint             & aBasisCounter )
        {
            // create temporary array containing all nodes
            Basis* tNodes[ 9 ] = { nullptr };

            // copy my own nodes into this array
            tNodes[   0 ] = mNodes[   0 ];
            tNodes[   2 ] = mNodes[   1 ];
            tNodes[   6 ] = mNodes[   3 ];
            tNodes[   8 ] = mNodes[   2 ];

            // get pointer to neighbor
            Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

            // test if neighbor 0 exists
            if ( tNeighbor != nullptr )
            {
                // test if nodes on edge 0 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor = tNeighbor->get_background_element();

                    // get pointer to child 2
                    Element* tChild = aAllElementsOnProc( tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[   1 ] = tChild->get_basis(   2 );
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
                    Background_Element_Base* tBackNeighbor = tNeighbor->get_background_element();

                    // get pointer to child 0
                    Element* tChild = aAllElementsOnProc( tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[   5 ] = tChild->get_basis(   3 );
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
                    Background_Element_Base* tBackNeighbor = tNeighbor->get_background_element();

                    // get pointer to child 0
                    Element* tChild = aAllElementsOnProc( tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[   7 ] = tChild->get_basis(   1 );
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
                    Background_Element_Base* tBackNeighbor = tNeighbor->get_background_element();

                    // get pointer to child 1
                    Element* tChild = aAllElementsOnProc( tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[   3 ] = tChild->get_basis(   2 );
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
            tAnchor[ 0 ] = 2 * tElIJ[ 0 ];
            tAnchor[ 1 ] = 2 * tElIJ[ 1 ];

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
                 ++aBasisCounter;
             }

            // test if node 3 exists
            if ( tNodes[ 3 ] == nullptr )
            {
                 // calculate position of node 3
                 tIJ[ 0 ] = tAnchor[ 0 ];
                 tIJ[ 1 ] = tAnchor[ 1 ] + 1;

                 // create node 3
                 tNodes[ 3 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 4
             tIJ[ 0 ] = tAnchor[ 0 ] + 1;
             tIJ[ 1 ] = tAnchor[ 1 ] + 1;

             // create node 4
             tNodes[ 4 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 5 exists
            if ( tNodes[ 5 ] == nullptr )
            {
                 // calculate position of node 5
                 tIJ[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJ[ 1 ] = tAnchor[ 1 ] + 1;

                 // create node 5
                 tNodes[ 5 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 7 exists
            if ( tNodes[ 7 ] == nullptr )
            {
                 // calculate position of node 7
                 tIJ[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJ[ 1 ] = tAnchor[ 1 ] + 2;

                 // create node 7
                 tNodes[ 7 ] =  new Lagrange_Node< 2 >( tIJ, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // pointer to child
             Element* tChild;

             // get pointer to child 0
             tChild = aAllElementsOnProc( mElement->get_child( 0 )->get_memory_index() );

             // init basis container for child 0
             tChild->init_basis_container();

             // link child 0 to nodes
             tChild->insert_basis(   0, tNodes[   0 ] );
             tChild->insert_basis(   1, tNodes[   1 ] );
             tChild->insert_basis(   2, tNodes[   4 ] );
             tChild->insert_basis(   3, tNodes[   3 ] );

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

            // set flag that this element has been processed
            this->set_children_basis_flag();
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_QUAD4_HPP_ */

