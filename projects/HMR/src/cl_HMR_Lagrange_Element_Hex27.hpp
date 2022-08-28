/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Element_Hex27.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX27_HPP_
#define SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX27_HPP_

#include "cl_HMR_Lagrange_Element.hpp"
#include "cl_MTK_Cell_Info_Hex27.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------
        template<>
        inline
        void
        Lagrange_Element< 3, 27 >::set_cell_info()
        {
            std::shared_ptr<moris::mtk::Cell_Info> tCellInfo = std::make_shared<moris::mtk::Cell_Info_Hex27 >();

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
        std::string Lagrange_Element< 3, 27 >::get_gmsh_string()
        {
            // gmsh type - number of tags - physical tag - geometry tag
            std::string aString = "12 2 0 1";

            aString += " " + std::to_string(
                             this->get_basis(  0 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  1 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  2 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  3 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  4 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  5 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  6 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  7 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  8 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 11 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 12 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis(  9 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 13 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 10 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 14 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 15 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 16 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 19 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 17 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 18 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 21 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 25 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 23 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 24 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 26 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 22 )->get_memory_index() + 1 );

            aString += " " + std::to_string(
                             this->get_basis( 20 )->get_memory_index() + 1 );

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
        Lagrange_Element< 3, 27 >::get_vtk_type()
        {
            return 29;
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
        void Lagrange_Element< 3, 27 >::get_basis_indices_for_vtk( Matrix< DDLUMat > & aBasis )
        {
            // assemble nodes in correct order
           aBasis(  0 ) =  mNodes[  0 ]->get_memory_index();
           aBasis(  1 ) =  mNodes[  1 ]->get_memory_index();
           aBasis(  2 ) =  mNodes[  2 ]->get_memory_index();
           aBasis(  3 ) =  mNodes[  3 ]->get_memory_index();
           aBasis(  4 ) =  mNodes[  4 ]->get_memory_index();
           aBasis(  5 ) =  mNodes[  5 ]->get_memory_index();
           aBasis(  6 ) =  mNodes[  6 ]->get_memory_index();
           aBasis(  7 ) =  mNodes[  7 ]->get_memory_index();
           aBasis(  8 ) =  mNodes[  8 ]->get_memory_index();
           aBasis(  9 ) =  mNodes[  9 ]->get_memory_index();
           aBasis( 10 ) =  mNodes[ 10 ]->get_memory_index();
           aBasis( 11 ) =  mNodes[ 11 ]->get_memory_index();
           aBasis( 12 ) =  mNodes[ 16 ]->get_memory_index();
           aBasis( 13 ) =  mNodes[ 17 ]->get_memory_index();
           aBasis( 14 ) =  mNodes[ 18 ]->get_memory_index();
           aBasis( 15 ) =  mNodes[ 19 ]->get_memory_index();
           aBasis( 16 ) =  mNodes[ 12 ]->get_memory_index();
           aBasis( 17 ) =  mNodes[ 13 ]->get_memory_index();
           aBasis( 18 ) =  mNodes[ 14 ]->get_memory_index();
           aBasis( 19 ) =  mNodes[ 15 ]->get_memory_index();
           aBasis( 20 ) =  mNodes[ 23 ]->get_memory_index();
           aBasis( 21 ) =  mNodes[ 24 ]->get_memory_index();
           aBasis( 22 ) =  mNodes[ 25 ]->get_memory_index();
           aBasis( 23 ) =  mNodes[ 26 ]->get_memory_index();
           aBasis( 24 ) =  mNodes[ 21 ]->get_memory_index();
           aBasis( 25 ) =  mNodes[ 22 ]->get_memory_index();
           aBasis( 26 ) =  mNodes[ 20 ]->get_memory_index();
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
        void Lagrange_Element< 3, 27 >::get_ijk_of_basis( const uint   & aBasisNumber,
                                                                 luint * aIJK )
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
            aIJK[ 0 ] += 2*tElIJK[ 0 ];
            aIJK[ 1 ] += 2*tElIJK[ 1 ];
            aIJK[ 2 ] += 2*tElIJK[ 2 ];
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
        void Lagrange_Element< 3, 27 >::create_basis_on_level_zero( moris::Cell< Element * > & aAllElementsOnProc,
                                                                    luint                    & aBasisCounter )
        {
             // initialize container for nodes
             this->init_basis_container();

             // get pointer to neighbor 4
             Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

             // test if neighbor 4 exists
             if ( tNeighbor != NULL )
             {
                 // copy nodes from this neighbor
                 mNodes[  0 ] = tNeighbor->get_basis(  4 );
                 mNodes[  1 ] = tNeighbor->get_basis(  5 );
                 mNodes[  2 ] = tNeighbor->get_basis(  6 );
                 mNodes[  3 ] = tNeighbor->get_basis(  7 );
                 mNodes[  8 ] = tNeighbor->get_basis( 16 );
                 mNodes[  9 ] = tNeighbor->get_basis( 17 );
                 mNodes[ 10 ] = tNeighbor->get_basis( 18 );
                 mNodes[ 11 ] = tNeighbor->get_basis( 19 );
                 mNodes[ 21 ] = tNeighbor->get_basis( 22 );
             }

             // get pointer to neighbor 0
             tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

             // test if neighbor 0 exists
             if ( tNeighbor != NULL )
             {
                 // copy nodes from this neighbor
                 mNodes[  0 ] = tNeighbor->get_basis(  3 );
                 mNodes[  1 ] = tNeighbor->get_basis(  2 );
                 mNodes[  4 ] = tNeighbor->get_basis(  7 );
                 mNodes[  5 ] = tNeighbor->get_basis(  6 );
                 mNodes[  8 ] = tNeighbor->get_basis( 10 );
                 mNodes[ 12 ] = tNeighbor->get_basis( 15 );
                 mNodes[ 13 ] = tNeighbor->get_basis( 14 );
                 mNodes[ 16 ] = tNeighbor->get_basis( 18 );
                 mNodes[ 25 ] = tNeighbor->get_basis( 26 );
             }

             // get pointer to neighbor 3
             tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

             // test if neighbor 3 exists
             if ( tNeighbor != NULL )
             {
                 // copy nodes from this neighbor
                 mNodes[  0 ] = tNeighbor->get_basis(  1 );
                 mNodes[  3 ] = tNeighbor->get_basis(  2 );
                 mNodes[  4 ] = tNeighbor->get_basis(  5 );
                 mNodes[  7 ] = tNeighbor->get_basis(  6 );
                 mNodes[ 11 ] = tNeighbor->get_basis(  9 );
                 mNodes[ 12 ] = tNeighbor->get_basis( 13 );
                 mNodes[ 15 ] = tNeighbor->get_basis( 14 );
                 mNodes[ 19 ] = tNeighbor->get_basis( 17 );
                 mNodes[ 23 ] = tNeighbor->get_basis( 24 );
             }

             // loop over all nodes
             for( uint k=0; k<27; ++k )
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
        void Lagrange_Element< 3, 27 >::create_basis_for_children( moris::Cell< Element * > & aAllElementsOnProc,
                                                                   luint                    & aBasisCounter )
        {
            // create temporary array containing all nodes
            Basis* tNodes[ 125 ] = { nullptr };

            // copy my own nodes into this array
            tNodes[   0 ] = mNodes[   0 ];
            tNodes[   2 ] = mNodes[   8 ];
            tNodes[   4 ] = mNodes[   1 ];
            tNodes[  10 ] = mNodes[  11 ];
            tNodes[  12 ] = mNodes[  21 ];
            tNodes[  14 ] = mNodes[   9 ];
            tNodes[  20 ] = mNodes[   3 ];
            tNodes[  22 ] = mNodes[  10 ];
            tNodes[  24 ] = mNodes[   2 ];
            tNodes[  50 ] = mNodes[  12 ];
            tNodes[  52 ] = mNodes[  25 ];
            tNodes[  54 ] = mNodes[  13 ];
            tNodes[  60 ] = mNodes[  23 ];
            tNodes[  62 ] = mNodes[  20 ];
            tNodes[  64 ] = mNodes[  24 ];
            tNodes[  70 ] = mNodes[  15 ];
            tNodes[  72 ] = mNodes[  26 ];
            tNodes[  74 ] = mNodes[  14 ];
            tNodes[ 100 ] = mNodes[   4 ];
            tNodes[ 102 ] = mNodes[  16 ];
            tNodes[ 104 ] = mNodes[   5 ];
            tNodes[ 110 ] = mNodes[  19 ];
            tNodes[ 112 ] = mNodes[  22 ];
            tNodes[ 114 ] = mNodes[  17 ];
            tNodes[ 120 ] = mNodes[   7 ];
            tNodes[ 122 ] = mNodes[  18 ];
            tNodes[ 124 ] = mNodes[   6 ];

            // get pointer to neighbor
            Element* tNeighbor = this->get_neighbor( aAllElementsOnProc, 0 );

            // test if neighbor 0 exists
            if ( tNeighbor != NULL )
            {
                // test if nodes on face 0 exist
                if ( tNeighbor->children_have_basis() )
                {
                    // get pointer to background element
                    Background_Element_Base* tBackNeighbor = tNeighbor->get_background_element();

                    // get pointer to child 2
                    Element* tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[   1 ] = tChild->get_basis(  10 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[   3 ] = tChild->get_basis(  10 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  25 ] = tChild->get_basis(  15 );
                    tNodes[  26 ] = tChild->get_basis(  26 );
                    tNodes[  27 ] = tChild->get_basis(  14 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  28 ] = tChild->get_basis(  26 );
                    tNodes[  29 ] = tChild->get_basis(  14 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  51 ] = tChild->get_basis(  18 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  53 ] = tChild->get_basis(  18 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  75 ] = tChild->get_basis(  15 );
                    tNodes[  76 ] = tChild->get_basis(  26 );
                    tNodes[  77 ] = tChild->get_basis(  14 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  78 ] = tChild->get_basis(  26 );
                    tNodes[  79 ] = tChild->get_basis(  14 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[ 101 ] = tChild->get_basis(  18 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[ 103 ] = tChild->get_basis(  18 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 1 );

            // test if neighbor 1 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   9 ] = tChild->get_basis(  11 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  19 ] = tChild->get_basis(  11 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  29 ] = tChild->get_basis(  12 );
                    tNodes[  34 ] = tChild->get_basis(  23 );
                    tNodes[  39 ] = tChild->get_basis(  15 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  44 ] = tChild->get_basis(  23 );
                    tNodes[  49 ] = tChild->get_basis(  15 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  59 ] = tChild->get_basis(  19 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[  69 ] = tChild->get_basis(  19 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[  79 ] = tChild->get_basis(  12 );
                    tNodes[  84 ] = tChild->get_basis(  23 );
                    tNodes[  89 ] = tChild->get_basis(  15 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  94 ] = tChild->get_basis(  23 );
                    tNodes[  99 ] = tChild->get_basis(  15 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[ 109 ] = tChild->get_basis(  19 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[ 119 ] = tChild->get_basis(  19 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 2 );

            // test if neighbor 2 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  21 ] = tChild->get_basis(   8 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  23 ] = tChild->get_basis(   8 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  45 ] = tChild->get_basis(  12 );
                    tNodes[  46 ] = tChild->get_basis(  25 );
                    tNodes[  47 ] = tChild->get_basis(  13 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  48 ] = tChild->get_basis(  25 );
                    tNodes[  49 ] = tChild->get_basis(  13 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[  71 ] = tChild->get_basis(  16 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  73 ] = tChild->get_basis(  16 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[  95 ] = tChild->get_basis(  12 );
                    tNodes[  96 ] = tChild->get_basis(  25 );
                    tNodes[  97 ] = tChild->get_basis(  13 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[  98 ] = tChild->get_basis(  25 );
                    tNodes[  99 ] = tChild->get_basis(  13 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[ 121 ] = tChild->get_basis(  16 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[ 123 ] = tChild->get_basis(  16 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 3 );

            // test if neighbor 3 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   5 ] = tChild->get_basis(   9 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  15 ] = tChild->get_basis(   9 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  25 ] = tChild->get_basis(  13 );
                    tNodes[  30 ] = tChild->get_basis(  24 );
                    tNodes[  35 ] = tChild->get_basis(  14 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  40 ] = tChild->get_basis(  24 );
                    tNodes[  45 ] = tChild->get_basis(  14 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[  55 ] = tChild->get_basis(  17 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[  65 ] = tChild->get_basis(  17 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[  75 ] = tChild->get_basis(  13 );
                    tNodes[  80 ] = tChild->get_basis(  24 );
                    tNodes[  85 ] = tChild->get_basis(  14 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  90 ] = tChild->get_basis(  24 );
                    tNodes[  95 ] = tChild->get_basis(  14 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[ 105 ] = tChild->get_basis(  17 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[ 115 ] = tChild->get_basis(  17 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 4 );

            // test if neighbor 4 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   1 ] = tChild->get_basis(  16 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[   3 ] = tChild->get_basis(  16 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[   5 ] = tChild->get_basis(  19 );
                    tNodes[   6 ] = tChild->get_basis(  22 );
                    tNodes[   7 ] = tChild->get_basis(  17 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[   8 ] = tChild->get_basis(  22 );
                    tNodes[   9 ] = tChild->get_basis(  17 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[  11 ] = tChild->get_basis(  18 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[  13 ] = tChild->get_basis(  18 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  15 ] = tChild->get_basis(  19 );
                    tNodes[  16 ] = tChild->get_basis(  22 );
                    tNodes[  17 ] = tChild->get_basis(  17 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  18 ] = tChild->get_basis(  22 );
                    tNodes[  19 ] = tChild->get_basis(  17 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  21 ] = tChild->get_basis(  18 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  23 ] = tChild->get_basis(  18 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 5 );

            // test if neighbor 5 exists
            if ( tNeighbor != NULL )
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

                    tNodes[ 101 ] = tChild->get_basis(   8 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[ 103 ] = tChild->get_basis(   8 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[ 105 ] = tChild->get_basis(  11 );
                    tNodes[ 106 ] = tChild->get_basis(  21 );
                    tNodes[ 107 ] = tChild->get_basis(   9 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[ 108 ] = tChild->get_basis(  21 );
                    tNodes[ 109 ] = tChild->get_basis(   9 );

                    // get pointer to child 0
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 0 )->get_memory_index() );

                    tNodes[ 111 ] = tChild->get_basis(  10 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[ 113 ] = tChild->get_basis(  10 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[ 115 ] = tChild->get_basis(  11 );
                    tNodes[ 116 ] = tChild->get_basis(  21 );
                    tNodes[ 117 ] = tChild->get_basis(   9 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[ 118 ] = tChild->get_basis(  21 );
                    tNodes[ 119 ] = tChild->get_basis(   9 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[ 121 ] = tChild->get_basis(  10 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[ 123 ] = tChild->get_basis(  10 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 6 );

            // test if neighbor 6 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   1 ] = tChild->get_basis(  18 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[   3 ] = tChild->get_basis(  18 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 7 );

            // test if neighbor 7 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   9 ] = tChild->get_basis(  19 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  19 ] = tChild->get_basis(  19 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 8 );

            // test if neighbor 8 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  21 ] = tChild->get_basis(  16 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[  23 ] = tChild->get_basis(  16 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 9 );

            // test if neighbor 9 exists
            if ( tNeighbor != NULL )
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

                    tNodes[   5 ] = tChild->get_basis(  17 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  15 ] = tChild->get_basis(  17 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 10 );

            // test if neighbor 10 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  25 ] = tChild->get_basis(  14 );

                    // get pointer to child 7
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 7 )->get_memory_index() );

                    tNodes[  75 ] = tChild->get_basis(  14 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 11 );

            // test if neighbor 11 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  29 ] = tChild->get_basis(  15 );

                    // get pointer to child 6
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 6 )->get_memory_index() );

                    tNodes[  79 ] = tChild->get_basis(  15 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 12 );

            // test if neighbor 12 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  49 ] = tChild->get_basis(  12 );

                    // get pointer to child 4
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 4 )->get_memory_index() );

                    tNodes[  99 ] = tChild->get_basis(  12 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 13 );

            // test if neighbor 13 exists
            if ( tNeighbor != NULL )
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

                    tNodes[  45 ] = tChild->get_basis(  13 );

                    // get pointer to child 5
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 5 )->get_memory_index() );

                    tNodes[  95 ] = tChild->get_basis(  13 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 14 );

            // test if neighbor 14 exists
            if ( tNeighbor != NULL )
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

                    tNodes[ 101 ] = tChild->get_basis(  10 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[ 103 ] = tChild->get_basis(  10 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 15 );

            // test if neighbor 15 exists
            if ( tNeighbor != NULL )
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

                    tNodes[ 109 ] = tChild->get_basis(  11 );

                    // get pointer to child 2
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 2 )->get_memory_index() );

                    tNodes[ 119 ] = tChild->get_basis(  11 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 16 );

            // test if neighbor 16 exists
            if ( tNeighbor != NULL )
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

                    tNodes[ 121 ] = tChild->get_basis(   8 );

                    // get pointer to child 1
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 1 )->get_memory_index() );

                    tNodes[ 123 ] = tChild->get_basis(   8 );
                }
            }

            // get pointer to neighbor
            tNeighbor = this->get_neighbor( aAllElementsOnProc, 17 );

            // test if neighbor 17 exists
            if ( tNeighbor != NULL )
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

                    tNodes[ 105 ] = tChild->get_basis(   9 );

                    // get pointer to child 3
                    tChild = aAllElementsOnProc(
                        tBackNeighbor->get_child( 3 )->get_memory_index() );

                    tNodes[ 115 ] = tChild->get_basis(   9 );
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
            tAnchor[ 0 ] = 4 * tElIJK[ 0 ];
            tAnchor[ 1 ] = 4 * tElIJK[ 1 ];
            tAnchor[ 2 ] = 4 * tElIJK[ 2 ];

            // array containing node position;
            luint tIJK[ 3 ] = { 0, 0, 0 };

            // test if node 1 exists
            if ( tNodes[ 1 ] == NULL )
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
            if ( tNodes[ 3 ] == NULL )
            {
                 // calculate position of node 3
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 3
                 tNodes[ 3 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 5 exists
            if ( tNodes[ 5 ] == NULL )
            {
                 // calculate position of node 5
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 5
                 tNodes[ 5 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 6 exists
            if ( tNodes[ 6 ] == NULL )
            {
                 // calculate position of node 6
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 6
                 tNodes[ 6 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 7 exists
            if ( tNodes[ 7 ] == NULL )
            {
                 // calculate position of node 7
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 7
                 tNodes[ 7 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 8 exists
            if ( tNodes[ 8 ] == NULL )
            {
                 // calculate position of node 8
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 8
                 tNodes[ 8 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 9 exists
            if ( tNodes[ 9 ] == NULL )
            {
                 // calculate position of node 9
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 9
                 tNodes[ 9 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 11 exists
            if ( tNodes[ 11 ] == NULL )
            {
                 // calculate position of node 11
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 11
                 tNodes[ 11 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 13 exists
            if ( tNodes[ 13 ] == NULL )
            {
                 // calculate position of node 13
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 13
                 tNodes[ 13 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 15 exists
            if ( tNodes[ 15 ] == NULL )
            {
                 // calculate position of node 15
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 15
                 tNodes[ 15 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 16 exists
            if ( tNodes[ 16 ] == NULL )
            {
                 // calculate position of node 16
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 16
                 tNodes[ 16 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 17 exists
            if ( tNodes[ 17 ] == NULL )
            {
                 // calculate position of node 17
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 17
                 tNodes[ 17 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 18 exists
            if ( tNodes[ 18 ] == NULL )
            {
                 // calculate position of node 18
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 18
                 tNodes[ 18 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 19 exists
            if ( tNodes[ 19 ] == NULL )
            {
                 // calculate position of node 19
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 19
                 tNodes[ 19 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 21 exists
            if ( tNodes[ 21 ] == NULL )
            {
                 // calculate position of node 21
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 21
                 tNodes[ 21 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 23 exists
            if ( tNodes[ 23 ] == NULL )
            {
                 // calculate position of node 23
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ];

                 // create node 23
                 tNodes[ 23 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 25 exists
            if ( tNodes[ 25 ] == NULL )
            {
                 // calculate position of node 25
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 25
                 tNodes[ 25 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 26 exists
            if ( tNodes[ 26 ] == NULL )
            {
                 // calculate position of node 26
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 26
                 tNodes[ 26 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 27 exists
            if ( tNodes[ 27 ] == NULL )
            {
                 // calculate position of node 27
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 27
                 tNodes[ 27 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 28 exists
            if ( tNodes[ 28 ] == NULL )
            {
                 // calculate position of node 28
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 28
                 tNodes[ 28 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 29 exists
            if ( tNodes[ 29 ] == NULL )
            {
                 // calculate position of node 29
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 29
                 tNodes[ 29 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 30 exists
            if ( tNodes[ 30 ] == NULL )
            {
                 // calculate position of node 30
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 30
                 tNodes[ 30 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 31
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 31
             tNodes[ 31 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 32
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 32
             tNodes[ 32 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 33
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 33
             tNodes[ 33 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 34 exists
            if ( tNodes[ 34 ] == NULL )
            {
                 // calculate position of node 34
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 34
                 tNodes[ 34 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 35 exists
            if ( tNodes[ 35 ] == NULL )
            {
                 // calculate position of node 35
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 35
                 tNodes[ 35 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 36
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 36
             tNodes[ 36 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 37
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 37
             tNodes[ 37 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 38
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 38
             tNodes[ 38 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 39 exists
            if ( tNodes[ 39 ] == NULL )
            {
                 // calculate position of node 39
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 39
                 tNodes[ 39 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 40 exists
            if ( tNodes[ 40 ] == NULL )
            {
                 // calculate position of node 40
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 40
                 tNodes[ 40 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 41
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 41
             tNodes[ 41 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 42
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 42
             tNodes[ 42 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 43
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 1;

             // create node 43
             tNodes[ 43 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 44 exists
            if ( tNodes[ 44 ] == NULL )
            {
                 // calculate position of node 44
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 44
                 tNodes[ 44 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 45 exists
            if ( tNodes[ 45 ] == NULL )
            {
                 // calculate position of node 45
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 45
                 tNodes[ 45 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 46 exists
            if ( tNodes[ 46 ] == NULL )
            {
                 // calculate position of node 46
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 46
                 tNodes[ 46 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 47 exists
            if ( tNodes[ 47 ] == NULL )
            {
                 // calculate position of node 47
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 47
                 tNodes[ 47 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 48 exists
            if ( tNodes[ 48 ] == NULL )
            {
                 // calculate position of node 48
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 48
                 tNodes[ 48 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 49 exists
            if ( tNodes[ 49 ] == NULL )
            {
                 // calculate position of node 49
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 1;

                 // create node 49
                 tNodes[ 49 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 51 exists
            if ( tNodes[ 51 ] == NULL )
            {
                 // calculate position of node 51
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 51
                 tNodes[ 51 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 53 exists
            if ( tNodes[ 53 ] == NULL )
            {
                 // calculate position of node 53
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 53
                 tNodes[ 53 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 55 exists
            if ( tNodes[ 55 ] == NULL )
            {
                 // calculate position of node 55
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 55
                 tNodes[ 55 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 56
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 56
             tNodes[ 56 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 57
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 57
             tNodes[ 57 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 58
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 58
             tNodes[ 58 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 59 exists
            if ( tNodes[ 59 ] == NULL )
            {
                 // calculate position of node 59
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 59
                 tNodes[ 59 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 61
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 61
             tNodes[ 61 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 63
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 63
             tNodes[ 63 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 65 exists
            if ( tNodes[ 65 ] == NULL )
            {
                 // calculate position of node 65
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 65
                 tNodes[ 65 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 66
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 66
             tNodes[ 66 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 67
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 67
             tNodes[ 67 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 68
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 2;

             // create node 68
             tNodes[ 68 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 69 exists
            if ( tNodes[ 69 ] == NULL )
            {
                 // calculate position of node 69
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 69
                 tNodes[ 69 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 71 exists
            if ( tNodes[ 71 ] == NULL )
            {
                 // calculate position of node 71
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 71
                 tNodes[ 71 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 73 exists
            if ( tNodes[ 73 ] == NULL )
            {
                 // calculate position of node 73
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 2;

                 // create node 73
                 tNodes[ 73 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 75 exists
            if ( tNodes[ 75 ] == NULL )
            {
                 // calculate position of node 75
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 75
                 tNodes[ 75 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 76 exists
            if ( tNodes[ 76 ] == NULL )
            {
                 // calculate position of node 76
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 76
                 tNodes[ 76 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 77 exists
            if ( tNodes[ 77 ] == NULL )
            {
                 // calculate position of node 77
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 77
                 tNodes[ 77 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 78 exists
            if ( tNodes[ 78 ] == NULL )
            {
                 // calculate position of node 78
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 78
                 tNodes[ 78 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 79 exists
            if ( tNodes[ 79 ] == NULL )
            {
                 // calculate position of node 79
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 79
                 tNodes[ 79 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 80 exists
            if ( tNodes[ 80 ] == NULL )
            {
                 // calculate position of node 80
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 80
                 tNodes[ 80 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 81
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 81
             tNodes[ 81 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 82
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 82
             tNodes[ 82 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 83
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 1;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 83
             tNodes[ 83 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 84 exists
            if ( tNodes[ 84 ] == NULL )
            {
                 // calculate position of node 84
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 84
                 tNodes[ 84 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 85 exists
            if ( tNodes[ 85 ] == NULL )
            {
                 // calculate position of node 85
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 85
                 tNodes[ 85 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 86
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 86
             tNodes[ 86 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 87
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 87
             tNodes[ 87 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 88
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 2;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 88
             tNodes[ 88 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 89 exists
            if ( tNodes[ 89 ] == NULL )
            {
                 // calculate position of node 89
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 89
                 tNodes[ 89 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 90 exists
            if ( tNodes[ 90 ] == NULL )
            {
                 // calculate position of node 90
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 90
                 tNodes[ 90 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

             // calculate position of node 91
             tIJK[ 0 ] = tAnchor[ 0 ] + 1;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 91
             tNodes[ 91 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 92
             tIJK[ 0 ] = tAnchor[ 0 ] + 2;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 92
             tNodes[ 92 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

             // calculate position of node 93
             tIJK[ 0 ] = tAnchor[ 0 ] + 3;
             tIJK[ 1 ] = tAnchor[ 1 ] + 3;
             tIJK[ 2 ] = tAnchor[ 2 ] + 3;

             // create node 93
             tNodes[ 93 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

             // increment node counter
             ++aBasisCounter;

            // test if node 94 exists
            if ( tNodes[ 94 ] == NULL )
            {
                 // calculate position of node 94
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 94
                 tNodes[ 94 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 95 exists
            if ( tNodes[ 95 ] == NULL )
            {
                 // calculate position of node 95
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 95
                 tNodes[ 95 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 96 exists
            if ( tNodes[ 96 ] == NULL )
            {
                 // calculate position of node 96
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 96
                 tNodes[ 96 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 97 exists
            if ( tNodes[ 97 ] == NULL )
            {
                 // calculate position of node 97
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 97
                 tNodes[ 97 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 98 exists
            if ( tNodes[ 98 ] == NULL )
            {
                 // calculate position of node 98
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 98
                 tNodes[ 98 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 99 exists
            if ( tNodes[ 99 ] == NULL )
            {
                 // calculate position of node 99
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 3;

                 // create node 99
                 tNodes[ 99 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 101 exists
            if ( tNodes[ 101 ] == NULL )
            {
                 // calculate position of node 101
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 101
                 tNodes[ 101 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 103 exists
            if ( tNodes[ 103 ] == NULL )
            {
                 // calculate position of node 103
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ];
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 103
                 tNodes[ 103 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 105 exists
            if ( tNodes[ 105 ] == NULL )
            {
                 // calculate position of node 105
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 105
                 tNodes[ 105 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 106 exists
            if ( tNodes[ 106 ] == NULL )
            {
                 // calculate position of node 106
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 106
                 tNodes[ 106 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 107 exists
            if ( tNodes[ 107 ] == NULL )
            {
                 // calculate position of node 107
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 107
                 tNodes[ 107 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 108 exists
            if ( tNodes[ 108 ] == NULL )
            {
                 // calculate position of node 108
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 108
                 tNodes[ 108 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 109 exists
            if ( tNodes[ 109 ] == NULL )
            {
                 // calculate position of node 109
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 1;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 109
                 tNodes[ 109 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 111 exists
            if ( tNodes[ 111 ] == NULL )
            {
                 // calculate position of node 111
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 111
                 tNodes[ 111 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 113 exists
            if ( tNodes[ 113 ] == NULL )
            {
                 // calculate position of node 113
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 2;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 113
                 tNodes[ 113 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 115 exists
            if ( tNodes[ 115 ] == NULL )
            {
                 // calculate position of node 115
                 tIJK[ 0 ] = tAnchor[ 0 ];
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 115
                 tNodes[ 115 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 116 exists
            if ( tNodes[ 116 ] == NULL )
            {
                 // calculate position of node 116
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 116
                 tNodes[ 116 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 117 exists
            if ( tNodes[ 117 ] == NULL )
            {
                 // calculate position of node 117
                 tIJK[ 0 ] = tAnchor[ 0 ] + 2;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 117
                 tNodes[ 117 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 118 exists
            if ( tNodes[ 118 ] == NULL )
            {
                 // calculate position of node 118
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 118
                 tNodes[ 118 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 119 exists
            if ( tNodes[ 119 ] == NULL )
            {
                 // calculate position of node 119
                 tIJK[ 0 ] = tAnchor[ 0 ] + 4;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 3;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 119
                 tNodes[ 119 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 121 exists
            if ( tNodes[ 121 ] == NULL )
            {
                 // calculate position of node 121
                 tIJK[ 0 ] = tAnchor[ 0 ] + 1;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 121
                 tNodes[ 121 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

                 // increment node counter
                 ++aBasisCounter;
             }

            // test if node 123 exists
            if ( tNodes[ 123 ] == NULL )
            {
                 // calculate position of node 123
                 tIJK[ 0 ] = tAnchor[ 0 ] + 3;
                 tIJK[ 1 ] = tAnchor[ 1 ] + 4;
                 tIJK[ 2 ] = tAnchor[ 2 ] + 4;

                 // create node 123
                 tNodes[ 123 ] =  new Lagrange_Node< 3 >( tIJK, tLevel, tOwner );

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
             tChild->insert_basis(   1, tNodes[   2 ] );
             tChild->insert_basis(   2, tNodes[  12 ] );
             tChild->insert_basis(   3, tNodes[  10 ] );
             tChild->insert_basis(   4, tNodes[  50 ] );
             tChild->insert_basis(   5, tNodes[  52 ] );
             tChild->insert_basis(   6, tNodes[  62 ] );
             tChild->insert_basis(   7, tNodes[  60 ] );
             tChild->insert_basis(   8, tNodes[   1 ] );
             tChild->insert_basis(   9, tNodes[   7 ] );
             tChild->insert_basis(  10, tNodes[  11 ] );
             tChild->insert_basis(  11, tNodes[   5 ] );
             tChild->insert_basis(  12, tNodes[  25 ] );
             tChild->insert_basis(  13, tNodes[  27 ] );
             tChild->insert_basis(  14, tNodes[  37 ] );
             tChild->insert_basis(  15, tNodes[  35 ] );
             tChild->insert_basis(  16, tNodes[  51 ] );
             tChild->insert_basis(  17, tNodes[  57 ] );
             tChild->insert_basis(  18, tNodes[  61 ] );
             tChild->insert_basis(  19, tNodes[  55 ] );
             tChild->insert_basis(  20, tNodes[  31 ] );
             tChild->insert_basis(  21, tNodes[   6 ] );
             tChild->insert_basis(  22, tNodes[  56 ] );
             tChild->insert_basis(  23, tNodes[  30 ] );
             tChild->insert_basis(  24, tNodes[  32 ] );
             tChild->insert_basis(  25, tNodes[  26 ] );
             tChild->insert_basis(  26, tNodes[  36 ] );

             // get pointer to child 1
             tChild = aAllElementsOnProc(
                 mElement->get_child( 1 )->get_memory_index() );

             // init basis container for child 1
             tChild->init_basis_container();

             // link child 1 to nodes
             tChild->insert_basis(   0, tNodes[   2 ] );
             tChild->insert_basis(   1, tNodes[   4 ] );
             tChild->insert_basis(   2, tNodes[  14 ] );
             tChild->insert_basis(   3, tNodes[  12 ] );
             tChild->insert_basis(   4, tNodes[  52 ] );
             tChild->insert_basis(   5, tNodes[  54 ] );
             tChild->insert_basis(   6, tNodes[  64 ] );
             tChild->insert_basis(   7, tNodes[  62 ] );
             tChild->insert_basis(   8, tNodes[   3 ] );
             tChild->insert_basis(   9, tNodes[   9 ] );
             tChild->insert_basis(  10, tNodes[  13 ] );
             tChild->insert_basis(  11, tNodes[   7 ] );
             tChild->insert_basis(  12, tNodes[  27 ] );
             tChild->insert_basis(  13, tNodes[  29 ] );
             tChild->insert_basis(  14, tNodes[  39 ] );
             tChild->insert_basis(  15, tNodes[  37 ] );
             tChild->insert_basis(  16, tNodes[  53 ] );
             tChild->insert_basis(  17, tNodes[  59 ] );
             tChild->insert_basis(  18, tNodes[  63 ] );
             tChild->insert_basis(  19, tNodes[  57 ] );
             tChild->insert_basis(  20, tNodes[  33 ] );
             tChild->insert_basis(  21, tNodes[   8 ] );
             tChild->insert_basis(  22, tNodes[  58 ] );
             tChild->insert_basis(  23, tNodes[  32 ] );
             tChild->insert_basis(  24, tNodes[  34 ] );
             tChild->insert_basis(  25, tNodes[  28 ] );
             tChild->insert_basis(  26, tNodes[  38 ] );

             // get pointer to child 2
             tChild = aAllElementsOnProc(
                 mElement->get_child( 2 )->get_memory_index() );

             // init basis container for child 2
             tChild->init_basis_container();

             // link child 2 to nodes
             tChild->insert_basis(   0, tNodes[  10 ] );
             tChild->insert_basis(   1, tNodes[  12 ] );
             tChild->insert_basis(   2, tNodes[  22 ] );
             tChild->insert_basis(   3, tNodes[  20 ] );
             tChild->insert_basis(   4, tNodes[  60 ] );
             tChild->insert_basis(   5, tNodes[  62 ] );
             tChild->insert_basis(   6, tNodes[  72 ] );
             tChild->insert_basis(   7, tNodes[  70 ] );
             tChild->insert_basis(   8, tNodes[  11 ] );
             tChild->insert_basis(   9, tNodes[  17 ] );
             tChild->insert_basis(  10, tNodes[  21 ] );
             tChild->insert_basis(  11, tNodes[  15 ] );
             tChild->insert_basis(  12, tNodes[  35 ] );
             tChild->insert_basis(  13, tNodes[  37 ] );
             tChild->insert_basis(  14, tNodes[  47 ] );
             tChild->insert_basis(  15, tNodes[  45 ] );
             tChild->insert_basis(  16, tNodes[  61 ] );
             tChild->insert_basis(  17, tNodes[  67 ] );
             tChild->insert_basis(  18, tNodes[  71 ] );
             tChild->insert_basis(  19, tNodes[  65 ] );
             tChild->insert_basis(  20, tNodes[  41 ] );
             tChild->insert_basis(  21, tNodes[  16 ] );
             tChild->insert_basis(  22, tNodes[  66 ] );
             tChild->insert_basis(  23, tNodes[  40 ] );
             tChild->insert_basis(  24, tNodes[  42 ] );
             tChild->insert_basis(  25, tNodes[  36 ] );
             tChild->insert_basis(  26, tNodes[  46 ] );

             // get pointer to child 3
             tChild = aAllElementsOnProc(
                 mElement->get_child( 3 )->get_memory_index() );

             // init basis container for child 3
             tChild->init_basis_container();

             // link child 3 to nodes
             tChild->insert_basis(   0, tNodes[  12 ] );
             tChild->insert_basis(   1, tNodes[  14 ] );
             tChild->insert_basis(   2, tNodes[  24 ] );
             tChild->insert_basis(   3, tNodes[  22 ] );
             tChild->insert_basis(   4, tNodes[  62 ] );
             tChild->insert_basis(   5, tNodes[  64 ] );
             tChild->insert_basis(   6, tNodes[  74 ] );
             tChild->insert_basis(   7, tNodes[  72 ] );
             tChild->insert_basis(   8, tNodes[  13 ] );
             tChild->insert_basis(   9, tNodes[  19 ] );
             tChild->insert_basis(  10, tNodes[  23 ] );
             tChild->insert_basis(  11, tNodes[  17 ] );
             tChild->insert_basis(  12, tNodes[  37 ] );
             tChild->insert_basis(  13, tNodes[  39 ] );
             tChild->insert_basis(  14, tNodes[  49 ] );
             tChild->insert_basis(  15, tNodes[  47 ] );
             tChild->insert_basis(  16, tNodes[  63 ] );
             tChild->insert_basis(  17, tNodes[  69 ] );
             tChild->insert_basis(  18, tNodes[  73 ] );
             tChild->insert_basis(  19, tNodes[  67 ] );
             tChild->insert_basis(  20, tNodes[  43 ] );
             tChild->insert_basis(  21, tNodes[  18 ] );
             tChild->insert_basis(  22, tNodes[  68 ] );
             tChild->insert_basis(  23, tNodes[  42 ] );
             tChild->insert_basis(  24, tNodes[  44 ] );
             tChild->insert_basis(  25, tNodes[  38 ] );
             tChild->insert_basis(  26, tNodes[  48 ] );

             // get pointer to child 4
             tChild = aAllElementsOnProc(
                 mElement->get_child( 4 )->get_memory_index() );

             // init basis container for child 4
             tChild->init_basis_container();

             // link child 4 to nodes
             tChild->insert_basis(   0, tNodes[  50 ] );
             tChild->insert_basis(   1, tNodes[  52 ] );
             tChild->insert_basis(   2, tNodes[  62 ] );
             tChild->insert_basis(   3, tNodes[  60 ] );
             tChild->insert_basis(   4, tNodes[ 100 ] );
             tChild->insert_basis(   5, tNodes[ 102 ] );
             tChild->insert_basis(   6, tNodes[ 112 ] );
             tChild->insert_basis(   7, tNodes[ 110 ] );
             tChild->insert_basis(   8, tNodes[  51 ] );
             tChild->insert_basis(   9, tNodes[  57 ] );
             tChild->insert_basis(  10, tNodes[  61 ] );
             tChild->insert_basis(  11, tNodes[  55 ] );
             tChild->insert_basis(  12, tNodes[  75 ] );
             tChild->insert_basis(  13, tNodes[  77 ] );
             tChild->insert_basis(  14, tNodes[  87 ] );
             tChild->insert_basis(  15, tNodes[  85 ] );
             tChild->insert_basis(  16, tNodes[ 101 ] );
             tChild->insert_basis(  17, tNodes[ 107 ] );
             tChild->insert_basis(  18, tNodes[ 111 ] );
             tChild->insert_basis(  19, tNodes[ 105 ] );
             tChild->insert_basis(  20, tNodes[  81 ] );
             tChild->insert_basis(  21, tNodes[  56 ] );
             tChild->insert_basis(  22, tNodes[ 106 ] );
             tChild->insert_basis(  23, tNodes[  80 ] );
             tChild->insert_basis(  24, tNodes[  82 ] );
             tChild->insert_basis(  25, tNodes[  76 ] );
             tChild->insert_basis(  26, tNodes[  86 ] );

             // get pointer to child 5
             tChild = aAllElementsOnProc(
                 mElement->get_child( 5 )->get_memory_index() );

             // init basis container for child 5
             tChild->init_basis_container();

             // link child 5 to nodes
             tChild->insert_basis(   0, tNodes[  52 ] );
             tChild->insert_basis(   1, tNodes[  54 ] );
             tChild->insert_basis(   2, tNodes[  64 ] );
             tChild->insert_basis(   3, tNodes[  62 ] );
             tChild->insert_basis(   4, tNodes[ 102 ] );
             tChild->insert_basis(   5, tNodes[ 104 ] );
             tChild->insert_basis(   6, tNodes[ 114 ] );
             tChild->insert_basis(   7, tNodes[ 112 ] );
             tChild->insert_basis(   8, tNodes[  53 ] );
             tChild->insert_basis(   9, tNodes[  59 ] );
             tChild->insert_basis(  10, tNodes[  63 ] );
             tChild->insert_basis(  11, tNodes[  57 ] );
             tChild->insert_basis(  12, tNodes[  77 ] );
             tChild->insert_basis(  13, tNodes[  79 ] );
             tChild->insert_basis(  14, tNodes[  89 ] );
             tChild->insert_basis(  15, tNodes[  87 ] );
             tChild->insert_basis(  16, tNodes[ 103 ] );
             tChild->insert_basis(  17, tNodes[ 109 ] );
             tChild->insert_basis(  18, tNodes[ 113 ] );
             tChild->insert_basis(  19, tNodes[ 107 ] );
             tChild->insert_basis(  20, tNodes[  83 ] );
             tChild->insert_basis(  21, tNodes[  58 ] );
             tChild->insert_basis(  22, tNodes[ 108 ] );
             tChild->insert_basis(  23, tNodes[  82 ] );
             tChild->insert_basis(  24, tNodes[  84 ] );
             tChild->insert_basis(  25, tNodes[  78 ] );
             tChild->insert_basis(  26, tNodes[  88 ] );

             // get pointer to child 6
             tChild = aAllElementsOnProc(
                 mElement->get_child( 6 )->get_memory_index() );

             // init basis container for child 6
             tChild->init_basis_container();

             // link child 6 to nodes
             tChild->insert_basis(   0, tNodes[  60 ] );
             tChild->insert_basis(   1, tNodes[  62 ] );
             tChild->insert_basis(   2, tNodes[  72 ] );
             tChild->insert_basis(   3, tNodes[  70 ] );
             tChild->insert_basis(   4, tNodes[ 110 ] );
             tChild->insert_basis(   5, tNodes[ 112 ] );
             tChild->insert_basis(   6, tNodes[ 122 ] );
             tChild->insert_basis(   7, tNodes[ 120 ] );
             tChild->insert_basis(   8, tNodes[  61 ] );
             tChild->insert_basis(   9, tNodes[  67 ] );
             tChild->insert_basis(  10, tNodes[  71 ] );
             tChild->insert_basis(  11, tNodes[  65 ] );
             tChild->insert_basis(  12, tNodes[  85 ] );
             tChild->insert_basis(  13, tNodes[  87 ] );
             tChild->insert_basis(  14, tNodes[  97 ] );
             tChild->insert_basis(  15, tNodes[  95 ] );
             tChild->insert_basis(  16, tNodes[ 111 ] );
             tChild->insert_basis(  17, tNodes[ 117 ] );
             tChild->insert_basis(  18, tNodes[ 121 ] );
             tChild->insert_basis(  19, tNodes[ 115 ] );
             tChild->insert_basis(  20, tNodes[  91 ] );
             tChild->insert_basis(  21, tNodes[  66 ] );
             tChild->insert_basis(  22, tNodes[ 116 ] );
             tChild->insert_basis(  23, tNodes[  90 ] );
             tChild->insert_basis(  24, tNodes[  92 ] );
             tChild->insert_basis(  25, tNodes[  86 ] );
             tChild->insert_basis(  26, tNodes[  96 ] );

             // get pointer to child 7
             tChild = aAllElementsOnProc(
                 mElement->get_child( 7 )->get_memory_index() );

             // init basis container for child 7
             tChild->init_basis_container();

             // link child 7 to nodes
             tChild->insert_basis(   0, tNodes[  62 ] );
             tChild->insert_basis(   1, tNodes[  64 ] );
             tChild->insert_basis(   2, tNodes[  74 ] );
             tChild->insert_basis(   3, tNodes[  72 ] );
             tChild->insert_basis(   4, tNodes[ 112 ] );
             tChild->insert_basis(   5, tNodes[ 114 ] );
             tChild->insert_basis(   6, tNodes[ 124 ] );
             tChild->insert_basis(   7, tNodes[ 122 ] );
             tChild->insert_basis(   8, tNodes[  63 ] );
             tChild->insert_basis(   9, tNodes[  69 ] );
             tChild->insert_basis(  10, tNodes[  73 ] );
             tChild->insert_basis(  11, tNodes[  67 ] );
             tChild->insert_basis(  12, tNodes[  87 ] );
             tChild->insert_basis(  13, tNodes[  89 ] );
             tChild->insert_basis(  14, tNodes[  99 ] );
             tChild->insert_basis(  15, tNodes[  97 ] );
             tChild->insert_basis(  16, tNodes[ 113 ] );
             tChild->insert_basis(  17, tNodes[ 119 ] );
             tChild->insert_basis(  18, tNodes[ 123 ] );
             tChild->insert_basis(  19, tNodes[ 117 ] );
             tChild->insert_basis(  20, tNodes[  93 ] );
             tChild->insert_basis(  21, tNodes[  68 ] );
             tChild->insert_basis(  22, tNodes[ 118 ] );
             tChild->insert_basis(  23, tNodes[  92 ] );
             tChild->insert_basis(  24, tNodes[  94 ] );
             tChild->insert_basis(  25, tNodes[  88 ] );
             tChild->insert_basis(  26, tNodes[  98 ] );

            // set flag that this element has been processed
            this->set_children_basis_flag();
        }

// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_LAGRANGE_ELEMENT_HEX27_HPP_ */

