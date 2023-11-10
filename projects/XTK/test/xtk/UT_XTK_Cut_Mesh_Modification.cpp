/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Mesh_Modification.cpp
 *
 */
#include "catch.hpp"
#include <iostream>
#include <chrono>
#include <thread>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Face_Registry.hpp"

// XTKL: Mesh Includes
#include "cl_MTK_Enums.hpp"
#include "fn_verify_tet_topology.hpp"

// XTKL: Linear Algebra Includes
#include "cl_XTK_Matrix_Base_Utilities.hpp"

#include "cl_Matrix.hpp"
#include "fn_bubble_sort.hpp"

namespace xtk
{

    /*
     * The face registry system is responsible for managing the creation of new faces
     * either by appending it to the end of the existing during the modification of
     * a child mesh. It returns a face index which the modification uses to ensures each new element
     * uses the same face index
     */
    TEST_CASE( "Face Registry", "[FACE_REGISTRY]" )
    {

        SECTION( "On a Tetrahedral Mesh" )
        {
            size_t tMax = std::numeric_limits< moris::moris_index >::max();
            // Starting with 2 standard tetrahedrals sharing face 3

            // Setup the problem --------------------------
            // Face to Node connectivity
            moris::Matrix< moris::IndexMat > tFaceToNodeConnectivity( { { 1, 2, 4 },
                    { 2, 3, 4 },
                    { 1, 3, 4 },
                    { 1, 2, 3 },
                    { 1, 2, 5 },
                    { 2, 3, 5 },
                    { 1, 3, 5 } } );

            // Face to Element connectivity
            moris::Matrix< moris::IndexMat > tFaceToElement( 7, 2, tMax );
            ( tFaceToElement )( 0, 0 ) = 0;
            ( tFaceToElement )( 1, 0 ) = 0;
            ( tFaceToElement )( 2, 0 ) = 0;
            ( tFaceToElement )( 3, 0 ) = 0;
            ( tFaceToElement )( 3, 1 ) = 1;
            ( tFaceToElement )( 4, 0 ) = 1;
            ( tFaceToElement )( 5, 0 ) = 1;
            ( tFaceToElement )( 6, 0 ) = 1;

            moris::Matrix< moris::IndexMat > tFaceParentIndices( { { 0, 1, 2, 3 } } );
            moris::Matrix< moris::IndexMat > tFaceParentRanks( { { 2, 2, 2, 2 } } );

            // Initialize the face registry with the given connectivity
            Face_Registry
                    tFaceRegistry( 6, 4, tFaceToNodeConnectivity, tFaceToElement, tFaceParentIndices, tFaceParentRanks );

            // Initialize Variables for testing purposes
            moris::Matrix< moris::IndexMat > tFaceIndices( 1, 4 );
            moris::Matrix< moris::IndexMat > tElementIndex( 1, 1 );
            moris::Matrix< moris::IndexMat > tExpectedParentRanks( 1, 23 );
            moris::Matrix< moris::IndexMat > tExpectedFaceToElement( 23, 2 );
            moris::Matrix< moris::IndexMat > tExpectedParentIndices( 1, 23 );
            moris::Matrix< moris::IndexMat > tSingleFaceToNodeIndices( 1, 3 );
            moris::Matrix< moris::IndexMat > tElementalFaceToNodeIndices( 4, 3 );
            moris::Matrix< moris::IndexMat > tExpectedElementalFaceIndices( 1, 4 );

            // Add Child Element 0
            tElementalFaceToNodeIndices( 0, 0 ) = 1;
            tElementalFaceToNodeIndices( 0, 1 ) = 6;
            tElementalFaceToNodeIndices( 0, 2 ) = 8;
            tElementalFaceToNodeIndices( 1, 0 ) = 1;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 7;
            tElementalFaceToNodeIndices( 2, 0 ) = 1;
            tElementalFaceToNodeIndices( 2, 1 ) = 7;
            tElementalFaceToNodeIndices( 2, 2 ) = 8;
            tElementalFaceToNodeIndices( 3, 0 ) = 6;
            tElementalFaceToNodeIndices( 3, 1 ) = 7;
            tElementalFaceToNodeIndices( 3, 2 ) = 8;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            tExpectedElementalFaceIndices.fill( tMax );

            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Replace face 0 with row 0
            // Replace face 3 with row 1
            // Replace face 2 with row 2
            // Append face row 3
            tFaceRegistry.replace_face( 0, 0, tElementalFaceToNodeIndices );
            tFaceRegistry.replace_face( 3, 1, tElementalFaceToNodeIndices );
            tFaceRegistry.replace_face( 2, 2, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 3, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 0;
            tExpectedElementalFaceIndices( 0, 1 ) = 3;
            tExpectedElementalFaceIndices( 0, 2 ) = 2;
            tExpectedElementalFaceIndices( 0, 3 ) = 7;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set Face to Element Information
            tElementIndex( 0, 0 ) = 0;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 2
            tElementalFaceToNodeIndices( 0, 0 ) = 4;
            tElementalFaceToNodeIndices( 0, 1 ) = 6;
            tElementalFaceToNodeIndices( 0, 2 ) = 7;
            tElementalFaceToNodeIndices( 1, 0 ) = 4;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 8;
            tElementalFaceToNodeIndices( 2, 0 ) = 4;
            tElementalFaceToNodeIndices( 2, 1 ) = 7;
            tElementalFaceToNodeIndices( 2, 2 ) = 8;
            tElementalFaceToNodeIndices( 3, 0 ) = 6;
            tElementalFaceToNodeIndices( 3, 1 ) = 7;
            tElementalFaceToNodeIndices( 3, 2 ) = 8;

            // Append with row 0
            // Append with row 1
            // Append with row 2
            // Face in row 3 exists already
            tFaceRegistry.append_face( 0, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 1, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 2, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 8;
            tExpectedElementalFaceIndices( 0, 1 ) = 9;
            tExpectedElementalFaceIndices( 0, 2 ) = 10;
            tExpectedElementalFaceIndices( 0, 3 ) = 7;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 2;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 3
            tElementalFaceToNodeIndices( 0, 0 ) = 2;
            tElementalFaceToNodeIndices( 0, 1 ) = 3;
            tElementalFaceToNodeIndices( 0, 2 ) = 6;
            tElementalFaceToNodeIndices( 1, 0 ) = 2;
            tElementalFaceToNodeIndices( 1, 1 ) = 3;
            tElementalFaceToNodeIndices( 1, 2 ) = 4;
            tElementalFaceToNodeIndices( 2, 0 ) = 3;
            tElementalFaceToNodeIndices( 2, 1 ) = 4;
            tElementalFaceToNodeIndices( 2, 2 ) = 6;
            tElementalFaceToNodeIndices( 3, 0 ) = 2;
            tElementalFaceToNodeIndices( 3, 1 ) = 4;
            tElementalFaceToNodeIndices( 3, 2 ) = 6;

            // Append faces 0,2,3 (face 1 already exists but needs to be reset)
            tFaceRegistry.append_face( 0, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 2, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 3, tElementalFaceToNodeIndices );

            CHECK( !tFaceRegistry.has_face_been_replaced( 1 ) );
            tFaceRegistry.reset_face_to_element( 1, 0 );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 11;
            tExpectedElementalFaceIndices( 0, 1 ) = 1;
            tExpectedElementalFaceIndices( 0, 2 ) = 12;
            tExpectedElementalFaceIndices( 0, 3 ) = 13;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 3;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 4
            tElementalFaceToNodeIndices( 0, 0 ) = 3;
            tElementalFaceToNodeIndices( 0, 1 ) = 4;
            tElementalFaceToNodeIndices( 0, 2 ) = 7;
            tElementalFaceToNodeIndices( 1, 0 ) = 3;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 7;
            tElementalFaceToNodeIndices( 2, 0 ) = 4;
            tElementalFaceToNodeIndices( 2, 1 ) = 6;
            tElementalFaceToNodeIndices( 2, 2 ) = 7;
            tElementalFaceToNodeIndices( 3, 0 ) = 3;
            tElementalFaceToNodeIndices( 3, 1 ) = 4;
            tElementalFaceToNodeIndices( 3, 2 ) = 6;

            // Append faces 0,1 (face 2,3 already exists)
            tFaceRegistry.append_face( 0, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 1, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 14;
            tExpectedElementalFaceIndices( 0, 1 ) = 15;
            tExpectedElementalFaceIndices( 0, 2 ) = 8;
            tExpectedElementalFaceIndices( 0, 3 ) = 12;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 4;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 1 (replaces element 1)
            tElementalFaceToNodeIndices( 0, 0 ) = 1;
            tElementalFaceToNodeIndices( 0, 1 ) = 7;
            tElementalFaceToNodeIndices( 0, 2 ) = 9;
            tElementalFaceToNodeIndices( 1, 0 ) = 1;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 7;
            tElementalFaceToNodeIndices( 2, 0 ) = 1;
            tElementalFaceToNodeIndices( 2, 1 ) = 6;
            tElementalFaceToNodeIndices( 2, 2 ) = 9;
            tElementalFaceToNodeIndices( 3, 0 ) = 6;
            tElementalFaceToNodeIndices( 3, 1 ) = 7;
            tElementalFaceToNodeIndices( 3, 2 ) = 9;

            // Replace face 6 with face 0 of element 1
            // Element Face 1 exists as global face 3
            // Replace face 4 with face 2 of element 1
            // Append Element Face 3
            tFaceRegistry.replace_face( 6, 0, tElementalFaceToNodeIndices );
            tFaceRegistry.replace_face( 4, 2, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 3, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 6;
            tExpectedElementalFaceIndices( 0, 1 ) = 3;
            tExpectedElementalFaceIndices( 0, 2 ) = 4;
            tExpectedElementalFaceIndices( 0, 3 ) = 16;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 1;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 5
            tElementalFaceToNodeIndices( 0, 0 ) = 5;
            tElementalFaceToNodeIndices( 0, 1 ) = 6;
            tElementalFaceToNodeIndices( 0, 2 ) = 9;
            tElementalFaceToNodeIndices( 1, 0 ) = 5;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 7;
            tElementalFaceToNodeIndices( 2, 0 ) = 6;
            tElementalFaceToNodeIndices( 2, 1 ) = 7;
            tElementalFaceToNodeIndices( 2, 2 ) = 9;
            tElementalFaceToNodeIndices( 3, 0 ) = 5;
            tElementalFaceToNodeIndices( 3, 1 ) = 7;
            tElementalFaceToNodeIndices( 3, 2 ) = 9;

            // Append element face 0,1,3
            tFaceRegistry.append_face( 0, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 1, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 3, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 17;
            tExpectedElementalFaceIndices( 0, 1 ) = 18;
            tExpectedElementalFaceIndices( 0, 2 ) = 16;
            tExpectedElementalFaceIndices( 0, 3 ) = 19;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 5;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 6
            tElementalFaceToNodeIndices( 0, 0 ) = 2;
            tElementalFaceToNodeIndices( 0, 1 ) = 3;
            tElementalFaceToNodeIndices( 0, 2 ) = 6;
            tElementalFaceToNodeIndices( 1, 0 ) = 2;
            tElementalFaceToNodeIndices( 1, 1 ) = 3;
            tElementalFaceToNodeIndices( 1, 2 ) = 5;
            tElementalFaceToNodeIndices( 2, 0 ) = 2;
            tElementalFaceToNodeIndices( 2, 1 ) = 5;
            tElementalFaceToNodeIndices( 2, 2 ) = 6;
            tElementalFaceToNodeIndices( 3, 0 ) = 3;
            tElementalFaceToNodeIndices( 3, 1 ) = 5;
            tElementalFaceToNodeIndices( 3, 2 ) = 6;

            // Append element face 2,3
            tFaceRegistry.append_face( 2, tElementalFaceToNodeIndices );
            tFaceRegistry.append_face( 3, tElementalFaceToNodeIndices );
            tFaceRegistry.reset_face_to_element( 5, 1 );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 11;
            tExpectedElementalFaceIndices( 0, 1 ) = 5;
            tExpectedElementalFaceIndices( 0, 2 ) = 20;
            tExpectedElementalFaceIndices( 0, 3 ) = 21;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 6;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            // Add Child Element 7
            tElementalFaceToNodeIndices( 0, 0 ) = 3;
            tElementalFaceToNodeIndices( 0, 1 ) = 5;
            tElementalFaceToNodeIndices( 0, 2 ) = 7;
            tElementalFaceToNodeIndices( 1, 0 ) = 3;
            tElementalFaceToNodeIndices( 1, 1 ) = 6;
            tElementalFaceToNodeIndices( 1, 2 ) = 7;
            tElementalFaceToNodeIndices( 2, 0 ) = 3;
            tElementalFaceToNodeIndices( 2, 1 ) = 5;
            tElementalFaceToNodeIndices( 2, 2 ) = 6;
            tElementalFaceToNodeIndices( 3, 0 ) = 5;
            tElementalFaceToNodeIndices( 3, 1 ) = 6;
            tElementalFaceToNodeIndices( 3, 2 ) = 7;

            // Append element face 0
            tFaceRegistry.append_face( 0, tElementalFaceToNodeIndices );

            // Expected face indices
            tExpectedElementalFaceIndices( 0, 0 ) = 22;
            tExpectedElementalFaceIndices( 0, 1 ) = 15;
            tExpectedElementalFaceIndices( 0, 2 ) = 21;
            tExpectedElementalFaceIndices( 0, 3 ) = 18;

            tFaceIndices = tFaceRegistry.get_face_indices( tElementalFaceToNodeIndices, true );
            CHECK( equal_to( tFaceIndices, tExpectedElementalFaceIndices ) );

            // Set face to elements
            tElementIndex( 0, 0 ) = 7;
            tFaceRegistry.set_face_to_element( tElementIndex, tFaceIndices );

            tFaceRegistry.end_modification();

            // Check the face to element
            tExpectedFaceToElement.fill( tMax );
            tExpectedFaceToElement( 0, 0 )  = 0;
            tExpectedFaceToElement( 1, 0 )  = 3;
            tExpectedFaceToElement( 2, 0 )  = 0;
            tExpectedFaceToElement( 3, 0 )  = 0;
            tExpectedFaceToElement( 3, 1 )  = 1;
            tExpectedFaceToElement( 4, 0 )  = 1;
            tExpectedFaceToElement( 5, 0 )  = 6;
            tExpectedFaceToElement( 6, 0 )  = 1;
            tExpectedFaceToElement( 7, 0 )  = 0;
            tExpectedFaceToElement( 7, 1 )  = 2;
            tExpectedFaceToElement( 8, 0 )  = 2;
            tExpectedFaceToElement( 8, 1 )  = 4;
            tExpectedFaceToElement( 9, 0 )  = 2;
            tExpectedFaceToElement( 10, 0 ) = 2;
            tExpectedFaceToElement( 11, 0 ) = 3;
            tExpectedFaceToElement( 11, 1 ) = 6;
            tExpectedFaceToElement( 12, 0 ) = 3;
            tExpectedFaceToElement( 12, 1 ) = 4;
            tExpectedFaceToElement( 13, 0 ) = 3;
            tExpectedFaceToElement( 14, 0 ) = 4;
            tExpectedFaceToElement( 15, 0 ) = 4;
            tExpectedFaceToElement( 15, 1 ) = 7;
            tExpectedFaceToElement( 16, 0 ) = 1;
            tExpectedFaceToElement( 16, 1 ) = 5;
            tExpectedFaceToElement( 17, 0 ) = 5;
            tExpectedFaceToElement( 18, 0 ) = 5;
            tExpectedFaceToElement( 18, 1 ) = 7;
            tExpectedFaceToElement( 19, 0 ) = 5;
            tExpectedFaceToElement( 20, 0 ) = 6;
            tExpectedFaceToElement( 21, 0 ) = 6;
            tExpectedFaceToElement( 21, 1 ) = 7;
            tExpectedFaceToElement( 22, 0 ) = 7;

            CHECK( equal_to( tExpectedFaceToElement, tFaceRegistry.get_face_to_element() ) );

            // Set the face ancestry (this can be done with each element)

            tFaceRegistry.set_face_ancestry( 0, 0, 2 );
            tFaceRegistry.set_face_ancestry( 1, 1, 2 );
            tFaceRegistry.set_face_ancestry( 2, 2, 2 );
            tFaceRegistry.set_face_ancestry( 3, 3, 2 );
            tFaceRegistry.set_face_ancestry( 4, 4, 2 );
            tFaceRegistry.set_face_ancestry( 5, 5, 2 );
            tFaceRegistry.set_face_ancestry( 6, 6, 2 );
            tFaceRegistry.set_face_ancestry( 7, 0, 3 );
            tFaceRegistry.set_face_ancestry( 8, 0, 3 );
            tFaceRegistry.set_face_ancestry( 9, 0, 2 );
            tFaceRegistry.set_face_ancestry( 10, 2, 2 );
            tFaceRegistry.set_face_ancestry( 11, 3, 2 );
            tFaceRegistry.set_face_ancestry( 12, 0, 3 );
            tFaceRegistry.set_face_ancestry( 13, 0, 2 );
            tFaceRegistry.set_face_ancestry( 14, 2, 2 );
            tFaceRegistry.set_face_ancestry( 15, 3, 2 );
            tFaceRegistry.set_face_ancestry( 16, 1, 3 );
            tFaceRegistry.set_face_ancestry( 17, 4, 2 );
            tFaceRegistry.set_face_ancestry( 18, 1, 3 );
            tFaceRegistry.set_face_ancestry( 19, 6, 2 );
            tFaceRegistry.set_face_ancestry( 20, 4, 2 );
            tFaceRegistry.set_face_ancestry( 21, 1, 3 );
            tFaceRegistry.set_face_ancestry( 22, 6, 2 );

            tExpectedParentRanks( 0, 0 )  = 2;
            tExpectedParentRanks( 0, 1 )  = 2;
            tExpectedParentRanks( 0, 2 )  = 2;
            tExpectedParentRanks( 0, 3 )  = 2;
            tExpectedParentRanks( 0, 4 )  = 2;
            tExpectedParentRanks( 0, 5 )  = 2;
            tExpectedParentRanks( 0, 6 )  = 2;
            tExpectedParentRanks( 0, 7 )  = 3;
            tExpectedParentRanks( 0, 8 )  = 3;
            tExpectedParentRanks( 0, 9 )  = 2;
            tExpectedParentRanks( 0, 10 ) = 2;
            tExpectedParentRanks( 0, 11 ) = 2;
            tExpectedParentRanks( 0, 12 ) = 3;
            tExpectedParentRanks( 0, 13 ) = 2;
            tExpectedParentRanks( 0, 14 ) = 2;
            tExpectedParentRanks( 0, 15 ) = 2;
            tExpectedParentRanks( 0, 16 ) = 3;
            tExpectedParentRanks( 0, 17 ) = 2;
            tExpectedParentRanks( 0, 18 ) = 3;
            tExpectedParentRanks( 0, 19 ) = 2;
            tExpectedParentRanks( 0, 20 ) = 2;
            tExpectedParentRanks( 0, 21 ) = 3;
            tExpectedParentRanks( 0, 22 ) = 2;

            CHECK( equal_to( tExpectedParentRanks, tFaceRegistry.get_face_inheritance_ranks() ) );

            tExpectedParentIndices( 0, 0 )  = 0;
            tExpectedParentIndices( 0, 1 )  = 1;
            tExpectedParentIndices( 0, 2 )  = 2;
            tExpectedParentIndices( 0, 3 )  = 3;
            tExpectedParentIndices( 0, 4 )  = 4;
            tExpectedParentIndices( 0, 5 )  = 5;
            tExpectedParentIndices( 0, 6 )  = 6;
            tExpectedParentIndices( 0, 7 )  = 0;
            tExpectedParentIndices( 0, 8 )  = 0;
            tExpectedParentIndices( 0, 9 )  = 0;
            tExpectedParentIndices( 0, 10 ) = 2;
            tExpectedParentIndices( 0, 11 ) = 3;
            tExpectedParentIndices( 0, 12 ) = 0;
            tExpectedParentIndices( 0, 13 ) = 0;
            tExpectedParentIndices( 0, 14 ) = 2;
            tExpectedParentIndices( 0, 15 ) = 3;
            tExpectedParentIndices( 0, 16 ) = 1;
            tExpectedParentIndices( 0, 17 ) = 4;
            tExpectedParentIndices( 0, 18 ) = 1;
            tExpectedParentIndices( 0, 19 ) = 6;
            tExpectedParentIndices( 0, 20 ) = 4;
            tExpectedParentIndices( 0, 21 ) = 1;
            tExpectedParentIndices( 0, 22 ) = 6;
            CHECK( equal_to( tExpectedParentIndices, tFaceRegistry.get_face_inheritance_indices() ) );
        }
    }
}    // namespace xtk
