/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Mesh_Side_Set_Input.cpp
 *
 */

#include "catch.hpp"

#include "mesh/cl_Mesh_Side_Set_Input.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

TEST_CASE( "Side Set Input", "[SIDE_SET_INPUT]" )
{

    /*
     * Functions tested
     * - add_element_id_and_side_ordinal
     * - set_side_set_name (including a throw check)
     * - get_num_of_sides
     * - get_num_side_sets (should be 1)
     * - get_side_set_name
     * - get_element_id
     * - get_side_ordinal
     * - has_sides
     */
    size_t                                         tNumSides     = 3;
    size_t                                         tStringLength = 24;
    mesh::Side_Set_Input< size_t, moris::DDSTMat > tSideSetInput( tNumSides, tStringLength );

    CHECK( !tSideSetInput.has_sides() );
    /*
     * Add element and side ordinal pairs (note adding 3 here to fill space allocated in constructor)
     */
    moris::Matrix< moris::DDSTMat > tElementIds( { { 1, 4, 6 } } );
    moris::Matrix< moris::DDSTMat > tSideOrdinals( { { 0, 3, 1 } } );

    tSideSetInput.add_element_id_and_side_ordinal( tElementIds, tSideOrdinals );

    CHECK( tSideSetInput.has_sides() );

    /*
     * Add the side set name (Note this is less than the length above)
     */

    std::string tSideSetName = "test_side_set";

    tSideSetInput.set_side_set_name( tSideSetName );

    /*
     * Do it again and expect a throw
     */
    CHECK_THROWS( tSideSetInput.set_side_set_name( tSideSetName ) );

    /*
     * Test accessing functions
     */

    CHECK( tSideSetName.compare( tSideSetInput.get_side_set_name() ) == 0 );

    CHECK( tSideSetInput.get_element_id( 0 ) == 1 );
    CHECK( tSideSetInput.get_element_id( 1 ) == 4 );
    CHECK( tSideSetInput.get_element_id( 2 ) == 6 );
    CHECK( tSideSetInput.get_side_ordinal( 0 ) == 0 );
    CHECK( tSideSetInput.get_side_ordinal( 1 ) == 3 );
    CHECK( tSideSetInput.get_side_ordinal( 2 ) == 1 );

    CHECK_THROWS( tSideSetInput.get_element_id( 4 ) );
    CHECK_THROWS( tSideSetInput.get_side_ordinal( 4 ) );

    CHECK( tSideSetInput.get_num_side_sets() == 1 );
}
