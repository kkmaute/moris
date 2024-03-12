/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Downward_Inheritance.cpp
 *
 */

#include <utility>

#include "cl_Vector.hpp"
#include "cl_XTK_Downward_Inheritance.hpp"
#include "catch.hpp"

TEST_CASE( "Downward Inheritance",
        "[XTK][INHERITANCE]" )
{
    /*
     * Tests the Downward inheritance structure to see if pairs are registered and stored correctly
     * and that existing pairs are not overwritten
     */
    moris::xtk::Downward_Inheritance< moris::size_t, moris::size_t > tInheritance( 5 );
    tInheritance.register_new_inheritance_pair( 1, 6 );
    tInheritance.register_new_inheritance_pair( 1, 7 );
    CHECK( tInheritance.has_inheritance( 1 ) );
    CHECK( !tInheritance.has_inheritance( 2 ) );
    CHECK( tInheritance.get_inheritance( 1 ) == 6 );    // Make sure 7 did not overwrite 6
}
