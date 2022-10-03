/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * assert.cpp
 *
 */

#include <catch.hpp>
#include "assert.hpp"
#include "cl_Communication_Tools.hpp" // COM/src

extern moris::Logger gLogger;

namespace moris
{

void test_moris_error( bool aCheck )
{
    MORIS_ERROR( aCheck , "Test output message for MORIS ERROR throw" );
}

void test_moris_assert( bool aCheck)
{
    MORIS_ASSERT( aCheck, "Test output message for MORIS ASSERT throw" );
}

TEST_CASE(
        "MORIS_ASSERT and MORIS_ERROR Tests",
        "[MORIS],[assert],[error]"      )
{

    SECTION("CHECK FOR THROWS")
    {
	// Verify that MORIS_ERROR throws. Otherwise all CHECK_THROW tests will fail.
        CHECK_THROWS( test_moris_error( false ) );

        // Verify that MORIS_ERROR throws only in debug mode. Otherwise all CHECK_THROW tests will fail.
#ifdef MORIS_HAVE_DEBUG
        CHECK_THROWS(test_moris_assert( false ));
#endif
    }

    SECTION("CHECK FOR NO THROWS")
    {
        // Verify that MORIS_ERROR does not throw.
        CHECK_NOTHROW( test_moris_error( true ) );

        // Verify that MORIS_ASSERT does not throw.
        CHECK_NOTHROW( test_moris_assert( true  ) );
    }

  }

 }
