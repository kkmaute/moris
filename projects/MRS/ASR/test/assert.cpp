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
TEST_CASE(
        "MORIS_ASSERT and MORIS_ERROR Tests",
        "[MORIS],[assert],[error]"      )
                                {

    SECTION("CHECK FOR THROWS")
    {
        // Verify that MORIS_ERROR throws. Otherwise all CHECK_THROW tests will fail.
        CHECK_THROWS(MORIS_ERROR( 1 == 0, "Test output message for MORIS ERROR throw" ));

        // Verify that MORIS_ERROR throws only in debug mode. Otherwise all CHECK_THROW tests will fail.
#ifdef MORIS_HAVE_DEBUG
        CHECK_THROWS(MORIS_ASSERT( 1 == 0, "Test output message for MORIS ASSERT throw" ));
#endif
    }

    SECTION("CHECK FOR NO THROWS")
    {
        // Verify that MORIS_ERROR does not throw.
        CHECK_NOTHROW(MORIS_ERROR( 1 == 1, "Test output message for MORIS ERROR throw" ));

        // Verify that MORIS_ASSERT does not throw.
        CHECK_NOTHROW(MORIS_ASSERT( 1 == 1, "Test output message for MORIS ASSERT throw" ));
    }

    SECTION("CHECK FOR THROWS WITH MULTIPLE ARGUMENTS")
    {
        // Verify that MORIS_ERROR throws. Otherwise all CHECK_THROW tests will fail.
        CHECK_THROWS(MORIS_ERROR( false, "Test output message for MORIS ERROR throw %-5i   ||   %-5.2e", 124, 0.00456789 ));

        // Verify that MORIS_ERROR throws only in debug mode. Otherwise all CHECK_THROW tests will fail.
#ifdef MORIS_HAVE_DEBUG
        CHECK_THROWS(MORIS_ASSERT( false, "Test output message for MORIS ERROR throw %-5i   ||   %-5.2e", 124, 0.00456789 ));
#endif
    }
                                }
}

