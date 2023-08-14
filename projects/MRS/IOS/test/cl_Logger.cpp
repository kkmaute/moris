/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Logger.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

#include "cl_Logger.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "MORIS_LOG_TEST",
        "[moris],[ios],[cl_Logger],[MORIS_LOG],[log_test]")
{ 
        MORIS_SECTION("                           Test Section Number  %-5i  \n", 1 );

        MORIS_SECTION( " Test check severity level warning " );

        MORIS_LOG_WARNING("Test check severity level warning  || %-5i   ", 44 );
}

