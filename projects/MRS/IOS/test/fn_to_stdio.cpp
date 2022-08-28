/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_to_stdio.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "ios.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::ios::to_stdio",
        "[moris],[ios],[fn_to_stdio],[to_stdio]")
{
    #include "fn_to_stdio.inc" // snippets IOS

    REQUIRE(read_mode == "r");
    REQUIRE(write_mode == "w");
    REQUIRE(append_mode == "a");
    REQUIRE(readwrite_mode == "r+");
    REQUIRE(trunc_mode == "w+");
    REQUIRE(update_mode == "a+");
}

