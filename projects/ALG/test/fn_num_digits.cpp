/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_num_digits.cpp
 *
 */

#include <catch.hpp>

// MORIS header files.
#include "algorithms.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::numeric::num_digits",
        "[moris],[numeric],[fn_num_digits],[num_digits]")
{
    #include "fn_num_digits.inc" // snippets ALG

    REQUIRE(num_digits_is_1 == 1);
    REQUIRE(num_digits_is_2 == 2);
    REQUIRE(num_digits_is_3 == 3);
    REQUIRE(num_digits_is_4 == 4);

    REQUIRE_FALSE(num_digits_is_324 == 1);
    REQUIRE_FALSE(num_digits_is_325 == 2);
    REQUIRE_FALSE(num_digits_is_326 == 3);
    REQUIRE_FALSE(num_digits_is_327 == 4);

    REQUIRE(neg_num_digits_is_3 == 3);
}

