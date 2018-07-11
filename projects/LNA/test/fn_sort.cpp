/*
 * fn_sort.cpp
 *
 *  Created on: Feb 20, 2017
 *      Author: doble
 */
// Third-party header files.
#include <catch.hpp>
#include <iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
TEST_CASE(
        "moris::sort",
        "[linalgebra],[sort]" )
{
#include "linalg/fn_sort.inc"


    REQUIRE(moris::equal_to(tMatSort(0,0),3));
    REQUIRE(moris::equal_to(tMatSort(1,0),4));
}



