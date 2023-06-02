/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Tuple.cpp
 *
 */

#include <iostream>
#include <string>

// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "core.hpp"
#include "cl_Tuple.hpp" // CON/src
#include "chronos.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::Tuple" )
{
    #include "containers/cl_tuple.inc"

    // check access to components of tuple
    REQUIRE(b.get< 0 >() == 2 );
    REQUIRE(b.get< 1 >() == 2.0 );
    REQUIRE(b.get< 2 >() == "hello" );

    // check initialize by argument
    moris::Tuple< moris::sint, moris::real, std::string > c( b );

    REQUIRE(c.get< 0 >() == 2 );
    REQUIRE(c.get< 1 >() == 2.0 );
    REQUIRE(c.get< 2 >() == "hello" );

    // check initialize by copy operator
    moris::Tuple< moris::sint, moris::real, std::string > d = b;

    REQUIRE(d.get< 0 >() == 2 );
    REQUIRE(d.get< 0 >() >  1 );
    REQUIRE(d.get< 1 >() == 2.0 );
    REQUIRE(d.get< 2 >() == "hello" );

    // check equality and inequality of tuples
    bool check1 = ( b == c);

    moris::Tuple< int, moris::real, std::string > bb( 3, 2.0, "hello" );
    moris::Tuple< int, moris::real, std::string > cc( 2, 5.0, "hello" );
    moris::Tuple< int, moris::real, std::string > dd( 2, 2.0, "goodbye" );

    bool check2 = (b == bb);
    bool check3 = (b == cc);
    bool check4 = (b == dd);

    REQUIRE( check1 == true);
    REQUIRE( check2 == false);
    REQUIRE( check3 == false);
    REQUIRE( check4 == false);

    moris::uint aInt21 = 4;
    moris::uint aInt22 = 5;
    moris::uint aInt23 = 6;

    moris::uint aInt24, aInt25, aInt26;

    // check std library functions
    std::tie(aInt24, aInt25, aInt26) = std::make_tuple(aInt21, aInt22, aInt23);

    REQUIRE( aInt24 == aInt21);
    REQUIRE( aInt25 == aInt22);
    REQUIRE( aInt26 == aInt23);
}

