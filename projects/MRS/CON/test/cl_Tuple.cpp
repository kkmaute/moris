// C++ header files.
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

    REQUIRE(b.get< 0 >() == 2 );
    REQUIRE(b.get< 1 >() == 2.0 );
    REQUIRE(b.get< 2 >() == "hello" );

    moris::Tuple< moris::sint, moris::real, std::string > c( b );

    REQUIRE(c.get< 0 >() == 2 );
    REQUIRE(c.get< 1 >() == 2.0 );
    REQUIRE(c.get< 2 >() == "hello" );

    moris::Tuple< moris::sint, moris::real, std::string > d = b;

    REQUIRE(d.get< 0 >() == 2 );
    REQUIRE(d.get< 0 >() >  1 );
    REQUIRE(d.get< 1 >() == 2.0 );
    REQUIRE(d.get< 2 >() == "hello" );

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


    int numOpts = 1000;


    moris::tic timer2;
    for(int i2 = 0; i2 < numOpts; ++i2)
    {
        moris::uint aInt21 = 4;
        moris::uint aInt22 = 5;
        moris::uint aInt23 = 6;

        moris::uint aInt24, aInt25, aInt26;

        std::tie(aInt24, aInt25, aInt26) = std::make_tuple(aInt21, aInt22, aInt23);
    }
    moris::real time2 = timer2.toc<moris::chronos::nanoseconds>().wall;


    moris::tic timer0;
    for(int i0 = 0; i0 < numOpts; ++i0)
    {
        moris::uint aInt01 = 4;
        moris::uint aInt02 = 5;
        moris::uint aInt03 = 6;

        moris::uint aInt04, aInt05, aInt06;

        moris::tie(aInt04, aInt05, aInt06) = moris::make_tuple(aInt01, aInt02, aInt03);
    }
    moris::real time0 = timer0.toc<moris::chronos::nanoseconds>().wall;


    moris::tic timer1;
    for(int i1 = 0; i1 < numOpts; ++i1)
    {
        moris::uint aInt11 = 4;
        moris::uint aInt12 = 5;
        moris::uint aInt13 = 6;

        moris::uint aInt14, aInt15, aInt16;

        std::tie(aInt14, aInt15, aInt16) = std::make_tuple(aInt11, aInt12, aInt13);
    }
    moris::real time1 = timer1.toc<moris::chronos::nanoseconds>().wall;

    //REQUIRE( ( time0 - time1) < 0.15*time1 );

}
