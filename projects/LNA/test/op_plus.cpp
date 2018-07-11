// Third-party header files.
#include <catch.hpp>
#include<iostream>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
// #include <complex>      // std::complex

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_plus",
        "[linalgebra],[op_plus]" )
{

#include "linalg/op_plus.inc"

    SECTION( "moris::Mat + moris::Mat + moris::Mat" )
    {
       REQUIRE( Dm( 0, 0 ) == 6.0 );
       REQUIRE( Dm( 0, 1 ) == 6.0 );
       REQUIRE( Dm( 0, 2 ) == 6.0 );

       REQUIRE( Dm( 1, 0 ) == 6.0 );
       REQUIRE( Dm( 1, 1 ) == 6.0 );
       REQUIRE( Dm( 1, 2 ) == 6.0 );

       REQUIRE( Dm( 2, 0 ) == 6.0 );
       REQUIRE( Dm( 2, 1 ) == 6.0 );
       REQUIRE( Dm( 2, 2 ) == 6.0 );
    }

    SECTION( "complex moris::Mat + complex moris::Mat + complex moris::Mat" )
    {
       REQUIRE( Di( 0,0 ).real() == 6.0 );
       REQUIRE( Di( 0,0 ).imag() == 6.0 );

       REQUIRE( Di( 0,1 ).real() == 6.0 );
       REQUIRE( Di( 0,1 ).imag() == 6.0 );

       REQUIRE( Di( 0,2 ).real() == 6.0 );
       REQUIRE( Di( 0,2 ).imag() == 6.0 );

       REQUIRE( Di( 1,0 ).real() == 6.0 );
       REQUIRE( Di( 1,0 ).imag() == 6.0 );

       REQUIRE( Di( 1,1 ).real() == 6.0 );
       REQUIRE( Di( 1,1 ).imag() == 6.0 );

       REQUIRE( Di( 1,2 ).real() == 6.0 );
       REQUIRE( Di( 1,2 ).imag() == 6.0 );

       REQUIRE( Di( 2,0 ).real() == 6.0 );
       REQUIRE( Di( 2,0 ).imag() == 6.0 );

       REQUIRE( Di( 2,1 ).real() == 6.0 );
       REQUIRE( Di( 2,1 ).imag() == 6.0 );

       REQUIRE( Di( 2,2 ).real() == 6.0 );
       REQUIRE( Di( 2,2 ).imag() == 6.0 );
    }
/*
    SECTION( "complex moris::Mat + complex moris::Mat + complex moris::Mat" )
    {
       REQUIRE( moris::equal_to( Di( 0,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 0,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 0,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 0,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 0,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 0,2 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 1,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 1,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 1,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 1,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 1,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 1,2 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 2,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 2,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 2,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 2,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Di( 2,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Di( 2,2 ).imag(),  6.0 ) );
    }
*/
    SECTION( "moris::Sp_Mat + moris::Sp_Mat + moris::Sp_Mat" )
    {
       REQUIRE( Ds( 0, 0 ) == 6.0 );
       REQUIRE( Ds( 0, 1 ) == 6.0 );
       REQUIRE( Ds( 0, 2 ) == 6.0 );

       REQUIRE( Ds( 1, 0 ) == 6.0 );
       REQUIRE( Ds( 1, 1 ) == 6.0 );
       REQUIRE( Ds( 1, 2 ) == 6.0 );

       REQUIRE( Ds( 2, 0 ) == 6.0 );
       REQUIRE( Ds( 2, 1 ) == 6.0 );
       REQUIRE( Ds( 2, 2 ) == 6.0 );
    }
/*
    SECTION( "complex moris::Sp_Mat + complex moris::Sp_Mat + complex moris::Sp_Mat" )
    {
       REQUIRE( Dsi( 0,0 ).real() == 6.0 );
       REQUIRE( Dsi( 0,0 ).imag() == 6.0 );

       REQUIRE( Dsi( 0,1 ).real() == 6.0 );
       REQUIRE( Dsi( 0,1 ).imag() == 6.0 );

       REQUIRE( Dsi( 0,2 ).real() == 6.0 );
       REQUIRE( Dsi( 0,2 ).imag() == 6.0 );

       REQUIRE( Dsi( 1,0 ).real() == 6.0 );
       REQUIRE( Dsi( 1,0 ).imag() == 6.0 );

       REQUIRE( Dsi( 1,1 ).real() == 6.0 );
       REQUIRE( Dsi( 1,1 ).imag() == 6.0 );

       REQUIRE( Dsi( 1,2 ).real() == 6.0 );
       REQUIRE( Dsi( 1,2 ).imag() == 6.0 );

       REQUIRE( Dsi( 2,0 ).real() == 6.0 );
       REQUIRE( Dsi( 2,0 ).imag() == 6.0 );

       REQUIRE( Dsi( 2,1 ).real() == 6.0 );
       REQUIRE( Dsi( 2,1 ).imag() == 6.0 );

       REQUIRE( Dsi( 2,2 ).real() == 6.0 );
       REQUIRE( Dsi( 2,2 ).imag() == 6.0 );
    }
    SECTION( "complex moris::Sp_Mat + complex moris::Sp_Mat + complex moris::Sp_Mat" )
    {
       REQUIRE( moris::equal_to( Dsi( 0,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 0,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 0,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 0,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 0,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 0,2 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 1,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 1,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 1,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 1,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 1,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 1,2 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 2,0 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 2,0 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 2,1 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 2,1 ).imag(),  6.0 ) );

       REQUIRE( moris::equal_to( Dsi( 2,2 ).real(),  6.0 ) );
       REQUIRE( moris::equal_to( Dsi( 2,2 ).imag(),  6.0 ) );
    }
*/
    SECTION( "moris::Col + moris::Col + moris::Col" )
    {
        REQUIRE( Dc( 0, 0 ) == 6.0 );
        REQUIRE( Dc( 1, 0 ) == 6.0 );
        REQUIRE( Dc( 2, 0 ) == 6.0 );

        REQUIRE( moris::equal_to( Dc( 0, 0 ), 6.0 ) );
        REQUIRE( moris::equal_to( Dc( 1, 0 ), 6.0 ) );
        REQUIRE( moris::equal_to( Dc( 2, 0 ), 6.0 ) );
    }

    SECTION( "complex moris::Col + complex moris::Col + complex moris::Col" )
    {
        REQUIRE( Dci( 0,0 ).real() == 6.0 );
        REQUIRE( Dci( 0,0 ).imag() == 6.0 );

        REQUIRE( Dci( 1,0 ).real() == 6.0 );
        REQUIRE( Dci( 1,0 ).imag() == 6.0 );

        REQUIRE( Dci( 2,0 ).real() == 6.0 );
        REQUIRE( Dci( 2,0 ).imag() == 6.0 );
    }

}
