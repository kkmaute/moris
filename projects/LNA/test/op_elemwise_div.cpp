// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_elemwise_div",
        "[linalgebra],[elemwise_div],[Mat]" )
{

#include "linalg/op_elemwise_div.inc"

    SECTION(
            "moris::Mat element-wise division with moris::Mat" )
    {
        REQUIRE( Cm( 0,0 ) == 3.0 ); REQUIRE( Cm( 0,1 ) == 2.0 );
        REQUIRE( Cm( 1,0 ) == 4.0 ); REQUIRE( Cm( 1,1 ) == 5.0 );
    }

    SECTION(
            "moris::Mat element-wise division with constant" )
    {
        REQUIRE( Csc( 0,0 ) == 3.0 );   REQUIRE( Csc( 0,1 ) == 4.0 );
        REQUIRE( Csc( 1,0 ) == 100.0 ); REQUIRE( Csc( 1,1 ) == 10.0 );
    }

    SECTION(
            "complex moris::Mat element-wise division with complex moris::Mat" )
    {
        moris::Mat< moris::cplx > A( 2, 2 );
        moris::Mat< moris::cplx > B( 2, 2 );

        A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0};
        A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0};

        B( 0, 0 ) = { 4.0, -3.0}; B( 0, 1 ) = { 2.0, 6.0};
        B( 1, 0 ) = {-3.0,  6.0}; B( 1, 1 ) = { 4.0, 6.0};

        moris::Mat< moris::cplx > C( 2, 2 );
        C = A / B;

        REQUIRE( moris::equal_to( C( 0,0 ).real(),  0.24 ) );
        REQUIRE( moris::equal_to( C( 0,0 ).imag(),  0.68 ) );
        REQUIRE( moris::equal_to( C( 0,1 ).real(),  0.95 ) );
        REQUIRE( moris::equal_to( C( 0,1 ).imag(), -0.35 ) );
        REQUIRE( moris::equal_to( C( 1,0 ).imag(), -0.20 ) );
    }

    SECTION(
            "complex moris::Mat element-wise division with complex constant" )
    {
        moris::Mat< moris::cplx > A( 2, 2 );

        A( 0, 0 ) = { 3.0, 2.0}; A( 0, 1 ) = { 3.0, 2.0};
        A( 1, 0 ) = { 3.0, 2.0}; A( 1, 1 ) = { 3.0, 2.0};

        moris::cplx B = {4.0, -3.0};

        moris::Mat< moris::cplx > C( 2, 2 );
        C = A / B;

        REQUIRE( moris::equal_to( C( 0,0 ).real(), 0.24 ) );
        REQUIRE( moris::equal_to( C( 1,1 ).imag(), 0.68 ) );
    }
}
