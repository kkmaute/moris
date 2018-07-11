// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::op_elemwise_mult",
        "[linalgebra],[elemwise_mult],[Mat]" )
{
    SECTION(
            "moris::Mat element-wise multiplication with moris::Mat" )
    {
        moris::Mat< moris::real > A( 2, 2 );
        moris::Mat< moris::real > B( 2, 2 );

        A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0;
        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0;

        B( 0, 0 ) = 3.0; B( 0, 1 ) = 6.0;
        B( 1, 0 ) = 1.0; B( 1, 1 ) = 7.0;

        moris::Mat< moris::real > C( 2, 2 );
        C = A % B;

        REQUIRE( C( 0,0 ) == 3.0 ); REQUIRE( C( 0,1 ) == 12.0 );
        REQUIRE( C( 1,0 ) == 4.0 ); REQUIRE( C( 1,1 ) == 35.0 );
    }

    SECTION(
            "complex moris::Mat elem-wise multiplic with complex moris::Mat" )
    {
        #include "linalg/op_elemwise_mult.inc"

        REQUIRE( moris::equal_to( C( 0,0 ).real(), -11.0 ) );
        REQUIRE( moris::equal_to( C( 0,0 ).imag(),  23.0 ) );
        REQUIRE( moris::equal_to( C( 0,1 ).real(),  0.0 ) );
        REQUIRE( moris::equal_to( C( 0,1 ).imag(),  2.0 ) );
    }
}
