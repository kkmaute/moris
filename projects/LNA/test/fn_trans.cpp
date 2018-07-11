// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::trans",
        "[linalgebra],[trans],[Mat]" )
{
    SECTION( "moris::Mat'" )
    {
        #include "linalg/fn_trans/trans_real.inc"

        REQUIRE( B( 0,0 ) == 1.0 ); REQUIRE( B( 0,1 ) == 4.0 ); REQUIRE( B( 0,2 ) == 7.0 );
        REQUIRE( B( 1,0 ) == 2.0 ); REQUIRE( B( 1,1 ) == 5.0 ); REQUIRE( B( 1,2 ) == 8.0 );
        REQUIRE( B( 2,0 ) == 3.0 ); REQUIRE( B( 2,1 ) == 6.0 ); REQUIRE( B( 2,2 ) == 9.0 );
    }

    SECTION( "complex moris::Mat'" )
    {
        #include "linalg/fn_trans/trans_complex.inc"

        REQUIRE( moris::equal_to( B( 0,1 ).real(),  2.0 ) );
        REQUIRE( moris::equal_to( B( 0,1 ).imag(), -1.0 ) );
        REQUIRE( moris::equal_to( B( 1,0 ).real(),  4.0 ) );
        REQUIRE( moris::equal_to( B( 1,0 ).imag(),  5.0 ) );
    }

    SECTION( "moris::Sp_Mat'" )
    {
        moris::Sp_Mat< moris::real > A( 3, 3 );

        A( 0, 0 ) = 1.0; A( 0, 1 ) = 2.0; A( 0, 2 ) = 3.0;
        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
        A( 2, 0 ) = 7.0; A( 2, 1 ) = 8.0; A( 2, 2 ) = 9.0;

        moris::Sp_Mat< moris::real > B( 3, 3 );
        B = moris::trans( A );

        REQUIRE( B( 0,0 ) == 1.0 ); REQUIRE( B( 0,1 ) == 4.0 ); REQUIRE( B( 0,2 ) == 7.0 );
        REQUIRE( B( 1,0 ) == 2.0 ); REQUIRE( B( 1,1 ) == 5.0 ); REQUIRE( B( 1,2 ) == 8.0 );
        REQUIRE( B( 2,0 ) == 3.0 ); REQUIRE( B( 2,1 ) == 6.0 ); REQUIRE( B( 2,2 ) == 9.0 );
    }

    SECTION( "complex moris::Sp_Mat'" )
    {
        moris::Sp_Mat< moris::cplx > A( 2, 2 );

        A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0};
        A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0};

        moris::Sp_Mat< moris::cplx > B( 2, 2 );
        B = moris::trans( A );

        std::complex<double> m01 = B( 0, 1 );
        std::complex<double> m10 = B( 1, 0 );

        REQUIRE( m01.real() ==  2.0 );
        REQUIRE( m01.imag() == -1.0 );
        REQUIRE( m10.real() ==  4.0 );
        REQUIRE( m10.imag() ==  5.0 );
    }
}
