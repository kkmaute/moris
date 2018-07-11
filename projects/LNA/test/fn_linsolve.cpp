// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::linsolve",
        "[linalgebra],[linsolve],[Mat]" )
{
    SECTION(
            "solve( moris::Mat, moris::Mat )" )
    {
        #include "linalg/fn_linsolve/linsolve_real.inc"

        REQUIRE( moris::equal_to ( C( 0,0 ),  4.5 ) );
        REQUIRE( moris::equal_to ( C( 1,0 ),  2.0 ) );
        REQUIRE( moris::equal_to ( C( 2,0 ), -4.5 ) );
    }

    SECTION(
            "solve( complex moris::Mat, complex moris::Mat )" )
    {
        #include "linalg/fn_linsolve/linsolve_complex.inc"

        REQUIRE( moris::equal_to( C( 0,0 ).real(),  4.529729729729730 ) );
        REQUIRE( moris::equal_to( C( 0,0 ).imag(), -3.778378378378378 ) );
        REQUIRE( moris::equal_to( C( 1,0 ).real(), -3.772972972972973 ) );
        REQUIRE( moris::equal_to( C( 1,0 ).imag(),  1.237837837837838 ) );
    }

    SECTION(
            "solve( moris::Sp_Mat, moris::Sp_Mat )" )
    {
//        moris::Sp_Mat< moris::real > A( 3, 3 );
//        moris::Mat< moris::real > B( 3, 1 );
//
//        A( 0, 0 ) = 3.0; A( 0, 1 ) = 0.0; A( 0, 2 ) = 3.0;
//        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
//        A( 2, 0 ) = 9.0; A( 2, 1 ) = 1.0; A( 2, 2 ) = 9.0;
//
//        B( 0, 0 ) = 0.0;
//        B( 1, 0 ) = 1.0;
//        B( 2, 0 ) = 2.0;
//
//        moris::Mat< moris::real > C = moris::solve( A, B, "superlu" );
//
//        moris::Mat< moris::real > D = moris::solve( A, B, "sparselu" );
//
//        moris::Mat< moris::real > E = moris::solve( A, B, "umfpack" );
//
//        REQUIRE( moris::equal_to ( C( 0,0 ),  4.5 ) );
//        REQUIRE( moris::equal_to ( C( 1,0 ),  2.0 ) );
//        REQUIRE( moris::equal_to ( C( 2,0 ), -4.5 ) );
//
//        REQUIRE( moris::equal_to ( D( 0,0 ),  4.5 ) );
//        REQUIRE( moris::equal_to ( D( 1,0 ),  2.0 ) );
//        REQUIRE( moris::equal_to ( D( 2,0 ), -4.5 ) );
//
//        REQUIRE( moris::equal_to ( E( 0,0 ),  4.5 ) );
//        REQUIRE( moris::equal_to ( E( 1,0 ),  2.0 ) );
//        REQUIRE( moris::equal_to ( E( 2,0 ), -4.5 ) );
    }

    SECTION(
            "solve( complex moris::Sp_Mat, complex moris::Sp_Mat )" )
    {
//        moris::Sp_Mat< moris::real > A( 3, 3 );
//        moris::Mat< moris::real > B( 3, 1 );
//
//        A( 0, 0 ) = 3.0; A( 0, 1 ) = 0.0; A( 0, 2 ) = 3.0;
//        A( 1, 0 ) = 4.0; A( 1, 1 ) = 5.0; A( 1, 2 ) = 6.0;
//        A( 2, 0 ) = 9.0; A( 2, 1 ) = 1.0; A( 2, 2 ) = 9.0;
//
//        B( 0, 0 ) = 0.0;
//        B( 1, 0 ) = 1.0;
//        B( 2, 0 ) = 2.0;
//
//        moris::Mat< moris::real > C = moris::solve( A, B, "superlu" );
//
//        std::complex<double> m00 = C( 0, 0 );
//        std::complex<double> m10 = C( 1, 0 );
//
//        REQUIRE( moris::equal_to ( m00.real(),  4.5 ) );
//        REQUIRE( moris::equal_to ( m00.imag(),  0.0 ) );
//        REQUIRE( moris::equal_to ( m10.real(),  2.0 ) );
//        REQUIRE( moris::equal_to ( m10.imag(),  0.0 ) );
//
//        moris::Mat< moris::real > D = moris::solve( A, B, "sparselu" );
//
//        std::complex<double> n00 = D( 0, 0 );
//        std::complex<double> n10 = D( 1, 0 );
//
//        REQUIRE( moris::equal_to ( n00.real(),  4.5 ) );
//        REQUIRE( moris::equal_to ( n00.imag(),  0.0 ) );
//        REQUIRE( moris::equal_to ( n10.real(),  2.0 ) );
//        REQUIRE( moris::equal_to ( n10.imag(),  0.0 ) );
    }
}
