// Third-party header files.
#include <catch.hpp>

#include <cmath>

// MORIS project header files.
#include "chronos.hpp"
#include "algorithms.hpp"
#include "linalg.hpp"

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::qr",
        "[linalgebra],[qr]" )
{
    SECTION( "qr(m-by-n mat)" )
    {
        #include "linalg/fn_qr.inc"


        REQUIRE( Q.numel() == 16 );
        REQUIRE( moris::equal_to( Q(0,0), -0.5) ); REQUIRE( moris::equal_to( Q(0,1),  0.5) ); REQUIRE( moris::equal_to( Q(0,2), -0.5) ); REQUIRE( moris::equal_to( Q(0,3), -0.5) );
        REQUIRE( moris::equal_to( Q(1,0), -0.5) ); REQUIRE( moris::equal_to( Q(1,1), -0.5) ); REQUIRE( moris::equal_to( Q(1,2),  0.5) ); REQUIRE( moris::equal_to( Q(1,3), -0.5) );
        REQUIRE( moris::equal_to( Q(2,0), -0.5) ); REQUIRE( moris::equal_to( Q(2,1), -0.5) ); REQUIRE( moris::equal_to( Q(2,2), -0.5) ); REQUIRE( moris::equal_to( Q(2,3),  0.5) );
        REQUIRE( moris::equal_to( Q(3,0), -0.5) ); REQUIRE( moris::equal_to( Q(3,1),  0.5) ); REQUIRE( moris::equal_to( Q(3,2),  0.5) ); REQUIRE( moris::equal_to( Q(3,3),  0.5) );

        moris::Mat< moris::real > I = Q * moris::trans(Q);

        REQUIRE( I.numel() == 16 );
        REQUIRE( std::abs( I(0,0) - 1.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(0,1) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(0,2) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(0,3) - 0.0 ) < 1.0e-12 );
        REQUIRE( std::abs( I(1,0) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(1,1) - 1.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(1,2) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(1,3) - 0.0 ) < 1.0e-12 );
        REQUIRE( std::abs( I(2,0) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(2,1) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(2,2) - 1.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(2,3) - 0.0 ) < 1.0e-12 );
        REQUIRE( std::abs( I(3,0) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(3,1) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(3,2) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(3,3) - 1.0 ) < 1.0e-12 );

        REQUIRE( R.numel() == 12 );
        REQUIRE( std::abs( R(0,0) - -2.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(0,1) - -3.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(0,2) - -2.0 ) < 1.0e-12 );
        REQUIRE( std::abs( R(1,0) -  0.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(1,1) - -5.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(1,2) -  2.0 ) < 1.0e-12 );
        REQUIRE( std::abs( R(2,0) -  0.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(2,1) -  0.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(2,2) - -4.0 ) < 1.0e-12 );
        REQUIRE( std::abs( R(3,0) -  0.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(3,1) -  0.0 ) < 1.0e-12 ); REQUIRE( std::abs( R(3,2) -  0.0 ) < 1.0e-12 );
    }

    SECTION( "qr(n-by-n mat)" )
    {
        moris::Mat< moris::real > a( 3, 3 );
        moris::Mat< moris::real > I;

        a( 0, 0 ) = 12.0; a( 0, 1 ) = -51.0; a( 0, 2 ) =   4.0;
        a( 1, 0 ) =  6.0; a( 1, 1 ) = 167.0; a( 1, 2 ) = -68.0;
        a( 2, 0 ) = -4.0; a( 2, 1 ) =  24.0; a( 2, 2 ) = -41.0;

        // initializes Q and R, then calls tie(Q,R)=qr(a)
        #include "linalg/fn_qr_tuple.inc"

        REQUIRE( Q.numel() == 9 );
        REQUIRE( std::abs( Q(0,0) - (-0.857142857142857) ) < 1.0e-12 ); REQUIRE( std::abs( Q(0,1) - ( 0.394285714285714) ) < 1.0e-12 ); REQUIRE( std::abs( Q(0,2) - ( 0.331428571428571) ) < 1.0e-12 );
        REQUIRE( std::abs( Q(1,0) - (-0.428571428571429) ) < 1.0e-12 ); REQUIRE( std::abs( Q(1,1) - (-0.902857142857143) ) < 1.0e-12 ); REQUIRE( std::abs( Q(1,2) - (-0.034285714285714) ) < 1.0e-12 );
        REQUIRE( std::abs( Q(2,0) - ( 0.285714285714286) ) < 1.0e-12 ); REQUIRE( std::abs( Q(2,1) - (-0.171428571428571) ) < 1.0e-12 ); REQUIRE( std::abs( Q(2,2) - ( 0.942857142857143) ) < 1.0e-12 );

        REQUIRE( R.numel() == 9 );
        REQUIRE( std::abs( R(0,0) - (-14.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(0,1) - ( -21.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(0,2) - (  14.0) ) < 1.0e-12 );
        REQUIRE( std::abs( R(1,0) - (  0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(1,1) - (-175.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(1,2) - (  70.0) ) < 1.0e-12 );
        REQUIRE( std::abs( R(2,0) - (  0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(2,1) - (   0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(2,2) - ( -35.0) ) < 1.0e-12 );

        I = moris::trans(Q) * Q;
        REQUIRE( I.numel() == 9 );
        REQUIRE( std::abs( I(0,0) - 1.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(0,1) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(0,2) - 0.0 ) < 1.0e-12 );
        REQUIRE( std::abs( I(1,0) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(1,1) - 1.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(1,2) - 0.0 ) < 1.0e-12 );
        REQUIRE( std::abs( I(2,0) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(2,1) - 0.0 ) < 1.0e-12 ); REQUIRE( std::abs( I(2,2) - 1.0 ) < 1.0e-12 );
    }

    SECTION( "qr(n-by-n mat)" )
    {
        moris::Mat< moris::real > a( 3, 3 );
        moris::Mat< moris::real > Q;

        a( 0, 0 ) = 12.0; a( 0, 1 ) = -51.0; a( 0, 2 ) =   4.0;
        a( 1, 0 ) =  6.0; a( 1, 1 ) = 167.0; a( 1, 2 ) = -68.0;
        a( 2, 0 ) = -4.0; a( 2, 1 ) =  24.0; a( 2, 2 ) = -41.0;

        moris::tie(Q,std::ignore) = moris::qr(a);

        REQUIRE( Q.numel() == 9 );
        REQUIRE( std::abs( Q(0,0) - (-0.857142857142857) ) < 1.0e-12 ); REQUIRE( std::abs( Q(0,1) - ( 0.394285714285714) ) < 1.0e-12 ); REQUIRE( std::abs( Q(0,2) - ( 0.331428571428571) ) < 1.0e-12 );
        REQUIRE( std::abs( Q(1,0) - (-0.428571428571429) ) < 1.0e-12 ); REQUIRE( std::abs( Q(1,1) - (-0.902857142857143) ) < 1.0e-12 ); REQUIRE( std::abs( Q(1,2) - (-0.034285714285714) ) < 1.0e-12 );
        REQUIRE( std::abs( Q(2,0) - ( 0.285714285714286) ) < 1.0e-12 ); REQUIRE( std::abs( Q(2,1) - (-0.171428571428571) ) < 1.0e-12 ); REQUIRE( std::abs( Q(2,2) - ( 0.942857142857143) ) < 1.0e-12 );

        moris::Tuple< moris::Mat< moris::real>, moris::Mat<moris::real>> myTuple = moris::qr(a);

        moris::Mat<moris::real> R = myTuple.get<1>();

        REQUIRE( R.numel() == 9 );
        REQUIRE( std::abs( R(0,0) - (-14.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(0,1) - ( -21.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(0,2) - (  14.0) ) < 1.0e-12 );
        REQUIRE( std::abs( R(1,0) - (  0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(1,1) - (-175.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(1,2) - (  70.0) ) < 1.0e-12 );
        REQUIRE( std::abs( R(2,0) - (  0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(2,1) - (   0.0) ) < 1.0e-12 ); REQUIRE( std::abs( R(2,2) - ( -35.0) ) < 1.0e-12 );
    }
}
