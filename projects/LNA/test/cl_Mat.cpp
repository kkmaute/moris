// Third-party header files.
#include <catch.hpp>
#include <iostream>
#include <stdlib.h>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "chronos.hpp"

// ----------------------------------------------------------------------------
using namespace moris;
TEST_CASE(
                "moris::Mat",
                "[linalgebra],[Mat]" )
{
    #include "linalg/cl_Mat/Mat.inc"
    #include "linalg/cl_Mat/Mat_cplx.inc"
    #include "linalg/cl_Mat/Mat_fill.inc"
    #include "linalg/cl_Mat/Mat_col.inc"
    #include "linalg/cl_Mat/Mat_cols.inc"
    #include "linalg/cl_Mat/Mat_n_cols.inc"
    #include "linalg/cl_Mat/Mat_row.inc"
    #include "linalg/cl_Mat/Mat_rows.inc"
    #include "linalg/cl_Mat/Mat_n_rows.inc"
    #include "linalg/cl_Mat/Mat_numel.inc"
    #include "linalg/cl_Mat/Mat_size.inc"
    #include "linalg/cl_Mat/Mat_length.inc"
    #include "linalg/cl_Mat/Mat_span.inc"
    #include "linalg/cl_Mat/Mat_min.inc"
    #include "linalg/cl_Mat/Mat_max.inc"


    SECTION( "moris::Mat random checks" )
    {
        REQUIRE( moris::equal_to( a( 0 ), 1.0 ) );
        REQUIRE( moris::equal_to( a( 1 ), 4.0 ) );
        REQUIRE( moris::equal_to( a( 2 ), 9.0 ) );

        REQUIRE( moris::equal_to( a( 6 ), 3.0 ) );
        REQUIRE( moris::equal_to( a( 7 ), 6.0 ) );
        REQUIRE( moris::equal_to( a( 8 ), 9.0 ) );

        REQUIRE_THROWS( a( 9 ) );

        REQUIRE( moris::equal_to( a( 0, 0 ), 1.0 ) );
        REQUIRE( moris::equal_to( a( 1, 0 ), 4.0 ) );
        REQUIRE( moris::equal_to( a( 2, 0 ), 9.0 ) );

        REQUIRE( ( m( 0, 0 ).real(), 1.0 ) );
        REQUIRE( ( m( 0, 1 ).imag(), 7.0 ) );
        REQUIRE( ( m( 0, 1 ).real(), 6.0 ) );
        REQUIRE( ( m( 1, 0 ).real(), 0.3 ) );
        REQUIRE( ( m( 1, 1 ).imag(), 2.0 ) );
    }

    SECTION( "moris::Mat copy constructor" )
    {
        moris::Mat< moris::real > CopyA( a );

        REQUIRE( moris::equal_to( CopyA( 0 ), 1.0 ) );
        REQUIRE( moris::equal_to( CopyA( 1 ), 4.0 ) );
        REQUIRE( moris::equal_to( CopyA( 2 ), 9.0 ) );

        REQUIRE( moris::equal_to( CopyA( 6 ), 3.0 ) );
        REQUIRE( moris::equal_to( CopyA( 7 ), 6.0 ) );
        REQUIRE( moris::equal_to( CopyA( 8 ), 9.0 ) );
    }

    SECTION( "moris::Mat::set_size and fill constructor" )
    {
        REQUIRE( rMat( 0, 0 ) == 123.5 );
        REQUIRE( rMat( 1, 0 ) == 123.5 );
        REQUIRE( rMat( 2, 0 ) == 123.5 );

        REQUIRE( rMat( 0, 1 ) == 123.5 );
        REQUIRE( rMat( 1, 1 ) == 123.5 );
        REQUIRE( rMat( 2, 1 ) == 123.5 );

        REQUIRE( rMat( 0, 2 ) == 123.5 );
        REQUIRE( rMat( 1, 2 ) == 123.5 );
        REQUIRE( rMat( 2, 2 ) == 123.5 );

        REQUIRE( iMat( 0, 0 ) == 10 );
        REQUIRE( iMat( 1, 0 ) == 10 );
        REQUIRE( iMat( 2, 0 ) == 10 );

        REQUIRE( iMat( 0, 1 ) == 10 );
        REQUIRE( iMat( 1, 1 ) == 10 );
        REQUIRE( iMat( 2, 1 ) == 10 );

        REQUIRE( iMat( 0, 2 ) == 10 );
        REQUIRE( iMat( 1, 2 ) == 10 );
        REQUIRE( iMat( 2, 2 ) == 10 );
    }

    SECTION( "moris::Mat::pointer constructor" )
    {
        moris::real* PointerA = new moris::real[6];

        PointerA[0] = 12.3;
        PointerA[1] = 2.35;
        PointerA[2] = 10.0;
        PointerA[3] = 14.1;
        PointerA[4] = -0.3;
        PointerA[5] = -0.3;

        moris::Mat< moris::real > MatA = moris::Mat< moris::real > (PointerA, 3, 2);

        REQUIRE( std::abs( PointerA[0] - MatA( 0, 0 ) ) < 1.0e-12 );
        REQUIRE( std::abs( PointerA[1] - 2.35 ) < 1.0e-12 );
        REQUIRE( std::abs( PointerA[3] - MatA( 0, 1 ) ) < 1.0e-12 );
        REQUIRE( std::abs( PointerA[5] - MatA( 2, 1 ) ) < 1.0e-12 );

        delete PointerA;
    }

    SECTION( "moris::Mat::fill real" )
    {
        REQUIRE( aMatr( 0, 0 ) == 123.5 );
        REQUIRE( aMatr( 1, 0 ) == 123.5 );
        REQUIRE( aMatr( 2, 0 ) == 123.5 );

        REQUIRE( aMatr( 0, 1 ) == 123.5 );
        REQUIRE( aMatr( 1, 1 ) == 123.5 );
        REQUIRE( aMatr( 2, 1 ) == 123.5 );

        REQUIRE( aMatr( 0, 2 ) == 123.5 );
        REQUIRE( aMatr( 1, 2 ) == 123.5 );
        REQUIRE( aMatr( 2, 2 ) == 123.5 );
    }

    SECTION( "moris::Mat::fill cplx" )
    {
        std::complex<double> aMatc_00 = aMatc( 0, 0 );
        std::complex<double> aMatc_10 = aMatc( 1, 0 );
        std::complex<double> aMatc_20 = aMatc( 2, 0 );

        std::complex<double> aMatc_01 = aMatc( 0, 1 );
        std::complex<double> aMatc_11 = aMatc( 1, 1 );
        std::complex<double> aMatc_21 = aMatc( 2, 1 );

        std::complex<double> aMatc_02 = aMatc( 0, 2 );
        std::complex<double> aMatc_12 = aMatc( 1, 2 );
        std::complex<double> aMatc_22 = aMatc( 2, 2 );

        REQUIRE( aMatc_00.real() == 123.5 );
        REQUIRE( aMatc_10.real() == 123.5 );
        REQUIRE( aMatc_20.real() == 123.5 );

        REQUIRE( aMatc_01.real() == 123.5 );
        REQUIRE( aMatc_11.real() == 123.5 );
        REQUIRE( aMatc_21.real() == 123.5 );

        REQUIRE( aMatc_02.real() == 123.5 );
        REQUIRE( aMatc_12.real() == 123.5 );
        REQUIRE( aMatc_22.real() == 123.5 );

        REQUIRE( aMatc_00.imag() == 1.0 );
        REQUIRE( aMatc_10.imag() == 1.0 );
        REQUIRE( aMatc_20.imag() == 1.0 );

        REQUIRE( aMatc_01.imag() == 1.0 );
        REQUIRE( aMatc_11.imag() == 1.0 );
        REQUIRE( aMatc_21.imag() == 1.0 );

        REQUIRE( aMatc_02.imag() == 1.0 );
        REQUIRE( aMatc_12.imag() == 1.0 );
        REQUIRE( aMatc_22.imag() == 1.0 );
    }

    SECTION( "moris::Mat::fill uint" )
    {
        REQUIRE( aMati( 0, 0 ) == 123 );
        REQUIRE( aMati( 1, 0 ) == 123 );
        REQUIRE( aMati( 2, 0 ) == 123 );

        REQUIRE( aMati( 0, 1 ) == 123 );
        REQUIRE( aMati( 1, 1 ) == 123 );
        REQUIRE( aMati( 2, 1 ) == 123 );

        REQUIRE( aMati( 0, 2 ) == 123 );
        REQUIRE( aMati( 1, 2 ) == 123 );
        REQUIRE( aMati( 2, 2 ) == 123 );
    }

    SECTION( "moris::Mat::col" )
    {
        REQUIRE( moris::equal_to( aCol( 0, 0 ) , 1.0 ) );
        REQUIRE( moris::equal_to( aCol( 1, 0 ) , 4.0 ) );
        REQUIRE( moris::equal_to( aCol( 2, 0 ) , 9.0 ) );

        a.col(0) += a.col(0);

        REQUIRE( moris::equal_to( a( 0, 0 ) ,  2.0 ) );
        REQUIRE( moris::equal_to( a( 1, 0 ) ,  8.0 ) );
        REQUIRE( moris::equal_to( a( 2, 0 ) , 18.0 ) );

    }

    SECTION( "moris::Mat::col cplx" )
    {
        REQUIRE( moris::equal_to( aColc( 0, 0 ).real() , 0.0 ) );  //check of real part
        REQUIRE( moris::equal_to( aColc( 0, 0 ).imag() , 1.0 ) );  //check of imag part
        moris::cplx first={0.0,1.0};
        REQUIRE( moris::equal_to( aColc( 0, 0 ), first ) );     //check of cplx number as a whole

        REQUIRE( moris::equal_to( aColc( 1, 0 ).real() , 6.0 ) );
        REQUIRE( moris::equal_to( aColc( 1, 0 ).imag() , 7.0 ) );
        moris::cplx second={6.0,7.0};
        REQUIRE( moris::equal_to( aColc( 1, 0 ), second ) );

        REQUIRE( moris::equal_to( aColc( 2, 0 ).real() , 2.0 ) );
        REQUIRE( moris::equal_to( aColc( 2, 0 ).imag() , 3.0 ) );
        moris::cplx third={2.0,3.0};
        REQUIRE( moris::equal_to( aColc( 2, 0 ), third ) );
    }


    SECTION( "moris::Mat::cols" )
    {
        REQUIRE( moris::equal_to( aCols( 0, 0 ) , 2.0 ) );
        REQUIRE( moris::equal_to( aCols( 0, 1 ) , 3.0 ) );
        REQUIRE( moris::equal_to( aCols( 1, 0 ) , 5.0 ) );
        REQUIRE( moris::equal_to( aCols( 1, 1 ) , 6.0 ) );
        REQUIRE( moris::equal_to( aCols( 2, 0 ) , 8.0 ) );
        REQUIRE( moris::equal_to( aCols( 2, 1 ) , 9.0 ) );

        a.cols(1,2) += 2*a.cols(1,2);

        REQUIRE( moris::equal_to( a( 0, 1 ) , 6.0 ) );
        REQUIRE( moris::equal_to( a( 0, 2 ) , 9.0 ) );
        REQUIRE( moris::equal_to( a( 1, 1 ) , 15.0 ) );
        REQUIRE( moris::equal_to( a( 1, 2 ) , 18.0 ) );
        REQUIRE( moris::equal_to( a( 2, 1 ) , 24.0 ) );
        REQUIRE( moris::equal_to( a( 2, 2 ) , 27.0 ) );
    }


    SECTION( "moris::Mat::n_cols" )
    {
        REQUIRE( moris::equal_to( aNumCols, 3 ) );
        REQUIRE( moris::equal_to( jNumCols, 1 ) );
    }

    SECTION( "moris::Mat::cols cplx" )
    {
        REQUIRE( moris::equal_to( aColsc( 0, 0 ).real() , 2.0 ) );  //check of real part
        REQUIRE( moris::equal_to( aColsc( 0, 0 ).imag() , 3.0 ) );  //check of imag part
        moris::cplx first={2.0,3.0};
        REQUIRE( moris::equal_to( aColsc( 0, 0 ) , first ) );     //check of cplx number as a whole

        REQUIRE( moris::equal_to( aColsc( 1, 0 ).real() , 8.0 ) );
        REQUIRE( moris::equal_to( aColsc( 1, 0 ).imag() , 9.0 ) );
        moris::cplx second={8.0,9.0};
        REQUIRE( moris::equal_to( aColsc( 1, 0 ) , second ) );

        REQUIRE( moris::equal_to( aColsc( 1, 1 ).real() , 0.0 ) );
        REQUIRE( moris::equal_to( aColsc( 1, 1 ).imag() , 1.0 ) );
        moris::cplx third={0.0,1.0};
        REQUIRE( moris::equal_to( aColsc( 1, 1 ) , third ) );
    }

    SECTION( "moris::Mat::row" )
    {
        REQUIRE( moris::equal_to( aRow( 0, 0 ), 4.0 ) );
        REQUIRE( moris::equal_to( aRow( 0, 1 ), 5.0 ) );
        REQUIRE( moris::equal_to( aRow( 0, 2 ), 6.0 ) );
    }

    SECTION( "moris::Mat::row cplx" )
    {
        REQUIRE( moris::equal_to( aRowc( 0, 0 ).real() , 6.0 ) );  //check of real part
        REQUIRE( moris::equal_to( aRowc( 0, 0 ).imag() , 7.0 ) );  //check of imag part
        moris::cplx first={6.0,7.0};
        REQUIRE( moris::equal_to( aRowc( 0, 0 ), first ) );     //check of cplx number as a whole

        REQUIRE( moris::equal_to( aRowc( 0, 1 ).real() , 8.0 ) );
        REQUIRE( moris::equal_to( aRowc( 0, 1 ).imag() , 9.0 ) );
        moris::cplx second={8.0,9.0};
        REQUIRE( moris::equal_to( aRowc( 0, 1 ), second ) );

        REQUIRE( moris::equal_to( aRowc( 0, 2 ).real() , 0.0 ) );
        REQUIRE( moris::equal_to( aRowc( 0, 2 ).imag() , 1.0 ) );
        moris::cplx third={0.0,1.0};        moris::Mat< moris::cplx > aSpanc( 2, 2 );
        aSpanc = c( {1, 2}, {1, 2} );
        REQUIRE( moris::equal_to( aRowc( 0, 2 ), third ) );

    }

    SECTION( "moris::Mat::rows" )
    {
        REQUIRE( moris::equal_to( aRows( 0, 0) , 4.0 ) );
        REQUIRE( moris::equal_to( aRows( 0, 1) , 5.0 ) );
        REQUIRE( moris::equal_to( aRows( 0, 2) , 6.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 0) , 9.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 1) , 8.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 2) , 9.0 ) );
    }


    SECTION( "moris::Mat::n_rows" )
    {
        REQUIRE( moris::equal_to( aNumRows, 3 ) );
        REQUIRE( moris::equal_to( kNumRows, 1 ) );
    }

    SECTION( "moris::Mat::rows cplx" )
    {
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ).real() , 6.0 ) );  //check of real part
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ).imag() , 7.0 ) );  //check of imag part
        moris::cplx first={6.0,7.0};
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ) , first ) );     //check of cplx number as a whole

        REQUIRE( moris::equal_to( aRowsc( 1, 0 ).real() , 2.0 ) );
        REQUIRE( moris::equal_to( aRowsc( 1, 0 ).imag() , 3.0 ) );
        moris::cplx second={2.0,3.0};
        REQUIRE( moris::equal_to( aRowsc( 1, 0 ) , second ) );

        REQUIRE( moris::equal_to( aRowsc( 1, 1 ).real() , 4.0 ) );
        REQUIRE( moris::equal_to( aRowsc( 1, 1 ).imag() , 5.0 ) );
        moris::cplx third={4.0,5.0};
        REQUIRE( moris::equal_to( aRowsc( 1, 1 ) , third ) );
    }

    SECTION( "moris::Mat::numel" )
    {
        REQUIRE( moris::equal_to( aNumel, 9 ) );
        REQUIRE( moris::equal_to( kNumel, 1 ) );
    }

    SECTION( "moris::Mat::size" )
    {
        moris::Tuple< moris::size_t, moris::size_t > tuplea( 3, 4 );
        moris::Tuple< moris::size_t, moris::size_t > tupleb( 3, 3 );
        REQUIRE_FALSE(aSize == tuplea );
        REQUIRE(aSize == tupleb );
        REQUIRE( moris::equal_to( nSizeRows, 5 ) );
        REQUIRE( moris::equal_to( nSizeCols, 7 ) );
    }

    SECTION( "moris::Mat::length" )
    {
        REQUIRE( moris::equal_to( tRowVecLength, 4 ) );
        REQUIRE( moris::equal_to( tColVecLength, 7 ) );
    }

    SECTION( "moris::Mat::span" )
    {
        REQUIRE( moris::equal_to( aSpan( 0, 0 ), 5.0 ) );
        REQUIRE( moris::equal_to( aSpan( 0, 1 ), 6.0 ) );
        REQUIRE( moris::equal_to( aSpan( 1, 0 ), 8.0 ) );
        REQUIRE( moris::equal_to( aSpan( 1, 1 ), 9.0 ) );

        moris::Mat< moris::cplx > aSpanc( 2, 2 );
        aSpanc = c( {1, 2}, {1, 2} );
    }

    SECTION( "moris::Mat::span cplx" )
    {
        // Check of real part.
        REQUIRE( moris::equal_to( aSpanc( 0, 0 ).real(), 8.0 ) );

        // Check of imag part.
        REQUIRE( moris::equal_to( aSpanc( 0, 0 ).imag(), 9.0 ) );

        moris::cplx first = {8.0, 9.0};

        // Check of cplx number as a whole.
        REQUIRE( moris::equal_to( aSpanc( 0, 0 ), first ) );
        REQUIRE( moris::equal_to( aSpanc( 1, 0 ).real(), 4.0 ) );
        REQUIRE( moris::equal_to( aSpanc( 1, 0 ).imag(), 5.0 ) );

        moris::cplx second = {4.0, 5.0};
        REQUIRE( moris::equal_to( aSpanc( 1, 0 ), second ) );
        REQUIRE( moris::equal_to( aSpanc( 1, 1 ).real(), 6.0 ) );
        REQUIRE( moris::equal_to( aSpanc( 1, 1 ).imag(), 7.0 ) );

        moris::cplx third = {6.0, 7.0};
        REQUIRE( moris::equal_to( aSpanc( 1, 1 ), third ) );
    }

    SECTION( "moris::Mat::eye")
    {
        moris::Mat< moris::real > a(4,2);

        REQUIRE( a.size(0) == 4 ); REQUIRE( a.size(1) == 2);

        a.eye(3);

        REQUIRE( a.size(0) == 3 ); REQUIRE( a.size(1) == 3);
        REQUIRE( moris::equal_to( a(0,0), 1.0) ); REQUIRE( moris::equal_to( a(0,1), 0.0) ); REQUIRE( moris::equal_to( a(0,2), 0.0) );
        REQUIRE( moris::equal_to( a(1,0), 0.0) ); REQUIRE( moris::equal_to( a(1,1), 1.0) ); REQUIRE( moris::equal_to( a(1,2), 0.0) );
        REQUIRE( moris::equal_to( a(2,0), 0.0) ); REQUIRE( moris::equal_to( a(2,1), 0.0) ); REQUIRE( moris::equal_to( a(2,2), 1.0) );

        moris::Mat< moris::real > c;
        moris::Mat< moris::real > b = c.eye(3);

        REQUIRE( b.size(0) == 3 ); REQUIRE( b.size(1) == 3);
        REQUIRE( moris::equal_to( b(0,0), 1.0) ); REQUIRE( moris::equal_to( b(0,1), 0.0) ); REQUIRE( moris::equal_to( b(0,2), 0.0) );
        REQUIRE( moris::equal_to( b(1,0), 0.0) ); REQUIRE( moris::equal_to( b(1,1), 1.0) ); REQUIRE( moris::equal_to( b(1,2), 0.0) );
        REQUIRE( moris::equal_to( b(2,0), 0.0) ); REQUIRE( moris::equal_to( b(2,1), 0.0) ); REQUIRE( moris::equal_to( b(2,2), 1.0) );

        REQUIRE( c.size(0) == 3 ); REQUIRE( c.size(1) == 3);
        REQUIRE( moris::equal_to( c(0,0), 1.0) ); REQUIRE( moris::equal_to( c(0,1), 0.0) ); REQUIRE( moris::equal_to( c(0,2), 0.0) );
        REQUIRE( moris::equal_to( c(1,0), 0.0) ); REQUIRE( moris::equal_to( c(1,1), 1.0) ); REQUIRE( moris::equal_to( c(1,2), 0.0) );
        REQUIRE( moris::equal_to( c(2,0), 0.0) ); REQUIRE( moris::equal_to( c(2,1), 0.0) ); REQUIRE( moris::equal_to( c(2,2), 1.0) );
    }

    SECTION( "moris::Mat::set_size" )
    {
        moris::Mat< moris::real > b( 3, 3 );
        b( 0, 0 ) = 1.0; b( 0, 1 ) = 2.0; b( 0, 2 ) = 3.0;
        b( 1, 0 ) = 4.0; b( 1, 1 ) = 5.0; b( 1, 2 ) = 6.0;
        b( 2, 0 ) = 7.0; b( 2, 1 ) = 8.0; b( 2, 2 ) = 9.0;

        b.set_size(4 , 4);
        REQUIRE( b.numel() == 16 );
    }

    SECTION( "complex moris::Mat::set_size" )
    {
        moris::Mat< moris::cplx > A( 2, 2 );

        A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0};
        A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0};

        A.resize(2 , 2);
        REQUIRE( A.numel() == 4 );
    }

    SECTION( "moris::Mat::resize",
     "[make larger]" ) // make larger
    {
        moris::Mat< moris::real > d( 3, 3 );
        d( 0, 0 ) = 1.0; d( 0, 1 ) = 2.0; d( 0, 2 ) = 3.0;
        d( 1, 0 ) = 4.0; d( 1, 1 ) = 5.0; d( 1, 2 ) = 6.0;
        d( 2, 0 ) = 7.0; d( 2, 1 ) = 8.0; d( 2, 2 ) = 9.0;

        d.resize(4 , 4);
        REQUIRE( d.numel() == 16 );
        REQUIRE( d(0,0)   == 1.0 ); REQUIRE( d(0,1)   == 2.0 ); REQUIRE( d(0,2)   == 3.0 );
        REQUIRE( d(1,0)   == 4.0 ); REQUIRE( d(1,1)   == 5.0 ); REQUIRE( d(1,2)   == 6.0 );
        REQUIRE( d(2,0)   == 7.0 ); REQUIRE( d(2,1)   == 8.0 ); REQUIRE( d(2,2)   == 9.0 );
    }

    SECTION( "moris::Mat::resize",
             "[make smaller]" )
    {
        moris::Mat< moris::real > d( 3, 3 );
        d( 0, 0 ) = 1.0; d( 0, 1 ) = 2.0; d( 0, 2 ) = 3.0;
        d( 1, 0 ) = 4.0; d( 1, 1 ) = 5.0; d( 1, 2 ) = 6.0;
        d( 2, 0 ) = 7.0; d( 2, 1 ) = 8.0; d( 2, 2 ) = 9.0;

        d.resize(2 , 2);
        REQUIRE( d.numel() == 4 );
        REQUIRE( d(0,0)   == 1.0 ); REQUIRE( d(0,1)   == 2.0 );
        REQUIRE( d(1,0)   == 4.0 ); REQUIRE( d(1,1)   == 5.0 );
    }

    SECTION( "moris::Mat::resize",
             "[make complex smaller]")
    {
        moris::Mat< moris::cplx > A( 2, 2 );

        A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0};
        A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0};

        A.resize(4 , 3);
        REQUIRE( A.numel() == 12 );
        REQUIRE( moris::equal_to( A( 0,0 ).real(),  3.0 ) );
        REQUIRE( moris::equal_to( A( 0,0 ).imag(),  2.0 ) );

        REQUIRE( moris::equal_to( A( 0,1 ).real(),  4.0 ) );
        REQUIRE( moris::equal_to( A( 0,1 ).imag(),  5.0 ) );

        REQUIRE( moris::equal_to( A( 1,0 ).real(),  2.0 ) );
        REQUIRE( moris::equal_to( A( 1,0 ).imag(), -1.0 ) );

        REQUIRE( moris::equal_to( A( 1,1 ).real(), -6.0 ) );
        REQUIRE( moris::equal_to( A( 1,1 ).imag(), -3.0 ) );
    }

    SECTION( "moris::Mat::resize",
             "[make complex smaller]")
    {
        moris::Mat< moris::cplx > A( 3, 3 );

        A( 0, 0 ) = { 3.0,  2.0}; A( 0, 1 ) = { 4.0,  5.0}; A( 0, 2 ) = { 3.0, -4.0};
        A( 1, 0 ) = { 2.0, -1.0}; A( 1, 1 ) = {-6.0, -3.0}; A( 1, 2 ) = {-1.0,  3.0};
        A( 2, 0 ) = { 1.0, -3.0}; A( 2, 1 ) = {-5.0, -2.0}; A( 2, 2 ) = {-6.0,  2.0};

        A.resize(4 , 3);
        REQUIRE( A.numel() == 12 );
        REQUIRE( moris::equal_to( A( 0,0 ).real(),  3.0 ) );
        REQUIRE( moris::equal_to( A( 0,0 ).imag(),  2.0 ) );

        REQUIRE( moris::equal_to( A( 0,1 ).real(),  4.0 ) );
        REQUIRE( moris::equal_to( A( 0,1 ).imag(),  5.0 ) );

        REQUIRE( moris::equal_to( A( 1,0 ).real(),  2.0 ) );
        REQUIRE( moris::equal_to( A( 1,0 ).imag(), -1.0 ) );

        REQUIRE( moris::equal_to( A( 1,1 ).real(), -6.0 ) );
        REQUIRE( moris::equal_to( A( 1,1 ).imag(), -3.0 ) );
    }

    SECTION( "moris::Mat::copy_size" )
    {
        moris::Mat< moris::real > e( 2, 2 );
        moris::Mat< moris::real > f;
        moris::Mat< moris::real > g( 4, 5 );

        e.copy_size(g);
        f.copy_size(g);

        REQUIRE( e.numel() == 20 );
        REQUIRE( f.numel() == 20 );
    }

    SECTION( "complex moris::Mat::copy_size" )
    {

        moris::Mat< moris::cplx > p( 1, 1 );
        moris::Mat< moris::cplx > q( 2, 3 );

        REQUIRE( p.numel() == 1 );

        p.copy_size(q);

        REQUIRE( q.numel() == 6 );
    }

    SECTION( "append matrix")
    {
        moris::Mat< moris::real> A;
        A.eye(3);

        moris::Mat< moris::real> B;
        B.eye(3);

        A.resize(6,3); // needs to be resized first!
        A.rows(3,5) = B.rows(0,2);

        REQUIRE( moris::equal_to( A(0,0), 1.0 )); REQUIRE( moris::equal_to( A(0,1), 0.0 )); REQUIRE( moris::equal_to( A(0,2), 0.0 ));
        REQUIRE( moris::equal_to( A(2,0), 0.0 )); REQUIRE( moris::equal_to( A(2,1), 0.0 )); REQUIRE( moris::equal_to( A(2,2), 1.0 ));

        REQUIRE( moris::equal_to( A(3,0), 1.0 )); REQUIRE( moris::equal_to( A(3,1), 0.0 )); REQUIRE( moris::equal_to( A(3,2), 0.0 ));
        REQUIRE( moris::equal_to( A(5,0), 0.0 )); REQUIRE( moris::equal_to( A(5,1), 0.0 )); REQUIRE( moris::equal_to( A(5,2), 1.0 ));
    }

    SECTION( "real moris::Mat::min" )
    {
        moris::uint RowIndex = 0;
        moris::uint ColIndex = 0;

        REQUIRE( A.min() == 1.0 );
        REQUIRE( A.min( RowIndex, ColIndex ) == 1.0 );
        REQUIRE( RowIndex == 0 );
        REQUIRE( ColIndex == 0 );
    }

    SECTION( "real moris::Mat::mint1" )
    {
        moris::uint RowIndex;
        moris::uint ColIndex;
        moris::real MinVal;

        moris::tie(MinVal, RowIndex, ColIndex) = A.mint();

        REQUIRE( MinVal == 1.0 );
        REQUIRE( RowIndex == 0 );
        REQUIRE( ColIndex == 0 );
    }

    SECTION( "real moris::Mat::mint2" )
    {
        moris::uint RowIndex = 1;
        moris::uint ColIndex = 1;
        moris::real MinVal;

        moris::tie(MinVal, RowIndex, moris::tilde) = A.mint();

        REQUIRE( MinVal == 1.0 );
        REQUIRE( RowIndex == 0 );
        REQUIRE( ColIndex == 1 );
    }

    SECTION( "real moris::Mat::max" )
    {
        moris::uint RowIndex = 0;
        moris::uint ColIndex = 0;

        REQUIRE( B.max() == 9.0 );
        REQUIRE( B.max( RowIndex, ColIndex ) == 9.0 );
        REQUIRE( RowIndex == 2 );
        REQUIRE( ColIndex == 0 );
    }

    SECTION( "real moris::Mat::maxt" )
    {
        moris::uint rowInd = 0;
        moris::uint colInd = 0;
        moris::real maxVal;

        moris::tie(maxVal, rowInd, colInd) = B.maxt();

        REQUIRE( maxVal == 9.0 );
        REQUIRE( rowInd == 2 );
        REQUIRE( colInd == 0 );
    }

}


TEST_CASE(
        "moris::MatPerformance",
        "[linalgebra],[MatPerformance]" )
{
    moris::uint nmats = 3;
    moris::uint nruns = 10000;
    moris::Mat< moris::real > TimeMat( nmats, 4 );

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        moris::uint nrows = 10 + n * nmats;
        moris::uint ncols = 10 + n * nmats;
        moris::Mat< moris::real > A( nrows, ncols );
        moris::Mat< moris::real > B( nrows, ncols );
        moris::Mat< moris::real > C( nrows, ncols );
        moris::Mat< moris::real > D( nrows, ncols );
        moris::Mat< moris::real > E( nrows, ncols );
        moris::Mat< moris::real > F( nrows, ncols );
        moris::Mat< moris::real > G( nrows, ncols );
        moris::Mat< moris::real > H( nrows, ncols );
        moris::Mat< moris::real > H1( nrows, ncols );
        moris::Mat< moris::real > H2( nrows, ncols );

        for ( moris::uint i = 0; i < nrows; i++ )
        {
            for ( moris::uint j = 0; j < ncols; j++ )
            {
                A( i, j ) = i * j;
                B( i, j ) = i * j;
                C( i, j ) = i * j;
                D( i, j ) = i * j;
            }
        }

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            E = A + B + C + D;
        }

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            F = A - B - C - D;
        }

        moris::real wall_time_microseconds_tF = tF.toc<moris::chronos::microseconds>().wall;

        // operator multiplication
        moris::tic tG;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            G = A * B * C * D;
        }

        moris::real wall_time_microseconds_tG = tG.toc<moris::chronos::microseconds>().wall;

        // operator elementwise division
        moris::tic tH;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            H1 = A / B;
            H2 = H1 / C;
            H  = H2 / D;
        }

        moris::real wall_time_microseconds_tH = tH.toc<moris::chronos::microseconds>().wall;

        TimeMat( n, 0 ) = wall_time_microseconds_tE;
        TimeMat( n, 1 ) = wall_time_microseconds_tF;
        TimeMat( n, 2 ) = wall_time_microseconds_tG;
        TimeMat( n, 3 ) = wall_time_microseconds_tH;
    }

#ifdef MORIS_USE_ARMA
//    std::cout << "The results based of moris_arma_mat:" << std::endl;
//    std::cout << TimeMat.data() << std::endl;

    // first matrix
    REQUIRE( TimeMat( 0, 0 ) < 8.00e04 );
    REQUIRE( TimeMat( 0, 1 ) < 7.00e04 );
    REQUIRE( TimeMat( 0, 2 ) < 3.00e05 );
    REQUIRE( TimeMat( 0, 3 ) < 7.00e04 );

    // second matrix
    REQUIRE( TimeMat( 1, 0 ) < 1.00e05 );
    REQUIRE( TimeMat( 1, 1 ) < 1.00e05 );
    REQUIRE( TimeMat( 1, 2 ) < 7.00e05 );
    REQUIRE( TimeMat( 1, 3 ) < 1.00e05 );

    // third matrix
    REQUIRE( TimeMat( 2, 0 ) < 2.00e05 );
    REQUIRE( TimeMat( 2, 1 ) < 2.00e05 );
    REQUIRE( TimeMat( 2, 2 ) < 1.00e06 );
    REQUIRE( TimeMat( 2, 3 ) < 2.00e05 );

#elif  MORIS_USE_EIGEN
//    std::cout << "The results based of moris_eigen_mat:" << std::endl;
//    std::cout << TimeMat.data() << std::endl;

    // first matrix
    REQUIRE( TimeMat( 0, 0 ) < 1.00e06 );
    REQUIRE( TimeMat( 0, 1 ) < 1.00e06 );
    REQUIRE( TimeMat( 0, 2 ) < 1.00e06 );
    REQUIRE( TimeMat( 0, 3 ) < 2.00e05 );

    // second matrix
    REQUIRE( TimeMat( 1, 0 ) < 2.00e05 );
    REQUIRE( TimeMat( 1, 1 ) < 2.00e05 );
    REQUIRE( TimeMat( 1, 2 ) < 2.00e06 );
    REQUIRE( TimeMat( 1, 3 ) < 3.00e05 );

    // third matrixcredit.
    REQUIRE( TimeMat( 2, 0 ) < 3.00e05 );
    REQUIRE( TimeMat( 2, 1 ) < 3.00e05 );
    REQUIRE( TimeMat( 2, 2 ) < 2.50e06 );
    REQUIRE( TimeMat( 2, 3 ) < 4.00e05 );

#endif

}
