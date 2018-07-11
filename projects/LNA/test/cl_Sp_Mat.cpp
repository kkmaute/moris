// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"

/**
 * @todo add tests to check size and span, also eigen is not working for .row()
 * function. This need some routine to return sparse matrix's row for a given index.
 */

// ----------------------------------------------------------------------------

TEST_CASE(
                "moris::Sp_Mat::accessors",
                "[linalgebra],[Sp_Mat],[accessors]" )
{
    #include "linalg/cl_Sp_Mat/Sp_Mat.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_cplx.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_col.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_cols.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_n_cols.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_rows.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_n_rows.inc"
    #include "linalg/cl_Sp_Mat/Sp_Mat_numel.inc"

    SECTION( "moris::Sp_Mat::Mat real" )
    {
        REQUIRE( a( 0, 0 ) == 1.0 );
        REQUIRE( a( 1, 0 ) == 3.0 );
        REQUIRE( a( 2, 0 ) == 0.0 );

        REQUIRE( a( 0, 1 ) == 2.0 );
        REQUIRE( a( 1, 1 ) == 4.0 );
        REQUIRE( a( 2, 1 ) == 0.0 );

        REQUIRE( a( 0, 2 ) == 0.0 );
        REQUIRE( a( 1, 2 ) == 5.0 );
        REQUIRE( a( 2, 2 ) == 6.0 );
    }

    SECTION( "moris::Sp_Mat::Mat cplx" )
    {
        std::complex<double> m00 = m( 0, 0 );
        std::complex<double> m01 = m( 0, 1 );
        std::complex<double> m10 = m( 1, 0 );
        std::complex<double> m11 = m( 1, 1 );

        REQUIRE( m00.real() == 0.0  );
        REQUIRE( m01.real() == 6.0  );
        REQUIRE( m01.imag() == 7.0  );
        REQUIRE( m10.real() == 0.3  );
        REQUIRE( m11.imag() == 2.0  );
    }

    SECTION( "moris::Sp_Mat::col" )
    {
        REQUIRE( moris::equal_to( aCol( 0, 0 ) , 1.0 ) );
        REQUIRE( moris::equal_to( aCol( 1, 0 ) , 3.0 ) );
        REQUIRE( moris::equal_to( aCol( 2, 0 ) , 0.0 ) );
    }

    SECTION( "moris::Sp_Mat::col cplx" )
    {
        REQUIRE( moris::equal_to( aColc( 0, 0 ).real() , 0.0 ) );
        REQUIRE( moris::equal_to( aColc( 0, 0 ).imag() , 1.0 ) );
        moris::cplx first={0.0,1.0};
        REQUIRE( moris::equal_to( aColc( 0, 0 ), first ) );

        REQUIRE( moris::equal_to( aColc( 1, 0 ).real() , 6.0 ) );
        REQUIRE( moris::equal_to( aColc( 1, 0 ).imag() , 7.0 ) );
        moris::cplx second={6.0,7.0};
        REQUIRE( moris::equal_to( aColc( 1, 0 ), second ) );

        REQUIRE( moris::equal_to( aColc( 2, 0 ).real() , 2.0 ) );
        REQUIRE( moris::equal_to( aColc( 2, 0 ).imag() , 3.0 ) );
        moris::cplx third={2.0,3.0};
        REQUIRE( moris::equal_to( aColc( 2, 0 ), third ) );
    }

    SECTION( "moris::Sp_Mat::cols" )
    {
        REQUIRE( moris::equal_to( aCols( 0, 0 ) , 2.0 ) );
        REQUIRE( moris::equal_to( aCols( 0, 1 ) , 0.0 ) );
        REQUIRE( moris::equal_to( aCols( 1, 0 ) , 4.0 ) );
        REQUIRE( moris::equal_to( aCols( 1, 1 ) , 5.0 ) );
        REQUIRE( moris::equal_to( aCols( 2, 0 ) , 0.0 ) );
        REQUIRE( moris::equal_to( aCols( 2, 1 ) , 6.0 ) );
    }

    SECTION( "moris::Sp_Mat::n_cols" )
    {
        REQUIRE( moris::equal_to( aNumCols, 3 ) );
        REQUIRE( moris::equal_to( jNumCols, 1 ) );
    }

    SECTION( "moris::Sp_Mat::cols cplx" )
    {
        REQUIRE( moris::equal_to( aColsc( 0, 0 ).real() , 2.0 ) );
        REQUIRE( moris::equal_to( aColsc( 0, 0 ).imag() , 3.0 ) );
        moris::cplx first={2.0,3.0};
        REQUIRE( moris::equal_to( aColsc( 0, 0 ) , first ) );

        REQUIRE( moris::equal_to( aColsc( 1, 0 ).real() , 8.0 ) );
        REQUIRE( moris::equal_to( aColsc( 1, 0 ).imag() , 9.0 ) );
        moris::cplx second={8.0,9.0};
        REQUIRE( moris::equal_to( aColsc( 1, 0 ) , second ) );

        REQUIRE( moris::equal_to( aColsc( 1, 1 ).real() , 0.0 ) );
        REQUIRE( moris::equal_to( aColsc( 1, 1 ).imag() , 1.0 ) );
        moris::cplx third={0.0,1.0};
        REQUIRE( moris::equal_to( aColsc( 1, 1 ) , third ) );
    }

    SECTION( "moris::Sp_Mat::rows" )
    {
        REQUIRE( moris::equal_to( aRows( 0, 0) , 3.0 ) );
        REQUIRE( moris::equal_to( aRows( 0, 1) , 4.0 ) );
        REQUIRE( moris::equal_to( aRows( 0, 2) , 5.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 0) , 0.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 1) , 0.0 ) );
        REQUIRE( moris::equal_to( aRows( 1, 2) , 6.0 ) );
    }

    SECTION( "moris::Sp_Mat::n_rows" )
    {
        REQUIRE( moris::equal_to( aNumRows, 3 ) );
        REQUIRE( moris::equal_to( kNumRows, 1 ) );
    }

    SECTION( "moris::Sp_Mat::rows cplx" )
    {
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ).real() , 6.0 ) );
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ).imag() , 7.0 ) );
        moris::cplx first={6.0,7.0};
        REQUIRE( moris::equal_to( aRowsc( 0, 0 ) , first ) );

        REQUIRE( moris::equal_to( aRowsc( 1, 0 ).real() , 2.0 ) );
        REQUIRE( moris::equal_to( aRowsc( 1, 0 ).imag() , 3.0 ) );
        moris::cplx second={2.0,3.0};
        REQUIRE( moris::equal_to( aRowsc( 1, 0 ) , second ) );

        REQUIRE( moris::equal_to( aRowsc( 1, 1 ).real() , 0.0 ) );
        REQUIRE( moris::equal_to( aRowsc( 1, 1 ).imag() , 5.0 ) );
        moris::cplx third={0.0,5.0};
        REQUIRE( moris::equal_to( aRowsc( 1, 1 ) , third ) );
    }

    SECTION( "moris::Sp_Mat::numel" )
    {
        REQUIRE( moris::equal_to( aNumel, 9 ) );
        REQUIRE( moris::equal_to( kNumel, 1 ) );
    }

    SECTION( "moris::Sp_Mat::get_nnz" )
    {
        #include "linalg/cl_Sp_Mat/Sp_Mat_nnz.inc"

        REQUIRE( NumNonZero == 4 );
    }

    SECTION( "moris::Sp_Mat::clear_sparsity" )
    {
        #include "linalg/cl_Sp_Mat/Sp_Mat_nnz.inc"

        A.clear_sparsity();

        REQUIRE( A.get_nnz() == 0 );
    }
}

TEST_CASE(
                "moris::Sp_Mat::constructors",
                "[linalgebra],[Sp_Mat],[constructors]" )
{
    SECTION( "moris::Sp_Mat::col" )
    {
        // All inputs in column fashion
        #include "linalg/cl_Sp_Mat/cl_Sp_Mat_form1.inc"

        REQUIRE( S1( 2, 0 ) == 1.0 );
        REQUIRE( S1( 1, 1 ) == 2.0 );
        REQUIRE( S1( 2, 1 ) == 3.0 );
        REQUIRE( S1( 3, 2 ) == 4.0 );
        REQUIRE( S1( 0, 3 ) == 5.0 );

        // All inputs in column fashion with user provided n_rows and n_cols
        #include "linalg/cl_Sp_Mat/cl_Sp_Mat_form2.inc"

        REQUIRE( S2( 2, 0 ) == 1.0 );
        REQUIRE( S2( 1, 1 ) == 2.0 );
        REQUIRE( S2( 2, 1 ) == 3.0 );
        REQUIRE( S2( 3, 2 ) == 4.0 );
        REQUIRE( S2( 0, 3 ) == 5.0 );

        // All inputs in column fashion with zero values for sparse matrix
        moris::Mat< moris::real > tValues( 5, 1 );

        tValues( 0, 0 ) = 0.0; tValues( 1, 0 ) = 0.0; tValues( 2, 0 ) = 0.0; tValues( 3, 0 ) = 0.0; tValues( 4, 0 ) = 0.0;

        moris::Sp_Mat< moris::real > S3( RowInd, ColInd, tValues );

        REQUIRE( S3( 2, 0 ) == 0.0 );
        REQUIRE( S3( 1, 1 ) == 0.0 );
        REQUIRE( S3( 2, 1 ) == 0.0 );
        REQUIRE( S3( 3, 2 ) == 0.0 );
        REQUIRE( S3( 0, 3 ) == 0.0 );

        // All inputs in column fashion with user provided n_rows and n_cols and zero values for sparse matrix
        moris::Sp_Mat< moris::real > S4( RowInd, ColInd, tValues, IElems , JElems );

        REQUIRE( S4( 2, 0 ) == 0.0 );
        REQUIRE( S4( 1, 1 ) == 0.0 );
        REQUIRE( S4( 2, 1 ) == 0.0 );
        REQUIRE( S4( 3, 2 ) == 0.0 );
        REQUIRE( S4( 0, 3 ) == 0.0 );

    }

    SECTION( "moris::Sp_Mat copy constructor" )
    {
        #include "linalg/cl_Sp_Mat/cl_Sp_Mat_form1.inc"

        moris::Sp_Mat< moris::real > CopyA( S1 );

        REQUIRE( CopyA( 2, 0 ) == 1.0 );
        REQUIRE( CopyA( 1, 1 ) == 2.0 );
        REQUIRE( CopyA( 2, 1 ) == 3.0 );
        REQUIRE( CopyA( 3, 2 ) == 4.0 );
        REQUIRE( CopyA( 0, 3 ) == 5.0 );
    }

    SECTION( "moris::Sp_Mat::row" )
    {
        // All inputs in row fashion
        moris::Mat< moris::uint > RowInd( 1, 5 );
        moris::Mat< moris::uint > ColInd( 1, 5 );
        moris::Mat< moris::real > Values( 1, 5 );
        moris::uint IElems = 5;
        moris::uint JElems = 5;

        RowInd( 0, 0 ) = 2;   RowInd( 0, 1 ) = 1;   RowInd( 0, 2 ) = 2;   RowInd( 0, 3 ) = 3;   RowInd( 0, 4 ) = 0;
        ColInd( 0, 0 ) = 0;   ColInd( 0, 1 ) = 1;   ColInd( 0, 2 ) = 1;   ColInd( 0, 3 ) = 2;   ColInd( 0, 4 ) = 3;
        Values( 0, 0 ) = 1.0; Values( 0, 1 ) = 2.0; Values( 0, 2 ) = 3.0; Values( 0, 3 ) = 4.0; Values( 0, 4 ) = 5.0;

        moris::Sp_Mat< moris::real > S1( RowInd, ColInd, Values );

        REQUIRE( S1( 2, 0 ) == 1.0 );
        REQUIRE( S1( 1, 1 ) == 2.0 );
        REQUIRE( S1( 2, 1 ) == 3.0 );
        REQUIRE( S1( 3, 2 ) == 4.0 );
        REQUIRE( S1( 0, 3 ) == 5.0 );

        // All inputs in row fashion with user provided n_rows and n_cols
        moris::Sp_Mat< moris::real > S2( RowInd, ColInd, Values, IElems , JElems );

        REQUIRE( S2( 2, 0 ) == 1.0 );
        REQUIRE( S2( 1, 1 ) == 2.0 );
        REQUIRE( S2( 2, 1 ) == 3.0 );
        REQUIRE( S2( 3, 2 ) == 4.0 );
        REQUIRE( S2( 0, 3 ) == 5.0 );

        // All inputs in row fashion with zero values for sparse matrix
        moris::Mat< moris::real > tValues( 1, 5 );

        tValues( 0, 0 ) = 0.0; tValues( 0, 1 ) = 0.0; tValues( 0, 2 ) = 0.0; tValues( 0, 3 ) = 0.0; tValues( 0, 4 ) = 0.0;

        moris::Sp_Mat< moris::real > S3( RowInd, ColInd, tValues );

        REQUIRE( S3( 2, 0 ) == 0.0 );
        REQUIRE( S3( 1, 1 ) == 0.0 );
        REQUIRE( S3( 2, 1 ) == 0.0 );
        REQUIRE( S3( 3, 2 ) == 0.0 );
        REQUIRE( S3( 0, 3 ) == 0.0 );

        // All inputs in row fashion with user provided n_rows and n_cols and zero values for sparse matrix
        moris::Sp_Mat< moris::real > S4( RowInd, ColInd, tValues, IElems , JElems );

        REQUIRE( S4( 2, 0 ) == 0.0 );
        REQUIRE( S4( 1, 1 ) == 0.0 );
        REQUIRE( S4( 2, 1 ) == 0.0 );
        REQUIRE( S4( 3, 2 ) == 0.0 );
        REQUIRE( S4( 0, 3 ) == 0.0 );

    }
}
