/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include "UT_LINALG_Test_Helpers.hpp"
#include <algorithm>
#include <numeric>
#include <sstream>
#include <type_traits>
#include "op_move.hpp"
#include "fn_trans.hpp"
#include "op_plus.hpp"
#include "op_equal_equal.hpp"
#include "op_not_equal.hpp"
#include "op_less.hpp"
#include "op_less_equal.hpp"
#include "op_greater.hpp"
#include "op_greater_equal.hpp"

namespace moris
{
    using linalg_test::check_matrix;

    TEST_CASE( "LINALG fixed and dynamic matrices support value conversion", "[linalgebra],[matrix_api]" )
    {
        Matrix< F31RMat > tFixed   = { { 1.0 }, { 2.0 }, { 3.0 } };
        Matrix< DDRMat >  tDynamic = tFixed;
        check_matrix( tDynamic, Matrix< DDRMat >{ { 1.0 }, { 2.0 }, { 3.0 } } );
        tFixed( 0 ) = 4.0;
        CHECK( tDynamic( 0 ) == 1.0 );

        // SDF passes fixed coordinates to functions accepting a dynamic matrix reference.
        const auto tSum = []( const Matrix< DDRMat >& aCoordinates ) {
            return aCoordinates( 0 ) + aCoordinates( 1 ) + aCoordinates( 2 );
        };
        CHECK( tSum( tFixed ) == 9.0 );
        tDynamic = tFixed;
        CHECK( tDynamic( 0 ) == 4.0 );
        tDynamic( 2 )       = 6.0;
        auto tDifference    = tFixed - tDynamic;
        auto tSumExpression = tDynamic + tFixed;
        tFixed( 1 )         = 5.0;
        check_matrix( Matrix< DDRMat >( tDifference ), Matrix< DDRMat >{ { 0.0 }, { 3.0 }, { -3.0 } } );
        check_matrix( Matrix< DDRMat >( tSumExpression ), Matrix< DDRMat >{ { 8.0 }, { 7.0 }, { 9.0 } } );
        check_matrix( Matrix< DDRMat >( tDynamic - tFixed ), Matrix< DDRMat >{ { 0.0 }, { -3.0 }, { 3.0 } } );
        check_matrix( Matrix< DDRMat >( tFixed + tDynamic ), Matrix< DDRMat >{ { 8.0 }, { 7.0 }, { 9.0 } } );
        tFixed = tDynamic;
        check_matrix( tFixed, Matrix< DDRMat >{ { 4.0 }, { 2.0 }, { 6.0 } } );
    }

    TEST_CASE( "LINALG scalar multiplication remains lazy in both operand orders", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat >        tA     = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat >& tConst = tA;
        auto                    tRight = tConst * 2.0;
        auto                    tLeft  = 3.0 * tConst;
        tA( 0, 0 )                     = 5.0;
        check_matrix( Matrix< DDRMat >( tRight ), Matrix< DDRMat >{ { 10.0, 4.0 }, { 6.0, 8.0 } } );
        check_matrix( Matrix< DDRMat >( tLeft ), Matrix< DDRMat >{ { 15.0, 6.0 }, { 9.0, 12.0 } } );

        Matrix< DDCMat > tComplex = { { cplx( 1.0, 2.0 ), cplx( 3.0, -1.0 ) } };
        check_matrix( Matrix< DDCMat >( tComplex * cplx( 0.0, 1.0 ) ),
                Matrix< DDCMat >{ { cplx( -2.0, 1.0 ), cplx( 1.0, 3.0 ) } } );
        check_matrix( Matrix< DDCMat >( tComplex + cplx( 1.0, -1.0 ) ),
                Matrix< DDCMat >{ { cplx( 2.0, 1.0 ), cplx( 4.0, -2.0 ) } } );
        check_matrix( Matrix< DDCMat >( cplx( 1.0, -1.0 ) - tComplex ),
                Matrix< DDCMat >{ { cplx( 0.0, -3.0 ), cplx( -2.0, 0.0 ) } } );
    }

    TEST_CASE( "LINALG row and column setters accept expressions and overlapping views", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat > tA = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };
        Matrix< DDRMat > tB = { { 10.0, 20.0, 30.0 }, { 40.0, 50.0, 60.0 }, { 70.0, 80.0, 90.0 } };

        tA.set_row( 1, tB.get_row( 0 ) + tB.get_row( 2 ) );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0, 3.0 }, { 80.0, 100.0, 120.0 }, { 7.0, 8.0, 9.0 } } );
        tA.set_column( 2, trans( tB.get_row( 1 ) ) );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0, 40.0 }, { 80.0, 100.0, 50.0 }, { 7.0, 8.0, 60.0 } } );

        // These assignments cross the destination row or column and require a snapshot.
        tA.set_row( 2, tA.get_column( 2 ) );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0, 40.0 }, { 80.0, 100.0, 50.0 }, { 40.0, 50.0, 60.0 } } );
        tA.set_column( 0, trans( tA.get_row( 2 ) ) );
        check_matrix( tA, Matrix< DDRMat >{ { 40.0, 2.0, 40.0 }, { 50.0, 100.0, 50.0 }, { 60.0, 50.0, 60.0 } } );

        Matrix< F33RMat > tFixed( 3, 3, 0.0 );
        tFixed.set_row( 0, tB.get_column( 1 ) );
        tFixed.set_column( 2, trans( tB.get_row( 0 ) ) );
        check_matrix( tFixed, Matrix< DDRMat >{ { 20.0, 50.0, 10.0 }, { 0.0, 0.0, 20.0 }, { 0.0, 0.0, 30.0 } } );
    }

    TEST_CASE( "LINALG iterators and column pointers expose column-major storage", "[linalgebra],[matrix_api]" )
    {
        using Real_Backend    = deduce_matrix< real, MDynamic, MDynamic >::Matrix_Type;
        using Complex_Backend = deduce_matrix< cplx, MDynamic, MDynamic >::Matrix_Type;
        STATIC_REQUIRE( ( std::is_same< Real_Backend, DDRMat >::value ) );
        STATIC_REQUIRE( ( std::is_same< Complex_Backend, DDCMat >::value ) );

        Matrix< DDRMat >        tA     = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
        const Matrix< DDRMat >& tConst = tA;
        CHECK( tA.size( 0 ) == 2 );
        CHECK( tA.size( 1 ) == 3 );
        CHECK( tA.capacity() == 6 * sizeof( real ) );
        CHECK( std::distance( tA.begin(), tA.end() ) == 6 );
        CHECK( std::distance( tConst.begin(), tConst.end() ) == 6 );
        CHECK( std::distance( tConst.cbegin(), tConst.cend() ) == 6 );
        CHECK( std::accumulate( tConst.begin(), tConst.end(), 0.0 ) == 21.0 );
        CHECK( tConst.data() == tConst.colptr( 0 ) );
        CHECK( tConst.colptr( 1 ) == tConst.data() + 2 );
        CHECK( tConst.colptr( 2 )[ 1 ] == 6.0 );
        std::fill( tA.begin(), tA.end(), 7.0 );
        tA.colptr( 1 )[ 0 ] = 8.0;
        check_matrix( tA, Matrix< DDRMat >{ { 7.0, 8.0, 7.0 }, { 7.0, 7.0, 7.0 } } );

        Matrix< DDRMat > tVector( 3 );
        tVector.fill( 2.0 );
        CHECK( tVector.n_rows() == 3 );
        CHECK( tVector.n_cols() == 1 );
        CHECK( tVector.length() == 3 );

        const Matrix< DDRMat > tEmpty;
        CHECK( tEmpty.begin() == tEmpty.end() );
        CHECK( tEmpty.cbegin() == tEmpty.cend() );
        CHECK( tEmpty.capacity() == 0 );
        CHECK( tEmpty.length() == 0 );
    }

    TEST_CASE( "LINALG resize preserves the rectangular overlap and zeros new entries", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat > tA = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
        tA.resize( 3, 2 );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0 }, { 4.0, 5.0 }, { 0.0, 0.0 } } );
        tA.resize( 1, 1 );
        CHECK( tA( 0 ) == 1.0 );
        tA.set_size( 2, 3, -2.0 );
        check_matrix( tA, Matrix< DDRMat >( 2, 3, -2.0 ) );
        tA.set_size( 0, 0 );
        CHECK( tA.numel() == 0 );
        tA.set_size( 2, 2 );
        tA.fill( 3.0 );
        check_matrix( tA, Matrix< DDRMat >( 2, 2, 3.0 ) );
    }

    TEST_CASE( "LINALG member reshape retains linear order and pads or truncates", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat > tA = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
        tA.reshape( 3, 2 );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 5.0 }, { 4.0, 3.0 }, { 2.0, 6.0 } } );
        tA.reshape( 2, 4 );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0, 3.0, 0.0 }, { 4.0, 5.0, 6.0, 0.0 } } );
        tA.reshape( 1, 3 );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 4.0, 2.0 } } );
        CHECK( tA.length() == 3 );
    }

    TEST_CASE( "LINALG member transpose matches the native Hermitian transpose", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDCMat > tA = { { cplx( 1, 2 ), cplx( 3, -4 ) } };
        tA.inplace_trans();
        check_matrix( tA, Matrix< DDCMat >{ { cplx( 1, -2 ) }, { cplx( 3, 4 ) } } );
        CHECK( tA.length() == 2 );
        tA.inplace_trans();
        check_matrix( tA, Matrix< DDCMat >{ { cplx( 1, 2 ), cplx( 3, -4 ) } } );
    }

    TEST_CASE( "LINALG const spans and row or column output overloads", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat >        tA     = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 }, { 7.0, 8.0, 9.0 } };
        const Matrix< DDRMat >& tConst = tA;
        check_matrix( Matrix< DDRMat >( tConst( { 1, 2 }, { 0, 1 } ) ), Matrix< DDRMat >{ { 4.0, 5.0 }, { 7.0, 8.0 } } );
        check_matrix( Matrix< DDRMat >( tConst( { 0, 1 }, 2 ) ), Matrix< DDRMat >{ { 3.0 }, { 6.0 } } );
        check_matrix( Matrix< DDRMat >( tConst( { 1, 2 } ) ), Matrix< DDRMat >{ { 4.0 }, { 7.0 } } );
        check_matrix( Matrix< DDRMat >( tConst.get_row( 1 ) ), Matrix< DDRMat >{ { 4.0, 5.0, 6.0 } } );
        check_matrix( Matrix< DDRMat >( tConst.get_column( 2 ) ), Matrix< DDRMat >{ { 3.0 }, { 6.0 }, { 9.0 } } );

        // The Armadillo output overloads require a correctly sized destination.
        Matrix< DDRMat > tRow( 1, 3 ), tColumn( 3, 1 );
        tConst.get_row( 1, tRow );
        tConst.get_column( 2, tColumn );
        check_matrix( tRow, Matrix< DDRMat >{ { 4.0, 5.0, 6.0 } } );
        check_matrix( tColumn, Matrix< DDRMat >{ { 3.0 }, { 6.0 }, { 9.0 } } );
        auto tRowView = tA.get_row( 0 );
        tRowView( 1 ) = 12.0;
        CHECK( tA( 0, 1 ) == 12.0 );
        auto tBlock    = tA( { 1, 2 }, { 1, 2 } );
        tBlock( 0, 0 ) = -5.0;
        CHECK( tA( 1, 1 ) == -5.0 );
    }

    TEST_CASE( "LINALG copying pointers and scalar assignments own their values", "[linalgebra],[matrix_api]" )
    {
        const real       tReadOnly[] = { 1.0, 4.0, 2.0, 5.0 };
        Matrix< DDRMat > tA( tReadOnly, 2, 2 );
        CHECK( tA.data() != tReadOnly );
        check_matrix( tA, Matrix< DDRMat >{ { 1.0, 2.0 }, { 4.0, 5.0 } } );
        real             tMutable[] = { 1.0, 4.0, 2.0, 5.0 };
        Matrix< DDRMat > tCopied( tMutable, 2, 2, true, false );
        tMutable[ 0 ] = -9.0;
        CHECK( tCopied( 0, 0 ) == 1.0 );
        tA += 2.0;
        tA -= 1.0;
        check_matrix( tA, Matrix< DDRMat >{ { 2.0, 3.0 }, { 5.0, 6.0 } } );
        tA = 8.0;
        REQUIRE( tA.n_rows() == 1 );
        REQUIRE( tA.n_cols() == 1 );
        CHECK( tA( 0 ) == 8.0 );
        const Matrix< DDRMat > tCopy = tCopied.copy();
        CHECK( tCopy.data() != tCopied.data() );
        tCopied( 0, 0 ) = 99.0;
        CHECK( tCopy( 0, 0 ) == 1.0 );
        CHECK( tCopy.matrix_data()( 1, 1 ) == 5.0 );
    }

    TEST_CASE( "LINALG extrema return coordinates as well as values", "[linalgebra],[matrix_api]" )
    {
        const Matrix< DDRMat > tA   = { { 4.0, -8.0, 2.0 }, { 11.0, 3.0, 7.0 } };
        uint                   tRow = 0, tColumn = 0;
        CHECK( tA.min( tRow, tColumn ) == -8.0 );
        CHECK( tRow == 0 );
        CHECK( tColumn == 1 );
        CHECK( tA.max( tRow, tColumn ) == 11.0 );
        CHECK( tRow == 1 );
        CHECK( tColumn == 0 );
        const Matrix< DDCMat > tComplex = { { cplx( 0, 3 ), cplx( 1 ) }, { cplx( 0, -5 ), cplx( 2 ) } };
        CHECK( tComplex.min() == cplx( 1 ) );
        CHECK( tComplex.max() == cplx( 0, -5 ) );
        CHECK( tComplex.min( tRow, tColumn ) == cplx( 1 ) );
        CHECK( tRow == 0 );
        CHECK( tColumn == 1 );
        CHECK( tComplex.max( tRow, tColumn ) == cplx( 0, -5 ) );
        CHECK( tRow == 1 );
        CHECK( tColumn == 0 );
    }

    TEST_CASE( "LINALG move helper and stream overloads", "[linalgebra],[matrix_api]" )
    {
        Matrix< DDRMat >       tA     = { { 1.0, 2.0, 3.0 } };
        const Matrix< DDRMat > tMoved = moris::move( tA );
        CHECK( tA.n_rows() == 0 );
        CHECK( tA.n_cols() == 0 );
        check_matrix( tMoved, Matrix< DDRMat >{ { 1.0, 2.0, 3.0 } } );
        std::ostringstream tValues, tPointer, tNull;
        tValues << tMoved;
        tPointer << &tMoved;
        tNull << static_cast< const Matrix< DDRMat >* >( nullptr );
        CHECK( tValues.str() == "1,2,3" );
        CHECK( tPointer.str() == tValues.str() );
        CHECK( tNull.str() == "NULL" );
    }
    TEST_CASE( "LINALG const matrix comparisons retain elementwise truth values", "[linalgebra],[matrix_api]" )
    {
        const Matrix< DDRMat > tA = { { -1.0, 2.0 } };
        const Matrix< DDRMat > tB = { { 0.0, 2.0 } };
        check_matrix( Matrix< DDBMat >( tA == tB ), Matrix< DDBMat >{ { false, true } } );
        check_matrix( Matrix< DDBMat >( tA != tB ), Matrix< DDBMat >{ { true, false } } );
        check_matrix( Matrix< DDBMat >( tA < tB ), Matrix< DDBMat >{ { true, false } } );
        check_matrix( Matrix< DDBMat >( tA <= tB ), Matrix< DDBMat >{ { true, true } } );
        check_matrix( Matrix< DDBMat >( tA > tB ), Matrix< DDBMat >{ { false, false } } );
        check_matrix( Matrix< DDBMat >( tA >= tB ), Matrix< DDBMat >{ { false, true } } );
    }

}    // namespace moris
