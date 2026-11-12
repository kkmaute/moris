/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */

#include <catch.hpp>

#include "cl_Matrix.hpp"
#include "fn_dot.hpp"
#include "fn_find.hpp"
#include "fn_diag_mat.hpp"
#include "fn_reshape.hpp"
#include "fn_vectorize.hpp"
#include "op_elemwise_mult.hpp"
#include "op_elemwise_div.hpp"
#include "op_greater_equal.hpp"
#include "op_equal_equal.hpp"
#include "op_not_equal.hpp"
#include "fn_sum.hpp"
#include "fn_prod.hpp"
#include "fn_trans.hpp"
#include "linalg_typedefs.hpp"

#include <type_traits>
#include <utility>

namespace moris
{
    namespace
    {
        void check_expression_result( const Matrix< DDRMat >& aActual,
                const Matrix< DDRMat >&                       aExpected )
        {
            REQUIRE( aActual.n_rows() == aExpected.n_rows() );
            REQUIRE( aActual.n_cols() == aExpected.n_cols() );

            for ( uint j = 0; j < aActual.n_cols(); ++j )
            {
                for ( uint i = 0; i < aActual.n_rows(); ++i )
                {
                    CHECK( aActual( i, j ) == Approx( aExpected( i, j ) ).margin( 1e-12 ) );
                }
            }
        }
    }    // namespace

    TEST_CASE( "LINALG scalar shifts compose with lazy backend expressions", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        // Keep the nested operand alive, as required by native Armadillo expressions.
        auto tScaled = 2.0 * tA;
        auto tShift  = tScaled + 3.0;
        tA( 0, 0 )   = 5.0;
        check_expression_result( tShift, { { 13.0, 7.0 }, { 9.0, 11.0 } } );
        check_expression_result( 3.0 + ( 2.0 * tA ), { { 13.0, 7.0 }, { 9.0, 11.0 } } );
        check_expression_result( ( 2.0 * tA ) - 3.0, { { 7.0, 1.0 }, { 3.0, 5.0 } } );
        check_expression_result( 3.0 - ( 2.0 * tA ), { { -7.0, -1.0 }, { -3.0, -5.0 } } );
        check_expression_result( tA + 1.0, { { 6.0, 3.0 }, { 4.0, 5.0 } } );
        check_expression_result( 1.0 + tA, { { 6.0, 3.0 }, { 4.0, 5.0 } } );
        check_expression_result( tA - 1.0, { { 4.0, 1.0 }, { 2.0, 3.0 } } );
        check_expression_result( 1.0 - tA, { { -4.0, -1.0 }, { -2.0, -3.0 } } );
    }

    TEST_CASE( "LINALG expressions retain backend types and const operands", "[linalgebra],[expressions]" )
    {
        using Wrapper = Matrix< DDRMat >;
        using Native  = std::decay_t< decltype( std::declval< Wrapper& >().matrix_data() ) >;
        using Sum     = decltype( std::declval< const Wrapper& >() + std::declval< const Wrapper& >() );
        using Product = decltype( std::declval< const Wrapper& >() * std::declval< const Wrapper& >() );

        // Operators must expose the backend's expression tree, not introduce
        // owning intermediates. Compare to native types, including constness.
        STATIC_REQUIRE( ( std::is_same< Sum,
                decltype( std::declval< const Native& >() + std::declval< const Native& >() ) >::value ) );
        STATIC_REQUIRE( ( std::is_same< Product,
                decltype( std::declval< const Native& >() * std::declval< const Native& >() ) >::value ) );
        STATIC_REQUIRE_FALSE( ( std::is_same< std::decay_t< Sum >, Native >::value ) );

        const Wrapper tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Wrapper tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
        const Wrapper tD = { { 2.0, -1.0 }, { 0.0, 3.0 } };

        check_expression_result( ( tA + tB ) + tD, { { 8.0, 7.0 }, { 10.0, 15.0 } } );
        check_expression_result( tD + ( tA + tB ), { { 8.0, 7.0 }, { 10.0, 15.0 } } );
        check_expression_result( ( tA + tB ) - tD, { { 4.0, 9.0 }, { 10.0, 9.0 } } );
        check_expression_result( tD - ( tA + tB ), { { -4.0, -9.0 }, { -10.0, -9.0 } } );
        check_expression_result( trans( tA ) * tD, { { 2.0, 8.0 }, { 4.0, 10.0 } } );
        check_expression_result( ( tA + tB ) * tD, { { 12.0, 18.0 }, { 20.0, 26.0 } } );
        check_expression_result( tD * ( tA + tB ), { { 2.0, 4.0 }, { 30.0, 36.0 } } );
    }

    TEST_CASE( "LINALG expression assignment reuses destination storage", "[linalgebra],[expressions]" )
    {
        const Matrix< DDRMat > tA( 64, 64, 2.0 );
        const Matrix< DDRMat > tB( 64, 64, 3.0 );
        Matrix< DDRMat >       tC( 64, 64, 0.0 );
        auto*                  tStorage = tC.data();

        // Moving an implicitly materialized wrapper into C changes its storage.
        // Same-size coefficient-wise evaluation should reuse the allocation.
        tC = 2.0 * tA + tB;
        CHECK( tC.data() == tStorage );
        CHECK( tC( 0, 0 ) == 7.0 );
        CHECK( tC( 63, 63 ) == 7.0 );

        tC += tA + tB;
        CHECK( tC.data() == tStorage );
        CHECK( tC( 63, 63 ) == 12.0 );

        tC -= 2.0 * tA + tB;
        CHECK( tC.data() == tStorage );
        CHECK( tC( 63, 63 ) == 5.0 );

        tC += tA;
        tC -= tB;
        CHECK( tC.data() == tStorage );
        CHECK( tC( 63, 63 ) == 4.0 );

        tC = moris::trans( tA + tB );
        CHECK( tC.data() == tStorage );
        CHECK( tC( 63, 63 ) == 5.0 );
    }

    TEST_CASE( "LINALG named expressions defer evaluation", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat >       tA( 2, 2, 1.0 );
        const Matrix< DDRMat > tB( 2, 2, 2.0 );
        auto                   tSum       = tA + tB;
        auto                   tTranspose = moris::trans( tSum );

        // Keep both operands and the intermediate expression alive. Changes
        // before evaluation must remain visible through the expression tree.
        tA( 0, 1 ) = 5.0;
        check_expression_result( tTranspose, { { 3.0, 3.0 }, { 7.0, 3.0 } } );
    }

    TEST_CASE( "LINALG assignment retains native product alias handling", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat >       tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat > tB = { { 2.0, 0.0 }, { 1.0, 3.0 } };

        // Do not impose noalias(): each backend must handle a destination that
        // also participates in the product, including compound assignments.
        tA = tA * tA;
        check_expression_result( tA, { { 7.0, 10.0 }, { 15.0, 22.0 } } );
        tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        tA += tA * tB;
        check_expression_result( tA, { { 5.0, 8.0 }, { 13.0, 16.0 } } );
        tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        tA -= tA * tB;
        check_expression_result( tA, { { -3.0, -4.0 }, { -7.0, -8.0 } } );
    }

    TEST_CASE( "LINALG self transpose assignment remains supported", "[linalgebra],[expressions]" )
    {
        // These assignments previously evaluated through an owning temporary.
        // Direct expression assignment must preserve their values and dimensions.
        Matrix< DDRMat > tSquare = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        tSquare                  = trans( tSquare );
        check_expression_result( tSquare, { { 1.0, 3.0 }, { 2.0, 4.0 } } );

        Matrix< DDRMat > tRectangle = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
        tRectangle                  = trans( tRectangle );
        check_expression_result( tRectangle, { { 1.0, 4.0 }, { 2.0, 5.0 }, { 3.0, 6.0 } } );

        Matrix< DDRMat > tVector = { { 1.0, 2.0, 3.0 } };
        tVector                  = trans( tVector );
        check_expression_result( tVector, { { 1.0 }, { 2.0 }, { 3.0 } } );

        Matrix< DDCMat > tComplex( 1, 2 );
        tComplex( 0 ) = cplx( 1.0, 2.0 );
        tComplex( 1 ) = cplx( 3.0, -4.0 );
        tComplex      = trans( tComplex );
        REQUIRE( tComplex.n_rows() == 2 );
        REQUIRE( tComplex.n_cols() == 1 );
        CHECK( tComplex( 0 ) == cplx( 1.0, 2.0 ) );
        CHECK( tComplex( 1 ) == cplx( 3.0, -4.0 ) );
    }

    TEST_CASE( "LINALG expression transpose and reductions compose", "[linalgebra],[expressions]" )
    {
        const Matrix< DDRMat > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
        const Matrix< DDRMat > tD = { { 2.0, -1.0 }, { 0.0, 3.0 } };

        // Retain the native lazy transpose rather than returning a dense matrix.
        using Transpose = decltype( moris::trans( tA + tB ) );
#ifdef MORIS_USE_ARMA
        STATIC_REQUIRE( ( std::is_same< Transpose, decltype( arma::strans( tA.matrix_data() + tB.matrix_data() ) ) >::value ) );
#else
        STATIC_REQUIRE( ( std::is_same< Transpose, decltype( ( tA.matrix_data() + tB.matrix_data() ).transpose() ) >::value ) );
#endif
        check_expression_result( moris::trans( tA + tB ) * tD, { { 12.0, 24.0 }, { 16.0, 28.0 } } );
        check_expression_result( trans( tA + tB ), { { 6.0, 10.0 }, { 8.0, 12.0 } } );
        // Integration weights compose a product, vectorization and unqualified transpose.
        check_expression_result( trans( vectorize( trans( tA ) * tD ) ), { { 2.0, 4.0, 8.0, 10.0 } } );
        CHECK( moris::sum( tA ) == Approx( 10.0 ) );
        CHECK( moris::sum( tA + tB ) == Approx( 36.0 ) );
        CHECK( moris::dot( tA, tB ) == Approx( 70.0 ) );
        CHECK( moris::dot( tA, tB + tB ) == Approx( 140.0 ) );
        CHECK( moris::dot( tB + tB, tA ) == Approx( 140.0 ) );
        CHECK( moris::dot( tA + tB, tA - tB ) == Approx( -144.0 ) );

        // Dimensioned reductions remain lazy and distinct from the scalar sum.
        check_expression_result( moris::sum( tA + tB, 0 ), { { 16.0, 20.0 } } );
        check_expression_result( moris::sum( tA + tB, 1 ), { { 14.0 }, { 22.0 } } );

        // A transpose must not conjugate a complex expression.
        Matrix< DDCMat > tComplexA( 2, 1 );
        Matrix< DDCMat > tComplexB( 2, 1 );
        tComplexA( 0 )               = cplx( 1.0, 2.0 );
        tComplexA( 1 )               = cplx( 3.0, -4.0 );
        tComplexB( 0 )               = cplx( 2.0, 1.0 );
        tComplexB( 1 )               = cplx( 4.0, 3.0 );
        Matrix< DDCMat > tTransposed = moris::trans( tComplexA + tComplexB );
        REQUIRE( tTransposed.n_rows() == 1 );
        REQUIRE( tTransposed.n_cols() == 2 );
        CHECK( tTransposed( 0, 0 ) == cplx( 3.0, 3.0 ) );
        CHECK( tTransposed( 0, 1 ) == cplx( 7.0, -1.0 ) );
    }
    TEST_CASE( "LINALG nested transpose aliases retain source coefficients", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat >       tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
        tA                        = trans( tA + tB );
        check_expression_result( tA, { { 6.0, 10.0 }, { 8.0, 12.0 } } );

        tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        tA = tA + trans( tA );
        check_expression_result( tA, { { 2.0, 5.0 }, { 5.0, 8.0 } } );
        tA += trans( tA + tB );
        check_expression_result( tA, { { 9.0, 17.0 }, { 16.0, 24.0 } } );

        Matrix< DDRMat > tRectangle = { { 1.0, 2.0, 3.0 }, { 4.0, 5.0, 6.0 } };
        tRectangle                  = trans( 2.0 * tRectangle );
        check_expression_result( tRectangle, { { 2.0, 8.0 }, { 4.0, 10.0 }, { 6.0, 12.0 } } );
    }

    TEST_CASE( "LINALG reordered expressions retain native nesting rules", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat >       tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        const Matrix< DDRMat > tB = { { 5.0, 6.0 }, { 7.0, 8.0 } };
        const Matrix< DDRMat > tD = { { 2.0, -1.0 }, { 0.0, 3.0 } };
        check_expression_result( moris::sum( ( tA + tB ) * tD, 0 ), { { 32.0, 44.0 } } );
        check_expression_result( moris::vectorize( ( tA + tB ) * tD, true ), { { 12.0, 18.0, 20.0, 26.0 } } );
        check_expression_result( moris::reshape( ( tA + tB ) * tD, 4, 1 ), { { 12.0 }, { 20.0 }, { 18.0 }, { 26.0 } } );

        auto tReshaped = moris::reshape( tA, 4, 1 );
        tA( 0, 0 )     = 9.0;
        check_expression_result( tReshaped, { { 9.0 }, { 3.0 }, { 2.0 }, { 4.0 } } );
        tA = moris::diag_mat( tA );
        check_expression_result( tA, { { 9.0, 0.0 }, { 0.0, 4.0 } } );
        tA = moris::reshape( tA, 4, 1 );
        check_expression_result( tA, { { 9.0 }, { 0.0 }, { 0.0 }, { 4.0 } } );
    }

}    // namespace moris

// Native expressions have no moris operand for argument-dependent lookup.
TEST_CASE( "LINALG native expression operators resolve outside moris", "[linalgebra],[expressions]" )
{
    using M          = moris::Matrix< moris::DDRMat >;
    const M tA       = { { 1.0, 2.0 }, { 3.0, 4.0 } };
    const M tB       = { { 5.0, 6.0 }, { 7.0, 8.0 } };
    const M tProduct = ( tA + tB ) % ( tA - tB );
    CHECK( tProduct( 0, 0 ) == -24.0 );
    CHECK( tProduct( 1, 1 ) == -48.0 );
    const M tQuotient = ( tA + tB ) / ( tA - tB );
    CHECK( tQuotient( 0, 0 ) == -1.5 );
    CHECK( tQuotient( 1, 1 ) == -3.0 );
    const moris::Matrix< moris::DDBMat > tComparison = ( tA + tB ) >= ( tA - tB );
    CHECK( tComparison( 0, 0 ) );
    CHECK( tComparison( 1, 1 ) );
    const moris::Matrix< moris::DDNIMat > tFound = moris::find( ( tA + tB ) >= ( tA - tB ) );
    CHECK( tFound.numel() == 4 );
    const moris::Matrix< moris::DDBMat > tEqual = ( tA + tB ) == ( tB + tA );
    REQUIRE( tEqual.n_rows() == 2 );
    REQUIRE( tEqual.n_cols() == 2 );
    CHECK( tEqual( 0, 0 ) );
    CHECK( tEqual( 1, 1 ) );
    const moris::Matrix< moris::DDBMat > tUnequal = ( tA + tB ) != ( tA - tB );
    REQUIRE( tUnequal.n_rows() == 2 );
    REQUIRE( tUnequal.n_cols() == 2 );
    CHECK( tUnequal( 0, 0 ) );
    CHECK( tUnequal( 1, 1 ) );
}

namespace moris
{
    TEST_CASE( "LINALG reductions retain deferred reads and destination aliases", "[linalgebra],[expressions]" )
    {
        Matrix< DDRMat > tA       = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        auto             tScaled  = 2.0 * tA;
        auto             tSum     = moris::sum( tScaled, 1 );
        auto             tProduct = moris::prod( tScaled );
        tA( 0, 0 )                = 5.0;
        check_expression_result( tSum, { { 14.0 }, { 14.0 } } );
        check_expression_result( tProduct, { { 60.0, 32.0 } } );

        tA = moris::sum( tA, 0 );
        check_expression_result( tA, { { 8.0, 6.0 } } );
        tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
        tA = moris::prod( tA );
        check_expression_result( tA, { { 3.0, 8.0 } } );
    }

#ifdef MORIS_USE_EIGEN
    TEST_CASE( "LINALG Eigen reduction copies share evaluated operands and retain lifetime", "[linalgebra],[expressions]" )
    {
        // A product snapshot must survive its operands and intermediate
        // expressions without being cloned when Eigen copies its functors.
        auto tSum = []() {
            const Matrix< DDRMat > tA         = { { 1.0, 2.0 }, { 3.0, 4.0 } };
            const Matrix< DDRMat > tB         = { { 2.0, 0.0 }, { 1.0, 3.0 } };
            auto                   tReduction = moris::sum( tA * tB, 0 );
            auto                   tCopy      = tReduction;
            CHECK( tCopy.functor().values.data() == tReduction.functor().values.data() );
            return tCopy;
        }();
        auto tCopy = tSum;
        CHECK( tCopy.functor().values.data() == tSum.functor().values.data() );
        auto tMoved = std::move( tCopy );
        check_expression_result( tMoved, { { 14.0, 18.0 } } );
        const Matrix< DDRMat > tWeights = { { 2.0, 3.0 } };
        check_expression_result( ( 0.5 * tSum ) % tWeights, { { 14.0, 27.0 } } );

        auto tProduct = []() {
            const Matrix< DDRMat > tA = { { 1.0, 2.0 }, { 3.0, 4.0 } };
            const Matrix< DDRMat > tB = { { 2.0, 0.0 }, { 1.0, 3.0 } };
            return moris::prod( tA * tB );
        }();
        auto tProductCopy = tProduct;
        CHECK( tProductCopy.functor().values.data() == tProduct.functor().values.data() );
        check_expression_result( tProductCopy, { { 40.0, 72.0 } } );
    }
#endif
}    // namespace moris
