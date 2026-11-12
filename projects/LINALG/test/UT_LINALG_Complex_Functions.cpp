/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt in the MORIS root.
 */
#include "UT_LINALG_Test_Helpers.hpp"
#include "fn_chol_u.hpp"
#include "fn_comp_abs.hpp"
#include "fn_ctrans.hpp"
#include "fn_det.hpp"
#include "fn_inv.hpp"
#include "fn_isfinite.hpp"
#include "fn_norm.hpp"
#include "fn_svd.hpp"
#include "fn_sylvester.hpp"
#include <limits>

namespace moris
{
    using linalg_test::check_matrix;

    TEST_CASE( "LINALG complex upper Cholesky and inverse preserve Hermitian algebra", "[linalgebra],[complex_functions]" )
    {
        const Matrix< DDCMat > tA     = { { cplx( 3 ), cplx( 0, 1 ) }, { cplx( 0, -1 ), cplx( 3 ) } };
        const Matrix< DDCMat > tUpper = chol_u( tA );
        CHECK( norm( ctrans( tUpper ) * tUpper - tA ) < 1e-12 );
        check_matrix( Matrix< DDCMat >( inv( tA ) ), Matrix< DDCMat >{ { cplx( 3.0 / 8 ), cplx( 0, -1.0 / 8 ) }, { cplx( 0, 1.0 / 8 ), cplx( 3.0 / 8 ) } } );
        CHECK( std::abs( det( tA ) - cplx( 8 ) ) < 1e-12 );
        const Matrix< DDCMat > tValues = { { cplx( 3, 4 ), cplx( 0, -2 ) } };
        check_matrix( Matrix< DDRMat >( comp_abs( tValues ) ), Matrix< DDRMat >{ { 5.0, 2.0 } } );
        CHECK( isfinite( tA ) );
        Matrix< DDCMat > tInvalid = tA;
        tInvalid( 0, 0 )          = cplx( 0, std::numeric_limits< real >::infinity() );
        CHECK_FALSE( isfinite( tInvalid ) );
    }

    TEST_CASE( "LINALG complex rectangular SVD uses adjoints and real singular values", "[linalgebra],[complex_functions]" )
    {
        const Matrix< DDCMat > tA = { { cplx( 0, 2 ), cplx( 0 ) }, { cplx( 0 ), cplx( 3 ) }, { cplx( 0 ), cplx( 0 ) } };
        Matrix< DDCMat >       tU, tV;
        Matrix< DDRMat >       tS, tValuesOnly;
        svd( tU, tS, tV, tA );
        check_matrix( tS, Matrix< DDRMat >{ { 3.0 }, { 2.0 } } );
        svd( tValuesOnly, tA );
        check_matrix( tValuesOnly, tS );
        Matrix< DDCMat > tSigma( 3, 2, cplx( 0 ) );
        tSigma( 0, 0 ) = tS( 0 );
        tSigma( 1, 1 ) = tS( 1 );
        CHECK( norm( tU * tSigma * ctrans( tV ) - tA ) < 1e-12 );
    }

    TEST_CASE( "LINALG complex rectangular Sylvester equation has the known solution", "[linalgebra],[complex_functions]" )
    {
        const Matrix< DDCMat > tA = { { cplx( 1, 1 ), cplx( 2 ) }, { cplx( 0 ), cplx( 3, -1 ) } };
        const Matrix< DDCMat > tB( 1, 1, cplx( 2 ) );
        const Matrix< DDCMat > tC = { { cplx( -6, -2 ) }, { cplx( -9, 7 ) } };
        const Matrix< DDCMat > tX = sylvester( tA, tB, tC );
        check_matrix( tX, Matrix< DDCMat >{ { cplx( 1, 1 ) }, { cplx( 2, -1 ) } } );
        CHECK( norm( tA * tX + tX * tB + tC ) < 1e-12 );
    }
}    // namespace moris
