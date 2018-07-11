// Third-party header files.
#include <catch.hpp>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#ifdef MORIS_USE_ARMA
#include <armadillo>
#endif
#ifdef MORIS_USE_EIGEN
#include <Eigen>
#endif

// MORIS project header files.
#include "algorithms.hpp"
#include "linalg.hpp"
#include "chronos.hpp"

#ifdef PERF_SP_MAT

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::Sp_MatPerformanceStudy",
        "[linalgebra],[Sp_MatPerformanceStudy]" )
{
    moris::uint nmats = 2;
    moris::uint nruns = 10;
    moris::uint matsize = 200;
    moris::uint nrows = 0;
    moris::uint ncols = 0;

    moris::Mat< moris::real > TimeMat_moris( nmats, 4 );
    moris::Mat< moris::real > TimeMat_pure( nmats, 4 );

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        nrows = matsize + n * matsize;
        ncols = matsize + n * matsize;
        moris::Sp_Mat< moris::real > A ( nrows, ncols );
        moris::Sp_Mat< moris::real > B ( nrows, ncols );
        moris::Sp_Mat< moris::real > C ( nrows, ncols );
        moris::Sp_Mat< moris::real > D ( nrows, ncols );
        moris::Sp_Mat< moris::real > E ( nrows, ncols );
        moris::Sp_Mat< moris::real > F ( nrows, ncols );
        moris::Sp_Mat< moris::real > G ( nrows, ncols );
        moris::Sp_Mat< moris::real > H ( nrows, ncols );
        moris::Sp_Mat< moris::real > H1( nrows, ncols );
        moris::Sp_Mat< moris::real > H2( nrows, ncols );
        moris::Sp_Mat< moris::real > At ( nrows, ncols );
        moris::Sp_Mat< moris::real > Ct ( nrows, ncols );

        for ( moris::uint i = 0; i < nrows; i++ )
        {
            for ( moris::uint j = 0; j < ncols; j++ )
            {
                A( i, j ) = i * j * 1.0;
                B( i, j ) = i * j * 1.0;
                C( i, j ) = i * j * 1.0;
                D( i, j ) = i * j * 1.0;
            }
        }

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = moris::trans( A );
            Ct = moris::trans( C );
            E = At + B + Ct + D;
        }
        tE.elapsed();

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = moris::trans( A );
            Ct = moris::trans( C );
            F = At - B - Ct - D;
        }
        tF.elapsed();

        moris::real wall_time_microseconds_tF = tF.toc<moris::chronos::microseconds>().wall;

        // operator multiplication
        moris::tic tG;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = moris::trans( A );
            Ct = moris::trans( C );
            G = At * B * Ct * D;
        }
        tG.elapsed();

        moris::real wall_time_microseconds_tG = tG.toc<moris::chronos::microseconds>().wall;

        // operator elementwise division
        moris::tic tH;
//        for ( moris::uint k = 0; k < nruns; k++ )
//        {
//            H1 = A / B;
//            H2 = H1 / C;
//            H  = H2 / D;
//        }
        tH.elapsed();

        moris::real wall_time_microseconds_tH = tH.toc<moris::chronos::microseconds>().wall;

        TimeMat_moris( n, 0 ) = wall_time_microseconds_tE;
        TimeMat_moris( n, 1 ) = wall_time_microseconds_tF;
        TimeMat_moris( n, 2 ) = wall_time_microseconds_tG;
        TimeMat_moris( n, 3 ) = wall_time_microseconds_tH;

    }

#ifdef MORIS_USE_ARMA

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        nrows = matsize + n * matsize;
        ncols = matsize + n * matsize;
        arma::SpMat< moris::real > A ( nrows, ncols );
        arma::SpMat< moris::real > B ( nrows, ncols );
        arma::SpMat< moris::real > C ( nrows, ncols );
        arma::SpMat< moris::real > D ( nrows, ncols );
        arma::SpMat< moris::real > E ( nrows, ncols );
        arma::SpMat< moris::real > F ( nrows, ncols );
        arma::SpMat< moris::real > G ( nrows, ncols );
        arma::SpMat< moris::real > H ( nrows, ncols );
        arma::SpMat< moris::real > H1( nrows, ncols );
        arma::SpMat< moris::real > H2( nrows, ncols );

        for ( moris::uint i = 0; i < nrows; i++ )
        {
            for ( moris::uint j = 0; j < ncols; j++ )
            {
                A( i, j ) = i * j * 1.0;
                B( i, j ) = i * j * 1.0;
                C( i, j ) = i * j * 1.0;
                D( i, j ) = i * j * 1.0;
            }
        }

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            E = arma::strans( A ) + B + arma::strans( C ) + D;
        }
        tE.elapsed();

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            F = arma::strans( A ) - B - arma::strans( C ) - D;
        }
        tF.elapsed();

        moris::real wall_time_microseconds_tF = tF.toc<moris::chronos::microseconds>().wall;

        // operator multiplication
        moris::tic tG;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            G = arma::strans( A ) * B * arma::strans( C ) * D;
        }
        tG.elapsed();

        moris::real wall_time_microseconds_tG = tG.toc<moris::chronos::microseconds>().wall;

        // operator elementwise division
        moris::tic tH;
//        for ( moris::uint k = 0; k < nruns; k++ )
//        {
//            H1 = A / B;
//            H2 = H1 / C;
//            H  = H2 / D;
//        }
        tH.elapsed();

        moris::real wall_time_microseconds_tH = tH.toc<moris::chronos::microseconds>().wall;

        TimeMat_pure( n, 0 ) = wall_time_microseconds_tE;
        TimeMat_pure( n, 1 ) = wall_time_microseconds_tF;
        TimeMat_pure( n, 2 ) = wall_time_microseconds_tG;
        TimeMat_pure( n, 3 ) = wall_time_microseconds_tH;

    }

#elif MORIS_USE_EIGEN

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        nrows = matsize + n * matsize;
        ncols = matsize + n * matsize;
        Eigen::SparseMatrix<  moris::real > A( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > B( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > C( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > D( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > E( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > F( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > G( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > H( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > H1( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > H2( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > At( nrows, ncols );
        Eigen::SparseMatrix<  moris::real > Ct( nrows, ncols );

        for ( moris::size_t i = 0; i < nrows; i++ )
        {
            for ( moris::size_t j = 0; j < ncols; j++ )
            {
                A.coeffRef( i, j ) = i * j * 1.0;
                B.coeffRef( i, j ) = i * j * 1.0;
                C.coeffRef( i, j ) = i * j * 1.0;
                D.coeffRef( i, j ) = i * j * 1.0;
            }
        }

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = A.transpose();
            Ct = C.transpose();
            E  = At + B + Ct + D;
        }
        tE.elapsed();

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = A.transpose();
            Ct = C.transpose();
            F = At - B - Ct - D;
        }
        tF.elapsed();

        moris::real wall_time_microseconds_tF = tF.toc<moris::chronos::microseconds>().wall;

        // operator multiplication
        moris::tic tG;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            G = A.transpose() * B * C.transpose() * D;
        }
        tG.elapsed();

        moris::real wall_time_microseconds_tG = tG.toc<moris::chronos::microseconds>().wall;

        // operator elementwise division
        moris::tic tH;
//        for ( moris::uint k = 0; k < nruns; k++ )
//        {
//            H1 = A.cwiseQuotient( B );
//            H2 = H1.cwiseQuotient( C );
//            H  = H2.cwiseQuotient( D );
//        }
        tH.elapsed();

        moris::real wall_time_microseconds_tH = tH.toc<moris::chronos::microseconds>().wall;

        TimeMat_pure( n, 0 ) = wall_time_microseconds_tE;
        TimeMat_pure( n, 1 ) = wall_time_microseconds_tF;
        TimeMat_pure( n, 2 ) = wall_time_microseconds_tG;
        TimeMat_pure( n, 3 ) = wall_time_microseconds_tH;

    }

#endif

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        for ( moris::uint m = 0; m < 4; m++ )
        {
            REQUIRE( (std::abs(TimeMat_pure( n, m ) - TimeMat_moris( n, m ) ) ) / TimeMat_pure( n, m ) <= 0.10 );
        }
    }
}
#endif
