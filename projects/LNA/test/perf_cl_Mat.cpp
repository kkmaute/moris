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

#ifdef PERF_MAT

// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::MatPerformanceStudy",
        "[linalgebra],[MatPerformanceStudy]" )
{
    moris::uint nmats = 2;
    moris::uint nruns = 10;
    moris::uint matsize = 300;
    moris::uint nrows = 0;
    moris::uint ncols = 0;

    moris::Mat< moris::real > TimeMat_moris( nmats, 4 );
    moris::Mat< moris::real > TimeMat_pure( nmats, 4 );

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        nrows = matsize + n * matsize;
        ncols = matsize + n * matsize;
        moris::Mat< moris::real > A ( nrows, ncols, 1.0 );
        moris::Mat< moris::real > B ( nrows, ncols, 1.0 );
        moris::Mat< moris::real > C ( nrows, ncols, 1.0 );
        moris::Mat< moris::real > D ( nrows, ncols, 1.0 );
        moris::Mat< moris::real > E ( nrows, ncols );
        moris::Mat< moris::real > F ( nrows, ncols );
        moris::Mat< moris::real > G ( nrows, ncols );
        moris::Mat< moris::real > H ( nrows, ncols );
        moris::Mat< moris::real > H1( nrows, ncols );
        moris::Mat< moris::real > H2( nrows, ncols );
        moris::Mat< moris::real > At ( nrows, ncols );
        moris::Mat< moris::real > Ct ( nrows, ncols );

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            E = moris::trans( A ) + B + moris::trans( C ) + D;
        }
        tE.elapsed();

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            F = moris::trans( A ) - B - moris::trans( C ) - D;
        }
        tF.elapsed();

        moris::real wall_time_microseconds_tF = tF.toc<moris::chronos::microseconds>().wall;

        // operator multiplication
        moris::tic tG;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            G = moris::trans( A ) * B * moris::trans( C ) * D;
        }
        tG.elapsed();

        moris::real wall_time_microseconds_tG = tG.toc<moris::chronos::microseconds>().wall;

        // operator elementwise division
        moris::tic tH;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = moris::trans( A );
            Ct = moris::trans( C );
            H1 = At / B;
            H2 = H1 / Ct;
            H  = H2 / D;
        }
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
        arma::Mat< moris::real > A ( nrows, ncols );
        arma::Mat< moris::real > B ( nrows, ncols );
        arma::Mat< moris::real > C ( nrows, ncols );
        arma::Mat< moris::real > D ( nrows, ncols );
        arma::Mat< moris::real > E ( nrows, ncols );
        arma::Mat< moris::real > F ( nrows, ncols );
        arma::Mat< moris::real > G ( nrows, ncols );
        arma::Mat< moris::real > H ( nrows, ncols );
        arma::Mat< moris::real > H1( nrows, ncols );
        arma::Mat< moris::real > H2( nrows, ncols );
        arma::Mat< moris::real > At ( nrows, ncols );
        arma::Mat< moris::real > Ct ( nrows, ncols );

        A.fill( 1.0 );
        B.fill( 1.0 );
        C.fill( 1.0 );
        D.fill( 1.0 );

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
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = arma::strans( A );
            Ct = arma::strans( C );
            H1 = At / B;
            H2 = H1 / Ct;
            H  = H2 / D;
        }
        tH.elapsed();

        moris::real wall_time_microseconds_tH = tH.toc<moris::chronos::microseconds>().wall;

        TimeMat_pure( n, 0 ) = wall_time_microseconds_tE;
        TimeMat_pure( n, 1 ) = wall_time_microseconds_tF;
        TimeMat_pure( n, 2 ) = wall_time_microseconds_tG;
        TimeMat_pure( n, 3 ) = wall_time_microseconds_tH;

    }

#elif MORIS_USE_EIGEN

    typedef Eigen::Matrix< moris::real, Eigen::Dynamic, Eigen::Dynamic > tMatrix;

    for ( moris::uint n = 0; n < nmats; n++ )
    {
        nrows = matsize + n * matsize;
        ncols = matsize + n * matsize;
        tMatrix A( nrows, ncols );
        tMatrix B( nrows, ncols );
        tMatrix C( nrows, ncols );
        tMatrix D( nrows, ncols );
        tMatrix E( nrows, ncols );
        tMatrix F( nrows, ncols );
        tMatrix G( nrows, ncols );
        tMatrix H( nrows, ncols );
        tMatrix H1( nrows, ncols );
        tMatrix H2( nrows, ncols );
        tMatrix At( nrows, ncols );
        tMatrix Ct( nrows, ncols );

        A.fill( 1.0 );
        B.fill( 1.0 );
        C.fill( 1.0 );
        D.fill( 1.0 );

        // operator plus
        moris::tic tE;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            E = A.transpose() + B + C.transpose() + D;
        }
        tE.elapsed();

        moris::real wall_time_microseconds_tE = tE.toc<moris::chronos::microseconds>().wall;

        // operator minus
        moris::tic tF;
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            F = A.transpose() - B - C.transpose() - D;
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
        for ( moris::uint k = 0; k < nruns; k++ )
        {
            At = A.transpose();
            Ct = C.transpose();
            H1 = At.cwiseQuotient( B );
            H2 = H1.cwiseQuotient( Ct );
            H  = H2.cwiseQuotient( D );
        }
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
