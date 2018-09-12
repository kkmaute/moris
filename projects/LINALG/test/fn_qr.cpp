// Third-party header files.
#include <catch.hpp>
#include <cmath>

// MORIS project header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"
#include "fn_qr.hpp"

#include "fn_print.hpp"
// ----------------------------------------------------------------------------

TEST_CASE(
        "moris::qr",
        "[linalgebra],[qr]" )
{

    SECTION( "qr(n-by-n mat)" )
    {
        moris::Matrix< moris::real,  moris::DDRMat>  tA = { { 5, 3, 4 } ,  { 1, 6, 3 } , { 1, 2, 4 } };
        moris::Matrix< moris::real,  moris::DDRMat>  tQ;
        moris::Matrix< moris::real,  moris::DDRMat>  tR;

        moris::qr( tQ, tR, tA );

        moris::Matrix< moris::real,  moris::DDRMat>  tB = tA - tQ*tR;

        moris::real tEpsilon = 1e-12;

        for( moris::uint j=0; j<tB.n_cols(); ++j )
        {
            for( moris::uint i=0; i<tB.n_rows(); ++i )
            {
                REQUIRE( std::abs( tB( i, j ) ) < tEpsilon );
            }
        }

    }
}
