// Third-party header files.
#include <catch.hpp>
#include <cmath>

// MORIS project header files.
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "typedefs.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"
#include "fn_lu.hpp"

#include "fn_print.hpp"
// ----------------------------------------------------------------------------

TEST_CASE("moris::lu", "[linalgebra],[lu]" )
{

    SECTION( "lu(n-by-n mat)" )
    {
        moris::Matrix< moris::real,  moris::DDRMat>  tA = { { 5, 3, 4 } ,  { 1, 6, 3 } , { 1, 2, 4 } };
        moris::Matrix< moris::real,  moris::DDRMat>  tL;
        moris::Matrix< moris::real,  moris::DDRMat>  tU;
        moris::Matrix< moris::real,  moris::DDRMat>  tP;

        moris::lu( tL, tU, tP, tA );

        moris::Matrix< moris::real,  moris::DDRMat>  tB = tP*tA - tL*tU;

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
