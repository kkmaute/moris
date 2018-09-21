#include <catch.hpp>
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_save_matrix_to_binary_file.hpp" //LNA/src
#include "fn_load_matrix_from_binary_file.hpp" //LNA/src
#include "op_times.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "fn_norm.hpp"

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange QUAD8", "[moris],[fem]" )
{

//------------------------------------------------------------------------------

        // step 1: load MATLAB precomputed data from binary files
        std::string tPrefix = std::getenv("MORISROOT");
        tPrefix = tPrefix + "/projects/FEM/INT/test/data/" ;

        // load point coordinates from file
        Matrix< DDRMat > tXi;
        load_matrix_from_binary_file( tXi,
                    "./../moris/projects/FEM/INT/test/data/points_2d.bin" );


        // load values from nodes from file
        Matrix< DDRMat > tPhiHat;
        load_matrix_from_binary_file( tPhiHat,
                "./../moris/projects/FEM/INT/test/data/lagrange_quad8_phihat.bin" );


        // load solutions for N*tPhiHat
        Matrix< DDRMat > tPhi;
           load_matrix_from_binary_file( tPhi,
                   "./../moris/projects/FEM/INT/test/data/lagrange_quad8_phi.bin" );

        // load solutions for dNdXi*tPhiHat
        Matrix< DDRMat > tdPhidXi;
        load_matrix_from_binary_file( tdPhidXi,
                   "./../moris/projects/FEM/INT/test/data/lagrange_quad8_dphidxi.bin" );

        // load solutions for d2NdXi2*tPhiHat
        Matrix< DDRMat > td2PhidXi2;
        load_matrix_from_binary_file( td2PhidXi2,
                "./../moris/projects/FEM/INT/test/data/lagrange_quad8_d2phidxi2.bin" );

//------------------------------------------------------------------------------

        // step 2 create function and interpolation matrices

        // create rule
        Interpolation_Rule tRule(
                mtk::Geometry_Type::QUAD,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::SERENDIPITY  );

        // create shape function object
        auto tFunction = tRule.create_space_time_interpolation_function();

        // create matrix that contains the shape function
        auto tN        = tFunction->create_matrix( 1, 0, 0 );

        // create matrix that contains the first derivative
        auto tdNdXi    = tFunction->create_matrix( 1, 1, 0 );

        // create matrix that contains the second derivative
        auto td2NdXi2  = tFunction->create_matrix( 1, 2, 0 );

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-12;

        // get number of points to test
        auto tNumberOfTestPoints = tXi.n_cols();

//------------------------------------------------------------------------------

        SECTION( "QUAD8: test for unity" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tN, tXi.get_column(k ) );

                // test unity
                tCheck = tCheck && ( std::abs( tN.sum() - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "QUAD8: test N" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tN, tXi.get_column(k ) );

                // test evaluated value
                Matrix< DDRMat > tError = tN*tPhiHat - tPhi.get_row( k );

                // test error
                tCheck = tCheck && ( norm(tError) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "QUAD8: test dNdXi" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_dNdXi( tdNdXi, tXi.get_column(k ) );

                // test evaluated value
                Matrix< DDRMat > tError = tdNdXi*tPhiHat- tdPhidXi.get_column( k );

                // test error
                tCheck = tCheck && ( norm(tError) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "QUAD8: test d2NdXi2" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_d2NdXi2( td2NdXi2, tXi.get_column(k ) );

                // test evaluated valueN

                Matrix< DDRMat > tError = td2NdXi2*tPhiHat- td2PhidXi2.get_column( k );

                // test error
                tCheck = tCheck && ( norm(tError) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        // tidy up
        delete tFunction;

//------------------------------------------------------------------------------
}
