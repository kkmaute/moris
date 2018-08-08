#include <catch.hpp>
#include "cl_FEM_Interpolation_Matrix.hpp" //FEM/INT/src

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src
#include "fn_save_matrix_to_binary_file.hpp" //LNA/src
#include "fn_load_matrix_from_binary_file.hpp" //LNA/src
#include "op_times.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange HEX20", "[moris],[fem]" )
{

//------------------------------------------------------------------------------

        // step 1: load MATLAB precomputed data from binary files

        // load point coordinates from file
        Mat< real > tXi;
        load_matrix_from_binary_file( tXi,
                    "./../moris/projects/FEM/INT/test/data/points_3d.bin" );


        // load values from nodes from file
        Mat< real > tPhiHat;
        load_matrix_from_binary_file( tPhiHat,
                "./../moris/projects/FEM/INT/test/data/lagrange_hex20_phihat.bin" );


        // load solutions for N*tPhiHat
        Mat< real > tPhi;
           load_matrix_from_binary_file( tPhi,
                   "./../moris/projects/FEM/INT/test/data/lagrange_hex20_phi.bin" );

        // load solutions for dNdXi*tPhiHat
        Mat< real > tdPhidXi;
        load_matrix_from_binary_file( tdPhidXi,
                   "./../moris/projects/FEM/INT/test/data/lagrange_hex20_dphidxi.bin" );

        // load solutions for d2NdXi2*tPhiHat
        Mat< real > td2PhidXi2;
        load_matrix_from_binary_file( td2PhidXi2,
                "./../moris/projects/FEM/INT/test/data/lagrange_hex20_d2phidxi2.bin" );

//------------------------------------------------------------------------------

        // step 2 create function and interpolation matrices

        // create rule
        Interpolation_Rule tRule(
                mtk::Geometry_Type::HEX,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::SERENDIPITY  );

        // create shape function object
        auto tFunction = tRule.create_space_time_interpolation_function();

        // create matrix that contains the shape function
        auto tN        = tFunction->create_matrix( 1, 0, 0, 1 );

        // create matrix that contains the first derivative
        auto tdNdXi    = tFunction->create_matrix( 1, 1, 0, 1 );

        // create matrix that contains the second derivative
        auto td2NdXi2  = tFunction->create_matrix( 1, 2, 0, 1 );

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-12;

        // get number of points to test
        auto tNumberOfTestPoints = tXi.n_cols();

//------------------------------------------------------------------------------

        SECTION( "HEX20: test for unity" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tN, tXi.cols( k,k ) );

                // test unity
                tCheck = tCheck && ( std::abs( tN.sum() - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "HEX20: test N" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tN, tXi.cols( k,k ) );

                // test evaluated value
                Mat< real > tError = tN*tPhiHat - tPhi( k );

                // test error
                tCheck = tCheck && ( tError.norm() < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "HEX20: test dNdXi" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_dNdXi( tdNdXi, tXi.cols( k,k ) );

                // test evaluated value
                Mat< real > tError = tdNdXi*tPhiHat- tdPhidXi.cols( k, k );

                // test error
                tCheck = tCheck && ( tError.norm() < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        SECTION( "HEX20: test d2NdXi2" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_d2NdXi2( td2NdXi2, tXi.cols( k,k ) );

                // test evaluated value

                Mat< real > tError = td2NdXi2*tPhiHat - td2PhidXi2.cols( k, k );

                // test error
                tCheck = tCheck && ( tError.norm() < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        // tidy up
        delete tFunction;

//------------------------------------------------------------------------------
}
