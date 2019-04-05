#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_save_matrix_to_binary_file.hpp" //LNA/src
#include "fn_load_matrix_from_binary_file.hpp" //LNA/src
#include "op_times.hpp" //LNA/src
#include "op_minus.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "fn_norm.hpp" //LNA/src
#include "fn_sum.hpp" //LNA/src

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange TRI3", "[moris],[fem],[Tri3LagInterpolation]" )
{

//------------------------------------------------------------------------------

        // step 2 create function and interpolation matrices

        // create rule
        Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        auto tFunction = tRule.create_space_interpolation_function();

        // create matrix that contains the shape function
        Matrix< DDRMat > tN;

        // create matrix that contains the first derivative
        Matrix< DDRMat > tdNdXi;

        // create matrix that contains the second derivative
        Matrix< DDRMat > td2NdXi2;

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // test poiny
        Matrix< DDRMat > tXi = {{ 0.5 }, { 0.6 }};

        // get number of points to test
        uint tNumberOfTestPoints = tXi.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TRI6: test for unity" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tN = tFunction->eval_N( tXi.get_column( k ) );

                // test unity
                tCheck = tCheck && ( std::abs( sum(tN) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        tdNdXi = tFunction->eval_dNdXi( tXi );
        td2NdXi2 = tFunction->eval_d2NdXi2( tXi );

//------------------------------------------------------------------------------

        SECTION( "TRI6: test dNdXi" )
        {

            bool tCheck = true;

            for( uint k=0; k < tNumberOfTestPoints; ++k )
            {
                Matrix< DDRMat > tXiTreated = tXi.get_column(k );
                tdNdXi = tFunction->eval_dNdXi( tXiTreated );

                for (  uint iDim = 0; iDim < 2; iDim++ )
                {
                    Matrix< DDRMat > tXiPert = tXiTreated;
                    real tPert = 1E-6 * tXiPert( iDim );
                    tXiPert( iDim ) = tXiPert( iDim ) + tPert;

                    // evaluate shape function at point k
                    tN = tFunction->eval_N( tXiTreated );
                    Matrix< DDRMat > tNPert = tFunction->eval_N( tXiPert );

                    Matrix< DDRMat > tdNdXi_iDim = ( tNPert - tN ) / tPert;

                    for ( uint iBase = 0; iBase < 3; iBase++ )
                    {
                        tCheck = tCheck && ( std::abs( tdNdXi_iDim( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

//        SECTION( "TRI6: test d2NdXi2" )
//        {
//            bool tCheck = true;
//
//            for( uint k=0; k < tNumberOfTestPoints; ++k )
//            {
//                Matrix< DDRMat > tXiTreated = tXi.get_column(k );
//                td2NdXi2 = tFunction->eval_d2NdXi2( tXiTreated );
//
//                for (  uint iDim = 0; iDim < 2; iDim++ )
//                {
//                    Matrix< DDRMat > tXiPert = tXiTreated;
//                    real tPert = 1E-6 * tXiPert( iDim );
//                    tXiPert( iDim ) = tXiPert( iDim ) + tPert;
//
//                    // evaluate shape function at point k
//                    tdNdXi = tFunction->eval_dNdXi( tXiTreated );
//                    Matrix< DDRMat > tdNdXiPert = tFunction->eval_dNdXi( tXiPert );
//
//                    Matrix< DDRMat > td2NdXi2_iDim = ( tdNdXiPert - tdNdXi ) / tPert;
//
//                    for ( uint iBase = 0; iBase < 3; iBase++ )
//                    {
//                        tCheck = tCheck && ( std::abs( td2NdXi2_iDim( iDim, iBase ) - td2NdXi2( iDim, iBase ) ) < tEpsilon );
//                        std::cout<<td2NdXi2_iDim( iDim, iBase )<<std::endl;
//                        std::cout<<td2NdXi2( iDim, iBase )<<std::endl;
//
//                        tCheck = tCheck && ( std::abs( td2NdXi2_iDim( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
//                        std::cout<<td2NdXi2_iDim( 2, iBase )<<std::endl;
//                        std::cout<<td2NdXi2( 2, iBase )<<std::endl;
//                    }
//                }
//            }
//            REQUIRE( tCheck );
//        }

//------------------------------------------------------------------------------

        // tidy up
        delete tFunction;

//------------------------------------------------------------------------------
}

TEST_CASE( "Lagrange TRI6", "[moris],[fem],[Tri6LagInterpolation]" )
{

//------------------------------------------------------------------------------

        // step 2 create function and interpolation matrices

        // create rule
        Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::QUADRATIC,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        auto tFunction = tRule.create_space_interpolation_function();

        // create matrix that contains the shape function
        Matrix< DDRMat > tN;

        // create matrix that contains the first derivative
        Matrix< DDRMat > tdNdXi;

        // create matrix that contains the second derivative
        Matrix< DDRMat > td2NdXi2;

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // test poiny
        Matrix< DDRMat > tXi = {{ 0.5 }, { 0.6 }};

        // get number of points to test
        uint tNumberOfTestPoints = tXi.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TRI6: test for unity" )
        {
            bool tCheck = true;
            for( uint k=0; k<tNumberOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tN = tFunction->eval_N( tXi.get_column( k ) );

                // test unity
                tCheck = tCheck && ( std::abs( sum(tN) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

        tdNdXi = tFunction->eval_dNdXi( tXi );
        td2NdXi2 = tFunction->eval_d2NdXi2( tXi );

//------------------------------------------------------------------------------

        SECTION( "TRI6: test dNdXi" )
        {


            bool tCheck = true;
            for( uint k=0; k < tNumberOfTestPoints; ++k )
            {
            	Matrix< DDRMat > tXiTreated = tXi.get_column(k );
            	tdNdXi = tFunction->eval_dNdXi( tXiTreated );

            	for (  uint iDim = 0; iDim < 2; iDim++ )
            	{

            		Matrix< DDRMat > tXiPert = tXiTreated;
            		real tPert = 1E-6 * tXiPert( iDim );
            		tXiPert( iDim ) = tXiPert( iDim ) + tPert;

                    // evaluate shape function at point k
                    tN = tFunction->eval_N( tXiTreated );
                    Matrix< DDRMat > tNPert = tFunction->eval_N( tXiPert );

                    Matrix< DDRMat > tdNdXi_iDim = ( tNPert - tN ) / tPert;

                    for ( uint iBase = 0; iBase < 6; iBase++ )
                    {
                    	tCheck = tCheck && ( std::abs( tdNdXi_iDim( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );

                    	std::cout<<std::abs( tdNdXi_iDim( iBase ) - tdNdXi( iDim, iBase ) )<<std::endl;
                    	//std::cout<<tdNdXi_iDim( iBase )<<std::endl;
                    	//std::cout<<tdNdXi( iDim, iBase )<<std::endl;
                    }
            	}
            }

            REQUIRE( tCheck );
        }

//------------------------------------------------------------------------------

//        SECTION( "QUAD4: test d2NdXi2" )
//        {
//            bool tCheck = true;
//            for( uint k=0; k<tNumberOfTestPoints; ++k )
//            {
//                // evaluate shape function at point k
//            	td2NdXi2 = tFunction->eval_d2NdXi2( tXi.get_column(k ) );
//
//                // test evaluated valueN
//                Matrix< DDRMat > tError = td2PhidXi2.get_column( k );
//                tError = tError - td2NdXi2*tPhiHat;
//
//                // test error
//                tCheck = tCheck && ( norm(tError) < tEpsilon );
//            }
//
//            REQUIRE( tCheck );
//        }

//------------------------------------------------------------------------------

        // tidy up
        delete tFunction;

//------------------------------------------------------------------------------
}
