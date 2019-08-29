#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_times.hpp" //LNA/src
#include "op_minus.hpp" //LNA/src
#include "fn_trans.hpp" //LNA/src
#include "fn_norm.hpp" //LNA/src
#include "fn_sum.hpp" //LNA/src

#include "cl_FEM_Interpolation_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Integrator.hpp" //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp" //FEM/INT/src

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange TRI3", "[moris],[fem],[Tri3LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TRI3
        Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TRI_1,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TRI3: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ), tN );

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TRI3: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;

                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint, tN );
                    tFunction->eval_N( tTestPointPert, tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TRI3: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the second order derivatives d2NdZeta2 at test point k
                tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TRI3: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint, tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }

TEST_CASE( "Lagrange TRI6", "[moris],[fem],[Tri6LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TRI6
        Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::QUADRATIC,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TRI_3,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TRI6: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ), tN );

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TRI6: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;

                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint,tN );
                    tFunction->eval_N( tTestPointPert,tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TRI6: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                 // evaluation of the second order derivatives d2NdZeta2 at test point k
                 tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TRI6: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint, tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }


TEST_CASE( "Lagrange TRI10", "[moris],[fem],[Tri10LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TRI10
        Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::CUBIC,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TRI_6,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TRI10: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ),tN );

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TRI10: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;

                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint, tN );
                    tFunction->eval_N( tTestPointPert ,tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TRI10: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                 // evaluation of the second order derivatives d2NdZeta2 at test point k
                 tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TRI10: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint, tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }

TEST_CASE( "Lagrange TET4", "[moris],[fem],[Tet4LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TET4
        Interpolation_Rule tRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TET_1,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TET4: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ),tN );

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TET4: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;
                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint, tN );
                    tFunction->eval_N( tTestPointPert, tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TET4: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                 // evaluation of the second order derivatives d2NdZeta2 at test point k
                 tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 8 , iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 8, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 3 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TET4: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint, tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }

TEST_CASE( "Lagrange TET10", "[moris],[fem],[Tet10LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TET10
        Interpolation_Rule tRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::QUADRATIC,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TET_4,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TET10: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ), tN );

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TET10: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;

                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint, tN );
                    tFunction->eval_N( tTestPointPert,tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TET10: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                 // evaluation of the second order derivatives d2NdZeta2 at test point k
                 tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 8 , iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 8, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 3 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TET10: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint, tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }


TEST_CASE( "Lagrange TET20", "[moris],[fem],[Tet20LagInterp]" )
{

//------------------------------------------------------------------------------

        // create an interpolation rule - space only lagrange cubic triangle TET20
        Interpolation_Rule tRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::CUBIC,
                                  Interpolation_Type::CONSTANT,
                                  mtk::Interpolation_Order::CONSTANT );

        // create shape function object
        Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

        // create matrix that contains the second derivative
        Matrix< DDRMat > td2NdXi2;

//------------------------------------------------------------------------------

        // define an epsilon environment
        double tEpsilon = 1E-4;

        // use the integration points as test points
        // create an integration rule
        Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                           Integration_Type::GAUSS,
                                           Integration_Order::TET_5,
                                           Integration_Type::GAUSS,
                                           Integration_Order:: BAR_1 );

        // create an integrator
        Integrator tIntegrator( tIntegrationRule );

        // get integration points
        Matrix< DDRMat > tZeta;
        tIntegrator.get_points( tZeta );

        // get number of points to test
        uint tNumOfTestPoints = tZeta.n_cols();

//------------------------------------------------------------------------------

        SECTION( "TET20: test for unity" )
        {
            bool tCheckPU = true;

            // create matrix that contains the shape function
            Matrix< DDRMat > tN;

            for( uint k=0; k<tNumOfTestPoints; ++k )
            {
                // evaluate shape function at point k
                tFunction->eval_N( tZeta.get_column( k ) , tN);

                // test unity
                tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
            }

            REQUIRE( tCheckPU );
        }

//------------------------------------------------------------------------------

        SECTION( "TET20: test dNdXi" )
        {

            // boolean to check evaluated dNdZeta
            bool tCheckdNdXi = true;

            // create matrix that contains the first order derivatives
            Matrix< DDRMat > tdNdXi;

            for( uint k=0; k < tNumOfTestPoints; ++k )
            {
                // unpack the test point k
                Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                // evaluation of the first order derivative dNdZeta at test point k
                tFunction->eval_dNdXi( tTestPoint, tdNdXi );

                for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                {
                    // perturbation of the test point
                    Matrix< DDRMat > tTestPointPert = tTestPoint;
                    real tPert = 1E-6 * tTestPointPert( iDim );
                    tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                    Matrix< DDRMat > tN;
                    Matrix< DDRMat > tNPert;

                    // evaluate shape functions at test point and perturbed test point
                    tFunction->eval_N( tTestPoint, tN );
                    tFunction->eval_N( tTestPointPert, tNPert );

                    // compute the first order derivatives wrt param coords by finite difference
                    Matrix< DDRMat > tdNdXi_FD = ( tNPert - tN ) / tPert;

                    // check evaluated derivatives against FD
                    for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                    {
                        tCheckdNdXi = tCheckdNdXi && ( std::abs( tdNdXi_FD( iBase ) - tdNdXi( iDim, iBase ) ) < tEpsilon );
                    }
                }
            }
            REQUIRE( tCheckdNdXi );
        }

//------------------------------------------------------------------------------

SECTION( "TET20: test d2NdXi2" )
{
             // boolean to check evaluated d2NdZeta2
             bool tCheck = true;

             // create matrix that contains the second order derivatives
             Matrix< DDRMat > td2NdXi2;

             // loop over the test points
             for( uint k=0; k < tNumOfTestPoints; ++k )
             {
                 // unpack the test point k
                 Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

                 // evaluation of the second order derivatives d2NdZeta2 at test point k
                 tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

                 for (  uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
                 {
                     // perturbation of the test point
                     Matrix< DDRMat > tTestPointPert = tTestPoint;
                     real tPert = 1E-6 * tTestPointPert( iDim );
                     tTestPointPert( iDim ) = tTestPointPert( iDim ) + tPert;

                     // evaluate the first derivatives of the shape functions at test point and perturbed test point
                     Matrix< DDRMat > tdNdXi;
                     tFunction->eval_dNdXi( tTestPoint, tdNdXi );
                     Matrix< DDRMat > tdNdXiPert;
                     tFunction->eval_dNdXi( tTestPointPert, tdNdXiPert );

                     // compute the second order derivatives by finite difference
                     Matrix< DDRMat > td2NdXi2_FD = ( tdNdXiPert - tdNdXi ) / tPert;

                     // check evaluated derivatives against FD
                     for ( uint iBase = 0; iBase < tFunction->get_number_of_bases(); iBase++ )
                     {
                         if ( iDim == 0 )
                         {
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 0, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 8 , iBase ) ) < tEpsilon );
                            tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 1 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 9, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 1, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 2 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 8, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 7, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 2, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                         }
                         else if( iDim == 3 )
                         {
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 0, iBase ) - td2NdXi2( 6, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 1, iBase ) - td2NdXi2( 5, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 2, iBase ) - td2NdXi2( 4, iBase ) ) < tEpsilon );
                             tCheck = tCheck && ( std::abs( td2NdXi2_FD( 3, iBase ) - td2NdXi2( 3, iBase ) ) < tEpsilon );
                         }
                     }
                 }
             }
             REQUIRE( tCheck );
}

SECTION( "TET20: test param coords" )
{

    // get the param points
    Matrix< DDRMat > tZetaCoords;
    tFunction->get_param_coords( tZetaCoords );

    // number of param points
    uint tNumOfParamPoints = tZetaCoords.n_cols();

    // create matrix that contains the second order derivatives
    Matrix< DDRMat > tN;

    // boolean to check param point
    bool tCheck = true;

    // loop over the param points
    for( uint k=0; k < tNumOfParamPoints; ++k )
    {
        // get the param point k
        Matrix< DDRMat > tParamPoint = tZetaCoords.get_column( k );

        // evaluate shape functions at param point k
        tFunction->eval_N( tParamPoint,tN );

        // check that kth shape function = 1
        tCheck = tCheck && ( std::abs( tN( k ) - 1.0 ) < tEpsilon );
    }
    REQUIRE( tCheck );
}

 //------------------------------------------------------------------------------

         // tidy up
         delete tFunction;

 //------------------------------------------------------------------------------
 }

TEST_CASE( "Lagrange TET4 integration", "[moris],[fem],[Tet4LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET4 in the physical space
    Matrix< DDRMat > tXHat = {{ 0.0,  0.0, 0.0 },
                              { 0.0, -1.0, 0.0 },
                              { 1.0,  0.0, 0.0 },
                              { 0.0,  0.0, 1.0 }};
    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 2.0 * 1.0 / 6.0;

    // define an interpolation rule for the TET4
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TET_5,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        real tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }
    REQUIRE( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
}

TEST_CASE( "Lagrange TET10 integration", "[moris],[fem],[Tet10LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET10 in the physical space
    real t12 = 1.0/2.0;
    Matrix< DDRMat > tXHat = {{ 0.0,  0.0, 0.0 },
                              { 0.0, -1.0, 0.0 },
                              { 1.0,  0.0, 0.0 },
                              { 0.0,  0.0, 1.0 },
                              { 0.0, -t12, 0.0 },
                              { t12, -t12, 0.0 },
                              { t12,  0.0, 0.0 },
                              { 0.0,  0.0, t12 },
                              { 0.0, -t12, t12 },
                              { t12,  0.0, t12 } };

    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 2 * 0.5 / 3.0;

    // define an interpolation rule for the TET10
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::QUADRATIC,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry interpolator and set coefficients
    Geometry_Interpolator tGeoInterpolator( tGeomRule );
    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TET_4
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TET_5,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0.0;
    real tWStar = 0.0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }

    bool tCheck = true;
    tCheck = tCheck && ( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
    REQUIRE( tCheck );
}

TEST_CASE( "Lagrange TET20 integration", "[moris],[fem],[Tet20LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET20 in the physical space
    real t13 = 1.0/3.0;
    real t23 = 2.0/3.0;
    Matrix< DDRMat > tXHat = {{ 0.0,  0.0, 0.0 },
                              { 0.0, -1.0, 0.0 },
                              { 1.0,  0.0, 0.0 },
                              { 0.0,  0.0, 1.0 },
                              { 0.0, -t13, 0.0 },
                              { 0.0, -t23, 0.0 },
                              { t13, -t23, 0.0 },
                              { t23, -t13, 0.0 },
                              { t13,  0.0, 0.0 },
                              { t23,  0.0, 0.0 },
                              { 0.0,  0.0, t13 },
                              { 0.0,  0.0, t23 },
                              { 0.0, -t23, t13 },
                              { 0.0, -t13, t23 },
                              { t23,  0.0, t13 },
                              { t13,  0.0, t23 },
                              { t13, -t13, 0.0 },
                              { 0.0, -t13, t13 },
                              { t13, -t13, t13 },
                              { t13,  0.0, t13 } };

    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 2 * 0.5 / 3.0;

    // define an interpolation rule for the TET20
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TET,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::CUBIC,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TET,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TET_15,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        real tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }
    REQUIRE( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
}

TEST_CASE( "Lagrange TRI3 integration", "[moris],[fem],[Tri3LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET4 in the physical space
    Matrix< DDRMat > tXHat = {{ 0.0,  0.0 },
                              { 1.0, -1.0 },
                              { 3.0,  0.0 }};
    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 3;

    // define an interpolation rule for the TRI3
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TRI_1,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
	tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        real tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }
    REQUIRE( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
}


TEST_CASE( "Lagrange TRI6 integration", "[moris],[fem],[Tri6LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET4 in the physical space
    Matrix< DDRMat > tXHat = {{ 0.0,  0.0 },
                              { 1.0, -1.0 },
                              { 3.0,  0.0 },
                              { 0.5, -0.5 },
                              { 2.0, -0.5 },
                              { 1.5,  0.0 }};
    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 3;

    // define an interpolation rule for the TRI10
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::QUADRATIC,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI6
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TRI_3,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        real tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }
    REQUIRE( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
}

TEST_CASE( "Lagrange TRI10 integration", "[moris],[fem],[Tri10LagInteg]" )
{
//------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET4 in the physical space
    Matrix< DDRMat > tXHat = {{ 0.0,          0.0 },
                              { 1.0,         -1.0 },
                              { 3.0,          0.0 },
                              { 1.0/3.0,     -1.0/3.0 },
                              { 2.0/3.0,     -2.0/3.0 },
                              { 1.0+2.0/3.0, -2.0/3.0 },
                              { 1.0+4.0/3.0, -1.0/3.0 },
                              { 2.0,          0.0 },
                              { 1.0,          0.0 },
                              { 4.0/3.0,     -0.5 }};

    Matrix< DDRMat > tTHat = {{0.0}, {2.0}};
    real tExpectedVolume = 3;

    // define an interpolation rule for the TRI10
    Interpolation_Rule tGeomRule( mtk::Geometry_Type::TRI,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::CUBIC,
                                  Interpolation_Type::LAGRANGE,
                                  mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    Integration_Rule tIntegrationRule( mtk::Geometry_Type::TRI,
                                       Integration_Type::GAUSS,
                                       Integration_Order::TRI_6,
                                       Integration_Type::GAUSS,
                                       Integration_Order:: BAR_1 );

    // create an integrator
    Integrator tIntegrator( tIntegrationRule );

    //get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0;

    // loop over integration points
    for( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        real tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }
    REQUIRE( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
}

