#include <catch.hpp>
// MRS/COR/src
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
// LINALG/src
#include "linalg_typedefs.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_sum.hpp"
// FEM/INT/src
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Integration_Rule.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "fn_FEM_Check.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "Lagrange TRI3", "[moris],[fem],[Tri3LagInterp]" )
{

    //------------------------------------------------------------------------------

    // create an interpolation rule - space only lagrange cubic triangle TRI3
    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction =
            tRule.create_space_interpolation_function();

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_1,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tZeta.get_column( k ), tN );

            // test unity
            tCheckPU = tCheckPU && ( std::abs( sum( tN ) - 1.0 ) < tEpsilon );
        }

        REQUIRE( tCheckPU );
    }

    SECTION( "TRI3: test dNdXi" )
    {
        // boolean to check evaluated dNdZeta
        bool tCheckdNdXi = true;

        // create matrix that contains the first order derivatives
        Matrix< DDRMat > tdNdXi;

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );


                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
        }
        REQUIRE( tCheck );
    }

    //------------------------------------------------------------------------------
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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
    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction =
            tRule.create_space_interpolation_function();

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_3,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );


                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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
    mtk::Interpolation_Rule tRule( mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_6,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tZeta.get_column( k ), tN );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );

                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 2 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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
    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_1,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tZeta.get_column( k ), tN );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );

                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 3 )
                {
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 3 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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
    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_4,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );

                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 3 )
                {
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 3 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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
    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create shape function object
    mtk::Interpolation_Function_Base* tFunction = tRule.create_space_interpolation_function();

    // create matrix that contains the second derivative
    Matrix< DDRMat > td2NdXi2;

    //------------------------------------------------------------------------------

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation
    real tPerturbation = 1e-6;

    // use the integration points as test points
    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_5,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // evaluate shape function at point k
            tFunction->eval_N( tZeta.get_column( k ), tN );

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

        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the first order derivative dNdZeta at test point k
            tFunction->eval_dNdXi( tTestPoint, tdNdXi );

            Matrix< DDRMat > tdNdXiFD( tdNdXi.n_rows(), tdNdXi.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tTestPoint.numel() - 1; iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tNMinus;
                tFunction->eval_N( tPertEvalPoint, tNMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tNPlus;
                tFunction->eval_N( tPertEvalPoint, tNPlus );

                // compute the first order derivatives wrt param coords by finite difference
                tdNdXiFD.get_row( iDim ) = ( tNPlus - tNMinus ) / ( 2.0 * tPerturbation );
            }
            // check evaluated derivatives against FD
            tCheckdNdXi = tCheckdNdXi && fem::check( tdNdXi, tdNdXiFD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfTestPoints; ++k )
        {
            // unpack the test point k
            Matrix< DDRMat > tTestPoint = tZeta.get_column( k );

            // evaluation of the second order derivatives d2NdZeta2 at test point k
            tFunction->eval_d2NdXi2( tTestPoint, td2NdXi2 );

            Matrix< DDRMat > td2NdXi2FD( td2NdXi2.n_rows(), td2NdXi2.n_cols(), 0.0 );

            for ( uint iDim = 0; iDim < tFunction->get_number_of_param_dimensions(); iDim++ )
            {
                // perturbed evaluation point
                Matrix< DDRMat > tPertEvalPoint = tTestPoint;
                tPertEvalPoint( iDim )          = tPertEvalPoint( iDim ) - tPerturbation;

                // evaluate N at point l - delta xi
                Matrix< DDRMat > tdNdXiMinus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiMinus );

                // perturbed evaluation point
                tPertEvalPoint         = tTestPoint;
                tPertEvalPoint( iDim ) = tPertEvalPoint( iDim ) + tPerturbation;

                // evaluate td2NdXi2 at point k + delta xi
                Matrix< DDRMat > tdNdXiPlus;
                tFunction->eval_dNdXi( tPertEvalPoint, tdNdXiPlus );

                // compute the second order derivatives by finite difference
                Matrix< DDRMat > td2NdXi2FDTemp = ( tdNdXiPlus - tdNdXiMinus ) / ( 2.0 * tPerturbation );

                if ( iDim == 0 )
                {
                    td2NdXi2FD.get_row( 0 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 1 )
                {
                    td2NdXi2FD.get_row( 9 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 1 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 2 )
                {
                    td2NdXi2FD.get_row( 8 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 7 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 2 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 3 );
                }
                else if ( iDim == 3 )
                {
                    td2NdXi2FD.get_row( 6 ) = td2NdXi2FDTemp.get_row( 0 );
                    td2NdXi2FD.get_row( 5 ) = td2NdXi2FDTemp.get_row( 1 );
                    td2NdXi2FD.get_row( 4 ) = td2NdXi2FDTemp.get_row( 2 );
                    td2NdXi2FD.get_row( 3 ) = td2NdXi2FDTemp.get_row( 3 );
                }
            }
            // check
            tCheck = tCheck && fem::check( td2NdXi2, td2NdXi2FD, tEpsilon );
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
        for ( uint k = 0; k < tNumOfParamPoints; ++k )
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

TEST_CASE( "Lagrange TET4 integration", "[moris],[fem],[Tet4LagInteg]" )
{
    //------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET4 in the physical space
    // clang-format off
    Matrix< DDRMat > tXHat = {
            { 0.0,  0.0, 0.0 },
            { 0.0, -1.0, 0.0 },
            { 1.0,  0.0, 0.0 },
            { 0.0,  0.0, 1.0 }};
    // clang-format on

    Matrix< DDRMat > tTHat = { { 0.0 }, { 2.0 } };

    real tExpectedVolume = 2.0 * 1.0 / 6.0;

    // define an interpolation rule for the TET4
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_5,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
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
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
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
    real t12 = 1.0 / 2.0;

    // clang-format off
    Matrix< DDRMat > tXHat = {
                { 0.0,  0.0, 0.0 },
                { 0.0, -1.0, 0.0 },
                { 1.0,  0.0, 0.0 },
                { 0.0,  0.0, 1.0 },
                { 0.0, -t12, 0.0 },
                { t12, -t12, 0.0 },
                { t12,  0.0, 0.0 },
                { 0.0,  0.0, t12 },
                { 0.0, -t12, t12 },
                { t12,  0.0, t12 } };
    // clang-format on

    Matrix< DDRMat > tTHat = { { 0.0 }, { 2.0 } };

    real tExpectedVolume = 2.0 * 0.5 / 3.0;

    // define an interpolation rule for the TET10
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry interpolator and set coefficients
    Geometry_Interpolator tGeoInterpolator( tGeomRule );
    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TET_4
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_5,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
    uint tNumOfIntegPoints = tIntegrator.get_number_of_points();

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );

    // init volume
    real tVolume = 0.0;
    real tWStar  = 0.0;

    // loop over integration points
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
    {
        // set integration point for geometry interpolator
        tGeoInterpolator.set_space_time( tIntegPoints.get_column( iGP ) );

        // compute integration point weight x detJ
        tWStar = tGeoInterpolator.det_J() * tIntegWeights( iGP );

        // add contribution to jacobian from evaluation point
        tVolume = tVolume + tWStar;
    }

    bool tCheck = true;
    tCheck      = tCheck && ( std::abs( tVolume - tExpectedVolume ) < tEpsilon );
    REQUIRE( tCheck );
}

TEST_CASE( "Lagrange TET20 integration", "[moris],[fem],[Tet20LagInteg]" )
{
    //------------------------------------------------------------------------------

    // define an epsilon environment
    double tEpsilon = 1E-12;

    // define a TET20 in the physical space
    real t13 = 1.0 / 3.0;
    real t23 = 2.0 / 3.0;

    // clang-format off
    Matrix< DDRMat > tXHat = {
                { 0.0,  0.0, 0.0 },
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
    // clang-format on

    Matrix< DDRMat > tTHat = { { 0.0 }, { 2.0 } };

    real tExpectedVolume = 2.0 * 0.5 / 3.0;

    // define an interpolation rule for the TET20
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TET,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TET,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TET_15,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
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
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
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
    // clang-format off
    Matrix< DDRMat > tXHat = {
            { 0.0,  0.0 },
            { 1.0, -1.0 },
            { 3.0,  0.0 }};
    // clang-format on

    Matrix< DDRMat > tTHat = { { 0.0 }, { 2.0 } };

    real tExpectedVolume = 3;

    // define an interpolation rule for the TRI3
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_1,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
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
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
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
    // clang-format off
    Matrix< DDRMat > tXHat = {
                { 0.0,  0.0 },
                { 1.0, -1.0 },
                { 3.0,  0.0 },
                { 0.5, -0.5 },
                { 2.0, -0.5 },
                { 1.5,  0.0 }};
    // clang-format on

    Matrix< DDRMat > tTHat = { { 0.0 }, { 2.0 } };

    real tExpectedVolume = 3;

    // define an interpolation rule for the TRI10
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI6
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_3,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
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
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
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
    // clang-format off
    Matrix< DDRMat > tXHat = {
                { 0.0,          0.0 },
                { 1.0,         -1.0 },
                { 3.0,          0.0 },
                { 1.0/3.0,     -1.0/3.0 },
                { 2.0/3.0,     -2.0/3.0 },
                { 1.0+2.0/3.0, -2.0/3.0 },
                { 1.0+4.0/3.0, -1.0/3.0 },
                { 2.0,          0.0 },
                { 1.0,          0.0 },
                { 4.0/3.0,     -0.5 }};
    // clang-format on

    Matrix< DDRMat > tTHat           = { { 0.0 }, { 2.0 } };
    real             tExpectedVolume = 3;

    // define an interpolation rule for the TRI10
    mtk::Interpolation_Rule tGeomRule(
            mtk::Geometry_Type::TRI,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::CUBIC,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create the element geometry intepolator
    Geometry_Interpolator tGeoInterpolator( tGeomRule );

    tGeoInterpolator.set_coeff( tXHat, tTHat );

    // create an integration rule - space only lagrange linear triangle TRI3
    mtk::Integration_Rule tIntegrationRule(
            mtk::Geometry_Type::TRI,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::TRI_6,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_1 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get number of integration points
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
    for ( uint iGP = 0; iGP < tNumOfIntegPoints; iGP++ )
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
