#include <string>
#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"              //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp" //FEM//INT//src
#include "cl_FEM_IWG_Test.hpp"           //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG", "[moris],[fem],[IWG]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::LINE,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    //create space coeff xHat
    Matrix< DDRMat > tXHat( 2, 1 );
    tXHat( 0, 0 ) = -1.0;
    tXHat( 1, 0 ) =  3.0;

    //create time coeff tHat
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat );

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::LINE,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    //set the number of interpolated fields
    uint tNumberOfFields = 1;

    //create a field interpolator
    Field_Interpolator * tFieldInterpolator = new Field_Interpolator( tNumberOfFields,
                                                                      tInterpolationRule,
                                                                      tGeomInterpolator );

    //get the number of basis for space time
    uint tNSpaceTimeBases = tFieldInterpolator->get_number_of_space_time_bases();

    //create field coeff uHat
    Matrix< DDRMat > tUHat( tNSpaceTimeBases , tNumberOfFields );
    tUHat( 0 ) = -1.0;
    tUHat( 1 ) =  7.0;
    tUHat( 2 ) = -6.0;
    tUHat( 3 ) = 82.0;

    //set the coefficients uHat
    tFieldInterpolator->set_coeff( tUHat );

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi( 1, 1 );
    tXi( 0, 0 ) =  0.35;
    Matrix< DDRMat > tTau( 1, 1 );
    tTau( 0, 0 ) = 0.70;
    Matrix< DDRMat > tParamPoint = { { tXi( 0 ) },
                                     { tTau( 0 ) } };

    //set the evaluation point xi, tau
    tFieldInterpolator->set_space_time( tParamPoint );

    // IWG
    //------------------------------------------------------------------------------
    // create an IWG
    IWG_Test tIWGField( tFieldInterpolator );

    // define an epsilon environment for comparison
    double tEpsilon = 1E-12;

    //------------------------------------------------------------------------------
    SECTION( "IWG : check IWG computed values - r = delta v * u" )
    {
        // residual from matlab
        Matrix< DDRMat > tResidualMatlab = {
        {2.244937500000001},
        {4.662562500000002},
        {12.721312500000002},
        {26.421187500000002}};

        // jacobian from matlab
        Matrix< DDRMat > tJacobianMatlab ={
        {0.002376562500000,   0.004935937500000,   0.013467187500000,   0.027970312500000},
        {0.004935937500000,   0.010251562500000,   0.027970312500000,   0.058092187500000},
        {0.013467187500000,   0.027970312500000,   0.076314062500000,   0.158498437500000},
        {0.027970312500000,   0.058092187500000,   0.158498437500000,   0.329189062500000}};

        // check independent evaluation of residual and jacobian

        //set the space and time derivative order for the field u
        uint tSpaceDerivativeOrder = 0;
        uint tTimeDerivativeOrder  = 0;

        // evaluate the residual from IWG
        Matrix< DDRMat > tResidual; tIWGField.compute_residual( tResidual,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // evaluate the jacobian from IWG
        Matrix< DDRMat > tJacobian; tIWGField.compute_jacobian( tJacobian,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // comparison value per value with matlab results
        bool tCheckJacobian = true; bool tCheckResidual = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual = tCheckResidual && ( ( tResidual( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian = tCheckJacobian && ( ( tJacobian( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual );
        REQUIRE( tCheckJacobian );

        // check simultaneous evaluation of residual and jacobian

        // evaluate the residual and the jacobian from IWG
        Matrix< DDRMat > tResidual2, tJacobian2;
        tIWGField.compute_jacobian_and_residual( tJacobian2,
                                                 tResidual2,
                                                 tSpaceDerivativeOrder,
                                                 tTimeDerivativeOrder);

        // comparison value per value with matlab results
        bool tCheckResidual2 = true; bool tCheckJacobian2 = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual2 = tCheckResidual2 && ( ( tResidual2( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian2 = tCheckJacobian2 && ( ( tJacobian2( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual2 );
        REQUIRE( tCheckJacobian2 );

    }

    //------------------------------------------------------------------------------
    SECTION( "IWG : check IWG computed values - r = delta v * dudx" )
    {
        // residual from matlab
        Matrix< DDRMat > tResidualMatlab = {
        {0.926250000000000},
        {1.923750000000000},
        {5.248750000000000},
        {10.901249999999999}};

        // jacobian from matlab
        Matrix< DDRMat > tJacobianMatlab ={
        {-0.001828125000000,   0.001828125000000,  -0.010359375000000,   0.010359375000000},
        {-0.003796875000000,   0.003796875000000,  -0.021515625000000,   0.021515625000000},
        {-0.010359375000000,   0.010359375000000,  -0.058703125000000,   0.058703125000000},
        {-0.021515625000000,   0.021515625000000,  -0.121921875000000,   0.121921875000000}};

        // check independent evaluation of residual and jacobian

        //set the space and time derivative order for the field u
        uint tSpaceDerivativeOrder = 1;
        uint tTimeDerivativeOrder  = 0;

        // evaluate the residual from IWG
        Matrix< DDRMat > tResidual; tIWGField.compute_residual( tResidual,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // evaluate the jacobian from IWG
        Matrix< DDRMat > tJacobian; tIWGField.compute_jacobian( tJacobian,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // comparison value per value with matlab results
        bool tCheckJacobian = true; bool tCheckResidual = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual = tCheckResidual && ( ( tResidual( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian = tCheckJacobian && ( ( tJacobian( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual );
        REQUIRE( tCheckJacobian );

        // check simultaneous evaluation of residual and jacobian

        // evaluate the residual and the jacobian from IWG
        Matrix< DDRMat > tResidual2, tJacobian2;
        tIWGField.compute_jacobian_and_residual( tJacobian2,
                                                 tResidual2,
                                                 tSpaceDerivativeOrder,
                                                 tTimeDerivativeOrder );

        // comparison value per value with matlab results
        bool tCheckResidual2 = true; bool tCheckJacobian2 = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual2 = tCheckResidual2 && ( ( tResidual2( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian2 = tCheckJacobian2 && ( ( tJacobian2( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual2 );
        REQUIRE( tCheckJacobian2 );
    }

    //------------------------------------------------------------------------------
    SECTION( "IWG : check IWG computed values - r = delta v * dudt" )
    {
        // residual from matlab
        Matrix< DDRMat > tResidualMatlab = {
        {0.477750000000000},
        {0.992250000000000},
        {2.707250000000000},
        {5.622750000000000}};

        // jacobian from matlab
        Matrix< DDRMat > tJacobianMatlab ={
        {-0.003168750000000,  -0.006581250000000,   0.003168750000000,   0.006581250000000},
        {-0.006581250000000,  -0.013668750000000,   0.006581250000000,   0.013668750000000},
        {-0.017956250000000,  -0.037293750000000,   0.017956250000000,   0.037293750000000},
        {-0.037293750000000,  -0.077456250000000,   0.037293750000000,   0.077456250000000}};

        // check independent evaluation of residual and jacobian

        // set the space and time derivative order for the field u
        uint tSpaceDerivativeOrder = 0;
        uint tTimeDerivativeOrder  = 1;

        // evaluate the residual from IWG
        Matrix< DDRMat > tResidual; tIWGField.compute_residual( tResidual,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // evaluate the jacobian from IWG
        Matrix< DDRMat > tJacobian; tIWGField.compute_jacobian( tJacobian ,
                                                                tSpaceDerivativeOrder,
                                                                tTimeDerivativeOrder );

        // comparison value per value with matlab results
        bool tCheckJacobian = true; bool tCheckResidual = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual = tCheckResidual && ( ( tResidual( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian = tCheckJacobian && ( ( tJacobian( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual );
        REQUIRE( tCheckJacobian );

        // check simultaneous evaluation of residual and jacobian

        // evaluate the residual and the jacobian from IWG
        Matrix< DDRMat > tResidual2, tJacobian2;
        tIWGField.compute_jacobian_and_residual( tJacobian2,
                                                 tResidual2,
                                                 tSpaceDerivativeOrder,
                                                 tTimeDerivativeOrder );

        // comparison value per value with matlab results
        bool tCheckResidual2 = true; bool tCheckJacobian2 = true;
        for ( moris::uint Ik = 0; Ik < tNSpaceTimeBases; Ik++ )
        {
            tCheckResidual2 = tCheckResidual2 && ( ( tResidual2( Ik ) - tResidualMatlab( Ik )) < tEpsilon );
            for ( moris::uint Jk = 0; Jk < tNSpaceTimeBases; Jk++ )
            {
                tCheckJacobian2 = tCheckJacobian2 && ( ( tJacobian2( Ik, Jk ) - tJacobianMatlab( Ik, Jk )) < tEpsilon );
            }
        }
        REQUIRE( tCheckResidual2 );
        REQUIRE( tCheckJacobian2 );
    }
}
