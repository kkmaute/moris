#include <string>
#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"    //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Hamilton-Jacobi", "[moris],[fem],[IWG_HJ]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator(tGeomInterpRule);

    //create space coeff xHat
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

    //create time coeff tHat
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat );

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
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
    tUHat( 4 ) = -1.0;
    tUHat( 5 ) =  7.0;
    tUHat( 6 ) = -6.0;
    tUHat( 7 ) = 82.0;

    //set the coefficients uHat
    tFieldInterpolator->set_coeff( tUHat );

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi( 2, 1 );
    tXi( 0, 0 ) =  0.35; tXi( 1, 0 ) = -0.25;
    Matrix< DDRMat > tTau( 1, 1 );
    tTau( 0, 0 ) = 0.0;

    //set the evaluation point xi, tau
    tFieldInterpolator->set_space_time( tXi, tTau );

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Hamilton_Jacobi_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        // create an IWG Hamilton Jacobi Bulk
        IWG_Hamilton_Jacobi_Bulk tIWGHJBulk( tFieldInterpolator );

        // check evaluation of the residual for IWG_Hamilton_Jacobi_Bulk ?
        //------------------------------------------------------------------------------
        // create a velocity field at evaluation point
        Matrix< DDRMat > tVelocityField( 2, 1, 1.0);

        // evaluate the residual from IWG_Hamilton_Jacobi_Bulk
        Matrix< DDRMat > tResidualHJBulk;
        tIWGHJBulk.compute_residual( tResidualHJBulk,
                                     tVelocityField );

        // check evaluation of the jacobian for IWG_Hamilton_Jacobi_Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Hamilton_Jacobi_Bulk
        Matrix< DDRMat > tJacobianHJBulk;
        tIWGHJBulk.compute_jacobian( tJacobianHJBulk,
                                     tVelocityField );

        //define a boolean for check
        bool tCheckJacobianBulk = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tUHatPert, tResidualHJBulkPert, tJacobianRow;
        real tPert;
        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            //set the perturbed values of uHat
            tUHatPert = tUHat;
            tPert = 1e-6 * tUHatPert( k );
            tUHatPert( k ) = tUHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolator->set_coeff( tUHatPert );

            // compute the perturbed residual
            tIWGHJBulk.compute_residual( tResidualHJBulkPert,
                                         tVelocityField );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualHJBulkPert - tResidualHJBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i=0; i<tNSpaceTimeBases; i++ )
            {
                tCheckJacobianBulk = tCheckJacobianBulk && ( tJacobianRow( i ) - tJacobianHJBulk( i, k ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );

        //set the coefficients uHat
        tFieldInterpolator->set_coeff( tUHat );
    }

    SECTION( "IWG_Hamilton_Jacobi_Bulk : check simultaneous evaluation of residual and jacobian" )
    {
        // define an epsilon environment
        double tEpsilon2 = 1E-12;

        // create an IWG Hamilton Jacobi Bulk
        IWG_Hamilton_Jacobi_Bulk tIWGHJBulk( tFieldInterpolator );

        // create a velocity field at evaluation point
        Matrix< DDRMat > tVelocityField( 2, 1, 1.0);

        // evaluate the residual from IWG_Hamilton_Jacobi_Bulk
        Matrix< DDRMat > tResidualHJBulk;
        tIWGHJBulk.compute_residual( tResidualHJBulk,
                                     tVelocityField );

        // evaluate the jacobian from IWG_Hamilton_Jacobi_Bulk
        Matrix< DDRMat > tJacobianHJBulk;
        tIWGHJBulk.compute_jacobian( tJacobianHJBulk,
                                     tVelocityField );

        // evaluate the residual and jacobian from IWG_Helmholtz_Bulk
        Matrix< DDRMat > tResidualHJBulk2, tJacobianHJBulk2;
        tIWGHJBulk.compute_jacobian_and_residual( tJacobianHJBulk2,
                                                  tResidualHJBulk2,
                                                  tVelocityField );

        // check simulatneous evaluation of the jacobian and residual
        // for IWG Hamilton Jacobi Bulk

        //define a boolean for check
        bool tCheckResidualBulkSim = true;
        bool tCheckJacobianBulkSim = true;

        for( uint k=0; k<tNSpaceTimeBases; k++ )
        {
            // check the value of the residual
            tCheckResidualBulkSim = tCheckResidualBulkSim && ( tResidualHJBulk2( k ) - tResidualHJBulk( k ) < tEpsilon2 );

            // check the value of the jacobian
           for( uint i=0; i<tNSpaceTimeBases; i++ )
           {
               tCheckJacobianBulkSim = tCheckJacobianBulkSim && ( tJacobianHJBulk2( i, k ) - tJacobianHJBulk( i, k ) < tEpsilon2 );
           }
        }
        REQUIRE( tCheckResidualBulkSim );
        REQUIRE( tCheckJacobianBulkSim );
    }

}
