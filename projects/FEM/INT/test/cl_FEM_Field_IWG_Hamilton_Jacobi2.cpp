/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_IWG_Hamilton_Jacobi2.cpp
 *
 */

#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                       //FEM//INT/src
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp"          //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Hamilton_Jacobi2", "[moris],[fem],[IWG_HJ2]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR );

    //create a space time geometry interpolator
    Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

    //create space coeff xHat
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

    //create time coeff tHat
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0; tTHat( 1 ) = 1.0;

    //set the coefficients xHat, tHat
    tGeomInterpolator->set_coeff( tXHat, tTHat );

    // interpolator for velocity field----------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tVNInterRule ( mtk::Geometry_Type::QUAD,
                                      mtk::Interpolation_Type::LAGRANGE,
                                      mtk::Interpolation_Order::LINEAR,
                                      mtk::Interpolation_Type::CONSTANT,
                                      mtk::Interpolation_Order::CONSTANT );
    //create a field interpolator
    uint tNVNFields = 1;
    Field_Interpolator* tVN = new Field_Interpolator( tNVNFields,
                                                      tVNInterRule,
                                                      tGeomInterpolator );

    // get the number of bases for vN
    uint tNVNBases = tVN->get_number_of_space_time_bases();

    // create and set the field coefficients vNHat
    Matrix< DDRMat > tVNHat( tNVNBases , tNVNFields );
    tVNHat( 0 ) = -1.0;
    tVNHat( 1 ) =  7.0;
    tVNHat( 2 ) = -6.0;
    tVNHat( 3 ) = 82.0;
    tVNHat( 4 ) = -1.0;
    tVNHat( 5 ) =  7.0;
    tVNHat( 6 ) = -6.0;
    tVNHat( 7 ) = 82.0;
    tVN->set_coeff( tVNHat );

    // interpolator for level set field----------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tPhiInterRule ( mtk::Geometry_Type::QUAD,
                                       mtk::Interpolation_Type::LAGRANGE,
                                       mtk::Interpolation_Order::LINEAR,
                                       mtk::Interpolation_Type::LAGRANGE,
                                       mtk::Interpolation_Order::LINEAR );

    //create a field interpolator
    uint tNPhiFields = 1;
    Field_Interpolator* tPhi = new Field_Interpolator( tNPhiFields,
                                                       tPhiInterRule,
                                                       tGeomInterpolator );

    // get the number of bases for phi
    uint tNPhiBases = tPhi->get_number_of_space_time_bases();

    // create and set the field coefficients vNHat
    Matrix< DDRMat > tPhiHat( tNVNBases , tNVNFields );
    tPhiHat( 0 ) = -1.0;
    tPhiHat( 1 ) =  7.0;
    tPhiHat( 2 ) = -6.0;
    tPhiHat( 3 ) = 82.0;
    tPhiHat( 4 ) = -1.0;
    tPhiHat( 5 ) =  7.0;
    tPhiHat( 6 ) = -6.0;
    tPhiHat( 7 ) = 82.0;
    tPhi->set_coeff( tPhiHat );

    // create and set the evaluation point xi, tau
    Matrix< DDRMat > tParamPoint = { {  0.35 },
                                     { -0.25 },
                                     {  0.0 } };
    tVN ->set_space_time( tParamPoint );
    tPhi->set_space_time( tParamPoint );

	// create a cell of field interpolators with tVN and tPhi
    Vector< Field_Interpolator * > tFieldInterpolators( 2 , nullptr);
    tFieldInterpolators( 1 ) = tVN;
    tFieldInterpolators( 0 ) = tPhi;

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Hamilton_Jacobi_Bulk : check residual and jacobian" )
    {
        // create an IWG Helmholtz Bulk2
        IWG_Hamilton_Jacobi_Bulk2 tIWGHJBulk;

        // set field interpolators
        tIWGHJBulk.set_field_interpolators( tFieldInterpolators );

        // check evaluation of the jacobian by FD --------------------------------------
        // evaluate the residual from IWG
        Matrix< DDRMat > tResidualHJBulk;
        tIWGHJBulk.compute_residual( tResidualHJBulk );
        //print( tResidualHJBulk, "r_phi" );

        // evaluate the jacobian from IWG
        Vector< Matrix< DDRMat > > tJacobiansHJBulk( 2 );
        tIWGHJBulk.compute_jacobian( tJacobiansHJBulk );
        //print( tJacobiansHJBulk( 0 ), "j_phi_phi" );
        //print( tJacobiansHJBulk( 1 ), "j_phi_vN" );

        //define a boolean for check
        bool tCheckJacobian = true;

        //evaluate the jacobian by FD
        for( uint k = 0; k < tNPhiBases; k++ )
        {
            //set the perturbed values of phiHat
            Matrix< DDRMat > tPhiHatPert = tPhiHat;
            real tPert = 1e-6 * tPhiHatPert( k );
            tPhiHatPert( k ) = tPhiHatPert( k ) + tPert;

            //set the coefficients phiHatPert
            tFieldInterpolators( 0 )->set_coeff( tPhiHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualHJBulkPert;
            tIWGHJBulk.compute_residual( tResidualHJBulkPert );

            // compute the jacobian by FD for the kth phiHat
            Matrix< DDRMat > tJacobianRow = ( tResidualHJBulkPert - tResidualHJBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNPhiBases; i++ )
            {
                //std::cout<<tJacobianRow( i )<<std::endl;
                //std::cout<<tJacobiansHJBulk( 0 )( i, k )<<std::endl;
                //std::cout<<"-----------"<<std::endl;
                tCheckJacobian = tCheckJacobian && ( std::abs( tJacobianRow( i ) - tJacobiansHJBulk( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobian );

        // reset the coefficients vNHat
        tFieldInterpolators( 0 )->set_coeff( tPhiHat );

        //reset boolean for check
        tCheckJacobian = true;

        //evaluate the jacobian by FD
        for( uint k = 0; k < tNVNBases; k++ )
        {
            //set the perturbed values of vNHat
            Matrix< DDRMat > tVNHatPert = tVNHat;
            real tPert = 1e-6 * tVNHatPert( k );
            tVNHatPert( k ) = tVNHatPert( k ) + tPert;

            //set the coefficients vNHatPert
            tFieldInterpolators( 1 )->set_coeff( tVNHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualHJBulkPert;
            tIWGHJBulk.compute_residual( tResidualHJBulkPert );

            // compute the jacobian by FD for the kth uHat
            Matrix< DDRMat > tJacobianRow = ( tResidualHJBulkPert - tResidualHJBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNVNBases; i++ )
            {
                //std::cout<<tJacobianRow( i )<<std::endl;
                //std::cout<<tJacobiansHJBulk( 1 )( i, k )<<std::endl;
                //std::cout<<"-----------"<<std::endl;
                tCheckJacobian = tCheckJacobian && ( std::abs( tJacobianRow( i ) - tJacobiansHJBulk( 1 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobian );

        //clean up
        delete tGeomInterpolator;

        uint tNumOfFieldInterpolators = tFieldInterpolators.size();
        for(uint i = 0; i < tNumOfFieldInterpolators; i++ )
        {
            delete tFieldInterpolators( i );
        }

//
//        //evaluate the jacobian by FD
//        for( uint k = 0; k < tNPhiBases; k++ )
//        {
//            //set the perturbed values of vNHat
//            Matrix< DDRMat > tPhiHatPert = tPhiHat;
//            real tPert = 1e-6 * tPhiHatPert( k );
//            tPhiHatPert( k ) = tPhiHatPert( k ) + tPert;
//
//            //set the coefficients vNHatPert
//            tFieldInterpolators( 1 )->set_coeff( tPhiHatPert );
//
//            // compute the perturbed residual
//            Matrix< DDRMat > tResidualHelmBulkPert;
//            tIWGHelmBulk.compute_residual( tResidualHelmBulkPert,
//                                           tFieldInterpolators );
//
//            // compute the jacobian by FD for the kth uHat
//            Matrix< DDRMat > tJacobianRow = ( tResidualHelmBulkPert - tResidualHelmBulk ) / tPert;
//
//            // check the value of the jacobian evaluated from the IWG with FD
//            for( uint i = 0; i < tNPhiBases; i++ )
//            {
//                //std::cout<<tJacobianRow( i )<<std::endl;
//                //std::cout<<tJacobiansHelmBulk( 1 )( i, k )<<std::endl;
//                //std::cout<<"-----------"<<std::endl;
//                tCheckJacobian = tCheckJacobian && ( std::abs( tJacobianRow( i ) - tJacobiansHelmBulk( 1 )( i, k ) ) < tEpsilon );
//            }
//        }
//        REQUIRE( tCheckJacobian );
    }
}
