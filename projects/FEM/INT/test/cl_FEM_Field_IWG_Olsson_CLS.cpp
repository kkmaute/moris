/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_IWG_Olsson_CLS.cpp
 *
 */

#include <string>
#include <catch.hpp>

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"         //FEM//INT//src
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"    //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Olsson_CLS", "[moris],[fem],[IWG_OCLS]" )
{
    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR);

    // create a space and a time geometry interpolator
    Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;

    // create time coeff tHat
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    // set the coefficients xHat, tHat
    tGeomInterpolator->set_coeff( tXHat, tTHat);

    // field interpolators
    //------------------------------------------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                            mtk::Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            mtk::Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    // create a scalar level set field interpolator
    Field_Interpolator * tPhi = new Field_Interpolator( 1,
                                                        tInterpolationRule,
                                                        tGeomInterpolator );

    // get the number of basis for space time
    uint tNBasesPhi = tPhi->get_number_of_space_time_bases();

    // create field coeff PhiHat
    Matrix< DDRMat > tPhiHat( tNBasesPhi , 1 );
    tPhiHat( 0 ) = -1.0;
    tPhiHat( 1 ) =  7.0;
    tPhiHat( 2 ) = -6.0;
    tPhiHat( 3 ) = 82.0;
    tPhiHat( 4 ) = -1.0;
    tPhiHat( 5 ) =  7.0;
    tPhiHat( 6 ) = -6.0;
    tPhiHat( 7 ) = 82.0;

    // set the coefficients phiHat
    tPhi->set_coeff( tPhiHat );

    // create a level set normal field interpolator
    Field_Interpolator * tNPhi = new Field_Interpolator( 2,
                                                         tInterpolationRule,
                                                         tGeomInterpolator );

    // get the number of basis for space time
    uint tNBasesNPhi = tNPhi->get_number_of_space_time_bases();

    // create field coeff nPhiHat
    Matrix< DDRMat > tNPhiHat( tNBasesNPhi , 2 );
    tNPhiHat( 0, 0 ) = -1.0; tNPhiHat( 0, 1 ) = -1.0;
    tNPhiHat( 1, 0 ) =  7.0; tNPhiHat( 1, 1 ) =  7.0;
    tNPhiHat( 2, 0 ) = -6.0; tNPhiHat( 2, 1 ) = -6.0;
    tNPhiHat( 3, 0 ) = 82.0; tNPhiHat( 3, 1 ) =  82.0;
    tNPhiHat( 4, 0 ) = -1.0; tNPhiHat( 4, 1 ) = -1.0;
    tNPhiHat( 5, 0 ) =  7.0; tNPhiHat( 5, 1 ) =  7.0;
    tNPhiHat( 6, 0 ) = -6.0; tNPhiHat( 6, 1 ) = -6.0;
    tNPhiHat( 7, 0 ) = 82.0; tNPhiHat( 7, 1 ) =  82.0;

    // set the coefficients nPhiHat
    tNPhi->set_coeff( tNPhiHat );

    // create a cell of field interpolators
    Cell< Field_Interpolator* > tFieldInterpolators( 2 );
    tFieldInterpolators( 0 ) = tPhi;
    tFieldInterpolators( 1 ) = tNPhi;

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi( 2, 1 );
    tXi( 0, 0 ) =  0.35; tXi( 1, 0 ) = -0.25;
    Matrix< DDRMat > tTau( 1, 1 );
    tTau( 0, 0 ) = 0.0;
    Matrix< DDRMat > tParamPoint = { { tXi( 0 ) },
                                     { tXi( 1 ) },
                                     { tTau( 0 ) } };

    //set the evaluation point xi, tau
    tFieldInterpolators( 0 )->set_space_time( tParamPoint );
    tFieldInterpolators( 1 )->set_space_time( tParamPoint );

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Olsson_CLS_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        // create an IWG Olsson CLS Bulk
        IWG_Olsson_CLS_Bulk tIWGOCLSBulk;

        // set field interpolators
        tIWGOCLSBulk.set_field_interpolators( tFieldInterpolators );

        // check evaluation of the residual for IWG Olsson CLS Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual from IWG_Olsson_CLS_Bulk
        Matrix< DDRMat > tResidualOCLSBulk;
        tIWGOCLSBulk.compute_residual( tResidualOCLSBulk );

        // check evaluation of the jacobian for IWG Olsson CLS Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Olsson_CLS_Bulk
        Cell< Matrix< DDRMat > > tJacobiansOCLSBulk( 2 );
        tIWGOCLSBulk.compute_jacobian( tJacobiansOCLSBulk );

        //define a boolean for check
        bool tJacobianBulk = true;

        // number of phi dofs
        uint tNDofsPhi = tNBasesPhi * 1;

        //evaluate the jacobian by FD
        for( uint k=0; k < tNDofsPhi; k++ )
        {
            //set the perturbed values of phiHat
            Matrix< DDRMat > tPhiHatPert = tPhiHat;
            real tPert = 1e-6 * tPhiHatPert( k );
            tPhiHatPert( k ) = tPhiHatPert( k ) + tPert;

            //set the coefficients phiHatPert
            tFieldInterpolators( 0 )->set_coeff( tPhiHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualOCLSBulkPert;
            tIWGOCLSBulk.compute_residual( tResidualOCLSBulkPert );

            // compute the jacobian by FD for the kth phiHat
            Matrix< DDRMat > tJacobianRow;
            tJacobianRow = ( tResidualOCLSBulkPert - tResidualOCLSBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNDofsPhi; i++ )
            {
                tJacobianBulk = tJacobianBulk && ( std::abs( tJacobianRow( i ) - tJacobiansOCLSBulk( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tJacobianBulk );
        // reset the coefficients phiHat
        tFieldInterpolators( 0 )->set_coeff( tPhiHat );

        // reset boolean for check
        tJacobianBulk = true;

        // number of nPhi dofs
        uint tNDofsNPhi = tNBasesNPhi * 2;

        //evaluate the jacobian by FD
        for( uint k=0; k < tNDofsNPhi; k++ )
        {
            //set the perturbed values of phiHat
        	Matrix< DDRMat > tNPhiHatPert = tNPhiHat;
            //define a perturbation of the dofs
            real tPert;
            if ( k < 8 )
            {
                tPert = 1e-6 * tNPhiHatPert( k, 0 );
                tNPhiHatPert( k, 0 ) = tNPhiHatPert( k, 0 ) + tPert;
            }
            else
            {
                tPert = 1e-6 * tNPhiHatPert( k-8, 1 );
                tNPhiHatPert( k-8, 1 ) = tNPhiHatPert( k-8, 1 ) + tPert;
            }

            //set the coefficients phiHatPert
            tFieldInterpolators( 1 )->set_coeff( tNPhiHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualOCLSBulkPert;
            tIWGOCLSBulk.compute_residual( tResidualOCLSBulkPert );

            // compute the jacobian by FD for the kth phiHat
            Matrix< DDRMat > tJacobianRow;
            tJacobianRow = ( tResidualOCLSBulkPert - tResidualOCLSBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNDofsPhi; i++ )
            {
//            	std::cout<<tJacobianRow( i )<<std::endl;
//            	std::cout<<tJacobiansOCLSBulk( 1 )( i, k )<<std::endl;
//            	std::cout<<"----------"<<std::endl;
                tJacobianBulk = tJacobianBulk && ( std::abs( tJacobianRow( i ) - tJacobiansOCLSBulk( 1 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tJacobianBulk );
    }

//    SECTION( "IWG_Olsson_CLS_Bulk : check simultaneous evaluation of residual and jacobian" )
//    {
//
//        // define an epsilon environment
//        double tEpsilon2 = 1E-12;
//
//        // set field upper and lower bounds
//        double tFieldUpperbound =  1.0;
//        double tFieldLowerbound = -1.0;
//
//        // set a method parameter epsilon
//        double tEpsilonParameter = 1.0;
//
//        // create an IWG Olsson CLS Bulk
//        IWG_Olsson_CLS_Bulk tIWGOCLSBulk( tFieldInterpolator,
//                                          tFieldUpperbound,
//                                          tFieldLowerbound,
//                                          tEpsilonParameter );
//
//        // check evaluation of the residual for IWG Olsson CLS Bulk ?
//        //------------------------------------------------------------------------------
//        // create a velocity field at evaluation point
//        Matrix< DDRMat > tFieldNormal( 2, 1, 1.0);
//
//        // evaluate the residual from IWG_Olsson_CLS_Bulk
//        Matrix< DDRMat > tResidualOCLSBulk;
//        tIWGOCLSBulk.compute_residual( tResidualOCLSBulk,
//                                       tFieldNormal );
//
//        // check evaluation of the jacobian for IWG Olsson CLS Bulk by FD
//        //------------------------------------------------------------------------------
//        // evaluate the jacobian from IWG_Olsson_CLS_Bulk
//        Matrix< DDRMat > tJacobianOCLSBulk;
//        tIWGOCLSBulk.compute_jacobian( tJacobianOCLSBulk,
//                                       tFieldNormal );
//
//        // evaluate the residual and jacobian from IWG_Olsson_CLS_Bulk
//        Matrix< DDRMat > tResidualOCLSBulk2, tJacobianOCLSBulk2;
//        tIWGOCLSBulk.compute_jacobian_and_residual( tJacobianOCLSBulk2,
//                                                    tResidualOCLSBulk2,
//                                                    tFieldNormal );
//
//        // check simulatneous evaluation of the jacobian and residual
//        // for IWG Olsson CLS Bulk
//
//        //define a boolean for check
//        bool tCheckResidualBulkSim = true;
//        bool tCheckJacobianBulkSim = true;
//
//        for( uint k=0; k<tNSpaceTimeBases; k++ )
//        {
//            // check the value of the residual
//            tCheckResidualBulkSim = tCheckResidualBulkSim && ( tResidualOCLSBulk2( k ) - tResidualOCLSBulk( k ) < tEpsilon2 );
//
//            // check the value of the jacobian
//           for( uint i=0; i<tNSpaceTimeBases; i++ )
//           {
//               tCheckJacobianBulkSim = tCheckJacobianBulkSim && ( tJacobianOCLSBulk2( i, k ) - tJacobianOCLSBulk( i, k ) < tEpsilon2 );
//           }
//        }
//        REQUIRE( tCheckResidualBulkSim );
//        REQUIRE( tCheckJacobianBulkSim );
//    }

}
