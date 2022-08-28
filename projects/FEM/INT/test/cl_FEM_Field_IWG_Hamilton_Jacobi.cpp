/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_IWG_Hamilton_Jacobi.cpp
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
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"    //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Hamilton-Jacobi", "[moris],[fem],[IWG_HJ]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR,
                                        mtk::Interpolation_Type::LAGRANGE,
                                        mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator(tGeomInterpRule);

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
    tGeomInterpolator->set_coeff( tXHat, tTHat );

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    mtk::Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                            mtk::Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR,
                                            mtk::Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    //create a field interpolator
    Field_Interpolator * tvN = new Field_Interpolator( 2,
                                                       tInterpolationRule,
                                                        tGeomInterpolator );
    //get the number of basis for space time
    //uint tNSpaceTimeBases = tFieldInterpolator->get_number_of_space_time_bases();
    uint tNSpaceTimeBasesvN = tvN->get_number_of_space_time_bases();

    //create field coeff vNHat
    Matrix< DDRMat > tvNHat( tNSpaceTimeBasesvN , 2 );
    tvNHat( 0, 0 ) = -1.0; tvNHat( 0, 1 ) = -1.0;
    tvNHat( 1, 0 ) =  7.0; tvNHat( 1, 1 ) = -1.0;
    tvNHat( 2, 0 ) = -6.0; tvNHat( 2, 1 ) = -1.0;
    tvNHat( 3, 0 ) = 82.0; tvNHat( 3, 1 ) = -1.0;
    tvNHat( 4, 0 ) = -1.0; tvNHat( 4, 1 ) = -1.0;
    tvNHat( 5, 0 ) =  7.0; tvNHat( 5, 1 ) = -1.0;
    tvNHat( 6, 0 ) = -6.0; tvNHat( 6, 1 ) = -1.0;
    tvNHat( 7, 0 ) = 82.0; tvNHat( 7, 1 ) = -1.0;

    //set the coefficients uHat
    tvN->set_coeff( tvNHat );

    //create a field interpolator
    Field_Interpolator * tPhi = new Field_Interpolator( 1,
                                                        tInterpolationRule,
                                                        tGeomInterpolator );
    //get the number of basis for space time
    uint tNSpaceTimeBasesPhi = tPhi->get_number_of_space_time_bases();

    //create field coeff uHat
    Matrix< DDRMat > tPhiHat( tNSpaceTimeBasesPhi , 1 );
    tPhiHat( 0 ) = -1.0;
    tPhiHat( 1 ) =  7.0;
    tPhiHat( 2 ) = -6.0;
    tPhiHat( 3 ) = 82.0;
    tPhiHat( 4 ) = -1.0;
    tPhiHat( 5 ) =  7.0;
    tPhiHat( 6 ) = -6.0;
    tPhiHat( 7 ) = 82.0;

    //set the coefficients uHat
    tPhi->set_coeff( tPhiHat );

    // create evaluation point xi, tau
    Matrix< DDRMat > tXi( 2, 1 );
    tXi( 0, 0 ) =  0.35; tXi( 1, 0 ) = -0.25;
    Matrix< DDRMat > tTau( 1, 1 );
    tTau( 0, 0 ) = 0.0;
    Matrix< DDRMat > tParamPoint = { { tXi( 0 ) },
                                     { tXi( 1 ) },
                                     { tTau( 0 ) } };

    //set the evaluation point xi, tau
    tvN->set_space_time( tParamPoint );
    tPhi->set_space_time( tParamPoint );

    Cell< Field_Interpolator* > tFieldInterpolators( 2 );
    tFieldInterpolators( 0 ) = tPhi;
    tFieldInterpolators( 1 ) = tvN;

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Hamilton_Jacobi_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        // create an IWG Hamilton Jacobi Bulk
        IWG_Hamilton_Jacobi_Bulk tIWGHJBulk;

        // set field interpolators
        tIWGHJBulk.set_field_interpolators( tFieldInterpolators );

        // check evaluation of the residual for IWG_Hamilton_Jacobi_Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual from IWG_Hamilton_Jacobi_Bulk
        Matrix< DDRMat > tResidualHJBulk;
        tIWGHJBulk.compute_residual( tResidualHJBulk,
                                     tFieldInterpolators );

        // check evaluation of the jacobian for IWG_Hamilton_Jacobi_Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Hamilton_Jacobi_Bulk
        Cell< Matrix< DDRMat > > tJacobiansHJBulk( 2 );
        tIWGHJBulk.compute_jacobian( tJacobiansHJBulk,
                                     tFieldInterpolators );

        //define a boolean for check
        bool tCheckJacobianBulk = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tPhiHatPert, tResidualHJBulkPert, tJacobianRow;
        real tPert;
        for( uint k=0; k<tNSpaceTimeBasesPhi; k++ )
        {
            //set the perturbed values of uHat
            tPhiHatPert = tPhiHat;
            tPert = 1e-6 * tPhiHatPert( k );
            tPhiHatPert( k ) = tPhiHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolators( 0 )->set_coeff( tPhiHatPert );

            // compute the perturbed residual
            tIWGHJBulk.compute_residual( tResidualHJBulkPert );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualHJBulkPert - tResidualHJBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i=0; i<tNSpaceTimeBasesPhi; i++ )
            {
                tCheckJacobianBulk = tCheckJacobianBulk && ( std::abs( tJacobianRow( i ) - tJacobiansHJBulk( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );

        //set the coefficients uHat
        tFieldInterpolators( 0 )->set_coeff( tPhiHat );

        //evaluate the jacobian by FD
        tCheckJacobianBulk = true;
        Matrix< DDRMat > tvNHatPert;
        for( uint k = 0; k < tNSpaceTimeBasesvN; k++ )
        {
            //set the perturbed values of vNHat
            tvNHatPert = tvNHat;
            if ( k < 8 )
            {
                tPert = 1e-6 * tvNHatPert( k, 0 );
                tvNHatPert( k, 0 ) = tvNHatPert( k, 0 ) + tPert;
            }
            else
            {
                tPert = 1e-6 * tvNHatPert( k, 1 );
                tvNHatPert( k, 1 ) = tvNHatPert( k, 1 ) + tPert;
            }

            //set the coefficients uHatPert
            tFieldInterpolators( 1 )->set_coeff( tvNHatPert );

            // compute the perturbed residual
            tIWGHJBulk.compute_residual( tResidualHJBulkPert );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualHJBulkPert - tResidualHJBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNSpaceTimeBasesvN; i++ )
            {
                tCheckJacobianBulk = tCheckJacobianBulk && ( std::abs( tJacobianRow( i ) - tJacobiansHJBulk( 1 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );

    }

//    SECTION( "IWG_Hamilton_Jacobi_Bulk : check simultaneous evaluation of residual and jacobian" )
//    {
//        // define an epsilon environment
//        double tEpsilon2 = 1E-12;
//
//        // create an IWG Hamilton Jacobi Bulk
//        IWG_Hamilton_Jacobi_Bulk tIWGHJBulk( tFieldInterpolator );
//
//        // create a velocity field at evaluation point
//        Matrix< DDRMat > tVelocityField( 2, 1, 1.0);
//
//        // evaluate the residual from IWG_Hamilton_Jacobi_Bulk
//        Matrix< DDRMat > tResidualHJBulk;
//        tIWGHJBulk.compute_residual( tResidualHJBulk,
//                                     tVelocityField );
//
//        // evaluate the jacobian from IWG_Hamilton_Jacobi_Bulk
//        Matrix< DDRMat > tJacobianHJBulk;
//        tIWGHJBulk.compute_jacobian( tJacobianHJBulk,
//                                     tVelocityField );
//
//        // evaluate the residual and jacobian from IWG_Helmholtz_Bulk
//        Matrix< DDRMat > tResidualHJBulk2, tJacobianHJBulk2;
//        tIWGHJBulk.compute_jacobian_and_residual( tJacobianHJBulk2,
//                                                  tResidualHJBulk2,
//                                                  tVelocityField );
//
//        // check simulatneous evaluation of the jacobian and residual
//        // for IWG Hamilton Jacobi Bulk
//
//        //define a boolean for check
//        bool tCheckResidualBulkSim = true;
//        bool tCheckJacobianBulkSim = true;
//
//        for( uint k=0; k<tNSpaceTimeBases; k++ )
//        {
//            // check the value of the residual
//            tCheckResidualBulkSim = tCheckResidualBulkSim && ( tResidualHJBulk2( k ) - tResidualHJBulk( k ) < tEpsilon2 );
//
//            // check the value of the jacobian
//           for( uint i=0; i<tNSpaceTimeBases; i++ )
//           {
//               tCheckJacobianBulkSim = tCheckJacobianBulkSim && ( tJacobianHJBulk2( i, k ) - tJacobianHJBulk( i, k ) < tEpsilon2 );
//           }
//        }
//        REQUIRE( tCheckResidualBulkSim );
//        REQUIRE( tCheckJacobianBulkSim );
//    }

}
