/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_IWG_Helmholtz.cpp
 *
 */

#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                       //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"          //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Helmholtz", "[moris],[fem],[IWG_Helmholtz]" )
{

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
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

    // field interpolator
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    mtk::Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                            mtk::Interpolation_Type::LAGRANGE,
                                            mtk::Interpolation_Order::LINEAR );

    //set the number of interpolated fields
    uint tNumberOfFields = 1;

    //create a field interpolator
    Field_Interpolator* tvN = new Field_Interpolator( tNumberOfFields,
                                                      tInterpolationRule,
                                                      tGeomInterpolator);

    //get the number of basis for space time
    uint tNSpaceTimeBases = tvN->get_number_of_space_time_bases();

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
    tvN->set_coeff( tUHat );

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

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_Helmholtz_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        Cell< Field_Interpolator * > tFieldInterpolators( 1 , nullptr);
        tFieldInterpolators( 0 ) = tvN;

        // create an IWG Helmholtz Bulk
        IWG_Helmholtz_Bulk tIWGHelmBulk;

        // set field interpolators
        tIWGHelmBulk.set_field_interpolators( tFieldInterpolators );

        // check evaluation of the residual for IWG Helmholtz Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual from IWG_Helmholtz_Bulk
        Matrix< DDRMat > tResidualHelmBulk;
        tIWGHelmBulk.compute_residual( tResidualHelmBulk );

        // check evaluation of the jacobian for IWG Helmholtz Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Helmholtz_Bulk
        Cell< Matrix< DDRMat > > tJacobiansHelmBulk( 1 );
        tIWGHelmBulk.compute_jacobian( tJacobiansHelmBulk );

        //define a boolean for check
        bool tCheckJacobianBulk = true;

        //evaluate the jacobian by FD
        Matrix< DDRMat > tUHatPert, tResidualHelmBulkPert, tJacobianRow;
        real tPert;
        for( uint k = 0; k < tNSpaceTimeBases; k++ )
        {
            //set the perturbed values of uHat
            tUHatPert = tUHat;
            tPert = 1e-4 * tUHatPert( k );
            tUHatPert( k ) = tUHatPert( k ) + tPert;

            //set the coefficients uHatPert
            tFieldInterpolators( 0 )->set_coeff( tUHatPert );

            // compute the perturbed residual
            tIWGHelmBulk.compute_residual( tResidualHelmBulkPert );

            // compute the jacobian by FD for the kth uHat
            tJacobianRow = ( tResidualHelmBulkPert - tResidualHelmBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNSpaceTimeBases; i++ )
            {
                tCheckJacobianBulk = tCheckJacobianBulk && ( std::abs( tJacobianRow( i ) - tJacobiansHelmBulk( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );
    }

//    SECTION( "IWG_Helmholtz_Interface : check residual and jacobian" )
//    {
//        // check evaluation of the residual for IWG Helmholtz Interface ?
//        //------------------------------------------------------------------------------
//        // create a Helmholtz filter length parameter
//        double tkappa = 2.0;
//
//        //create an IWG Helmholtz Interface
//        IWG_Helmholtz_Interface tIWGHelmInterface( tkappa );
//
//        // create a normal to the interface
//        uint tNormalDim = tGeomInterpolator.get_number_of_space_dimensions();
//        Matrix< DDRMat > tNormalInterface( tNormalDim, 1, 1.0 );
//
//        // evaluate the residual from IWG_Helmholtz_Interface
//        Matrix< DDRMat > tResidualHelmInterface;
//        tIWGHelmInterface.compute_residual( tResidualHelmInterface,
//                                            tFieldInterpolator,
//                                            tNormalInterface );
//
//        // check evaluation of the jacobian for IWG Helmholtz Interface by FD
//        //------------------------------------------------------------------------------
//        // evaluate the jacobian from IWG Helmholtz Interface
//        Matrix< DDRMat > tJacobianHelmInterface;
//        tIWGHelmInterface.compute_jacobian( tJacobianHelmInterface,
//                                            tFieldInterpolator,
//                                            tNormalInterface);
//
//        //define a boolean for check
//        bool tCheckJacobianInterface = true;
//
//        //evaluate the jacobian by FD
//        Matrix< DDRMat > tUHatPert, tResidualHelmInterfacePert, tJacobianRow;
//        real tPert;
//        for( uint k=0; k<tNSpaceTimeBases; k++ )
//        {
//            //set the perturbed values of uHat
//            tUHatPert = tUHat;
//            tPert = 1e-6 * tUHatPert( k );
//            tUHatPert( k ) = tUHatPert( k ) + tPert;
//
//            //set the coefficients uHatPert
//            //tFieldInterpolator->set_coeff( tUHatPert );
//            tFieldInterpolator.set_coeff( tUHatPert );
//
//            // compute the perturbed residual
//            tIWGHelmInterface.compute_residual( tResidualHelmInterfacePert,
//                                                tFieldInterpolator,
//                                                tNormalInterface );
//
//            // compute the jacobian by FD for the kth uHat
//            tJacobianRow = ( tResidualHelmInterfacePert - tResidualHelmInterface ) / tPert;
//
//            // check the value of the jacobian evaluated from the IWG with FD
//            for( uint i=0; i<tNSpaceTimeBases; i++ )
//            {
//                tCheckJacobianInterface = tCheckJacobianInterface && ( tJacobianRow( i ) - tJacobianHelmInterface( i, k ) < tEpsilon );
//            }
//        }
//        REQUIRE( tCheckJacobianInterface );
//
//        //set the coefficients uHat
//        //tFieldInterpolator->set_coeff( tUHat );
//        tFieldInterpolator.set_coeff( tUHat );
//    }
//
//    SECTION( "IWG_Helmholtz_Bulk : check simultaneous evaluation of residual and jacobian" )
//    {
//
//        // define an epsilon environment
//        double tEpsilon2 = 1E-12;
//
//        // IWG
//        //------------------------------------------------------------------------------
//        // create a Helmholtz filter length parameter
//        double tkappa = 2.0;
//
//        // create an IWG Helmholtz Bulk
//        IWG_Helmholtz_Bulk tIWGHelmBulk( tkappa );
//
//        // create an unfiltered value at evaluation point
//        double tUnfilteredVHat = 2.0;
//
//        // evaluate the residual from IWG_Helmholtz_Bulk
//        Matrix< DDRMat > tResidualHelmBulk;
//        tIWGHelmBulk.compute_residual( tResidualHelmBulk,
//                                       tFieldInterpolator,
//                                       tUnfilteredVHat );
//
//        // evaluate the jacobian from IWG_Helmholtz_Bulk
//        Matrix< DDRMat > tJacobianHelmBulk;
//        tIWGHelmBulk.compute_jacobian( tJacobianHelmBulk,
//                                       tFieldInterpolator );
//
//        // evaluate the residual and jacobian from IWG_Helmholtz_Bulk
//        Matrix< DDRMat > tResidualHelmBulk2, tJacobianHelmBulk2;
//        tIWGHelmBulk.compute_jacobian_and_residual( tJacobianHelmBulk2,
//                                                    tResidualHelmBulk2,
//                                                    tFieldInterpolator,
//                                                    tUnfilteredVHat );
//
//        // check simulatneous evaluation of the jacobian and residual
//        // for IWG Helmholtz Bulk
//
//        //define a boolean for check
//        bool tCheckResidualBulkSim = true;
//        bool tCheckJacobianBulkSim = true;
//
//        for( uint k=0; k<tNSpaceTimeBases; k++ )
//        {
//            // check the value of the residual
//            tCheckResidualBulkSim = tCheckResidualBulkSim && ( tResidualHelmBulk2( k ) - tResidualHelmBulk( k ) < tEpsilon2 );
//
//            // check the value of the jacobian
//           for( uint i=0; i<tNSpaceTimeBases; i++ )
//           {
//               tCheckJacobianBulkSim = tCheckJacobianBulkSim && ( tJacobianHelmBulk2( i, k ) - tJacobianHelmBulk( i, k ) < tEpsilon2 );
//           }
//        }
//        REQUIRE( tCheckResidualBulkSim );
//        REQUIRE( tCheckJacobianBulkSim );
//    }
//
//    SECTION( "IWG_Helmholtz_Interface : check simultaneous evaluation of residual and jacobian" )
//    {
//
//        // define an epsilon environment
//        double tEpsilon2 = 1E-12;
//
//        // IWG
//        //------------------------------------------------------------------------------
//        // create a Helmholtz filter length parameter
//        double tkappa = 2.0;
//
//        //create an IWG Helmholtz Interface
//        IWG_Helmholtz_Interface tIWGHelmInterface( tkappa );
//
//        // create a normal to the interface
//        uint tNormalDim = tGeomInterpolator.get_number_of_space_dimensions();
//        Matrix< DDRMat > tNormalInterface( tNormalDim, 1, 1.0 );
//
//        // evaluate the residual from IWG_Helmholtz_Interface
//        Matrix< DDRMat > tResidualHelmInterface;
//        tIWGHelmInterface.compute_residual( tResidualHelmInterface,
//                                            tFieldInterpolator,
//                                            tNormalInterface );
//
//        // evaluate the jacobian from IWG Helmholtz Interface
//        Matrix< DDRMat > tJacobianHelmInterface;
//        tIWGHelmInterface.compute_jacobian( tJacobianHelmInterface,
//                                            tFieldInterpolator,
//                                            tNormalInterface);
//
//        // evaluate the residual and jacobian from IWG_Helmholtz_Interface
//        Matrix< DDRMat > tResidualHelmInterface2, tJacobianHelmInterface2;
//        tIWGHelmInterface.compute_jacobian_and_residual( tJacobianHelmInterface2,
//                                                         tResidualHelmInterface2,
//                                                         tFieldInterpolator,
//                                                         tNormalInterface );
//
//        // check simulatneous evaluation of the jacobian and residual
//        // for IWG Helmholtz Bulk
//
//        //define a boolean for check
//        bool tCheckResidualInterfaceSim = true;
//        bool tCheckJacobianInterfaceSim = true;
//
//        for( uint k=0; k<tNSpaceTimeBases; k++ )
//        {
//            // check the value of the residual
//            tCheckResidualInterfaceSim = tCheckResidualInterfaceSim && ( tResidualHelmInterface2( k ) - tResidualHelmInterface( k ) < tEpsilon2 );
//
//            // check the value of the jacobian
//           for( uint i=0; i<tNSpaceTimeBases; i++ )
//           {
//               tCheckJacobianInterfaceSim = tCheckJacobianInterfaceSim && ( tJacobianHelmInterface2( i, k ) - tJacobianHelmInterface( i, k ) < tEpsilon2 );
//           }
//        }
//        REQUIRE( tCheckResidualInterfaceSim );
//        REQUIRE( tCheckJacobianInterfaceSim );
//    }
}
