/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Field_IWG_LSNormal.cpp
 *
 */

#include <string>
#include <catch.hpp>

#include "moris_typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src

#include "cl_FEM_Enums.hpp"                       //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"          //FEM//INT//src
#include "cl_FEM_IWG_LSNormal_Bulk.hpp"    //FEM//INT//src

#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_LSNormal", "[moris],[fem],[IWG_LSNormal]" )
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
    Field_Interpolator * tnPhi = new Field_Interpolator( 2,
                                                         tInterpolationRule,
                                                         tGeomInterpolator );
    //get the number of basis for space time
    uint tNSpaceTimeBasesnPhi = tnPhi->get_number_of_space_time_bases();

    //create field coeff vNHat
    Matrix< DDRMat > tnPhiHat( tNSpaceTimeBasesnPhi , 2 );
    tnPhiHat( 0, 0 ) = -1.0; tnPhiHat( 0, 1 ) = -1.0;
    tnPhiHat( 1, 0 ) =  7.0; tnPhiHat( 1, 1 ) = -1.0;
    tnPhiHat( 2, 0 ) = -6.0; tnPhiHat( 2, 1 ) = -1.0;
    tnPhiHat( 3, 0 ) = 82.0; tnPhiHat( 3, 1 ) = -1.0;
    tnPhiHat( 4, 0 ) = -1.0; tnPhiHat( 4, 1 ) = -1.0;
    tnPhiHat( 5, 0 ) =  7.0; tnPhiHat( 5, 1 ) = -1.0;
    tnPhiHat( 6, 0 ) = -6.0; tnPhiHat( 6, 1 ) = -1.0;
    tnPhiHat( 7, 0 ) = 82.0; tnPhiHat( 7, 1 ) = -1.0;

    //set the coefficients uHat
    tnPhi->set_coeff( tnPhiHat );

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
    tnPhi->set_space_time( tParamPoint );
    tPhi->set_space_time( tParamPoint );

    Cell< Field_Interpolator* > tFieldInterpolators( 2 );
    tFieldInterpolators( 0 ) = tnPhi;
    tFieldInterpolators( 1 ) = tPhi;

    // define an epsilon environment
    double tEpsilon = 1E-6;

    SECTION( "IWG_LSNormal_Bulk : check residual and jacobian" )
    {
        // IWG
        //------------------------------------------------------------------------------
        // create an IWG LSnormal Bulk
        IWG_LSNormal_Bulk tIWGLSNormalBulk;

        // set field interpolators
        tIWGLSNormalBulk.set_field_interpolators( tFieldInterpolators );

        // check evaluation of the residual for IWG_LSNormal_Bulk ?
        //------------------------------------------------------------------------------
        // evaluate the residual from IWG_LSNormal_Bulk
        Matrix< DDRMat > tResidualLSNormalBulk;
        tIWGLSNormalBulk.compute_residual( tResidualLSNormalBulk );

        // check evaluation of the jacobian j_nPhi_nPhi for IWG_LSNormal_Bulk by FD
        //------------------------------------------------------------------------------
        // evaluate the jacobian from IWG_Hamilton_Jacobi_Bulk
        Cell< Matrix< DDRMat > > tJacobiansLSNormalBulk( 2 );
        tIWGLSNormalBulk.compute_jacobian( tJacobiansLSNormalBulk );

        // number of nPhi dofs
        uint tNDofsnPhi = tNSpaceTimeBasesnPhi * 2;

        //define a boolean for check
        bool tCheckJacobianBulk = true;

        //evaluate the jacobian by FD
        for( uint k = 0; k < tNDofsnPhi; k++ )
        {
            //set the perturbed values of vNHat
        	Matrix< DDRMat > tnPhiHatPert = tnPhiHat;

            //define a perturbation of the dofs
            real tPert;
            if ( k < 8 )
            {
                tPert = 1e-6 * tnPhiHatPert( k, 0 );
                tnPhiHatPert( k, 0 ) = tnPhiHatPert( k, 0 ) + tPert;
            }
            else
            {
                tPert = 1e-6 * tnPhiHatPert( k-8, 1 );
                tnPhiHatPert( k-8, 1 ) = tnPhiHatPert( k-8, 1 ) + tPert;
            }

            //set the coefficients nPhiHatPert
            tFieldInterpolators( 0 )->set_coeff( tnPhiHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualLSNormalBulkPert;
            tIWGLSNormalBulk.compute_residual( tResidualLSNormalBulkPert );

            // compute the jacobian by FD for the kth nPhiHat
            Matrix< DDRMat > tJacobianRow;
            tJacobianRow = ( tResidualLSNormalBulkPert - tResidualLSNormalBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNDofsnPhi; i++ )
            {
                //std::cout<<tJacobianRow( i )<<std::endl;
                //std::cout<<tJacobiansLSNormalBulk( 0 )( i, k )<<std::endl;
                //std::cout<<"---------"<<std::endl;
                tCheckJacobianBulk = tCheckJacobianBulk
                                   && ( std::abs( tJacobianRow( i ) - tJacobiansLSNormalBulk( 0 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );

        // reset the coefficients nPhiHat
        tFieldInterpolators( 0 )->set_coeff( tnPhiHat );

        // number of phi dofs
        uint tNDofsPhi = tNSpaceTimeBasesPhi * 1;

        // reset the boolean for check
        tCheckJacobianBulk = true;

        //evaluate the jacobian by FD
        for( uint k = 0; k < tNDofsPhi; k++ )
        {
            // set the perturbed values of phiHat
        	Matrix< DDRMat > tPhiHatPert = tPhiHat;
            real tPert = 1e-6 * tPhiHatPert( k );
            tPhiHatPert( k ) = tPhiHatPert( k ) + tPert;

            // set the coefficients phiHatPert
            tFieldInterpolators( 1 )->set_coeff( tPhiHatPert );

            // compute the perturbed residual
            Matrix< DDRMat > tResidualLSNormalBulkPert;
            tIWGLSNormalBulk.compute_residual( tResidualLSNormalBulkPert );

            // compute the jacobian by FD for the kth phiHat
            Matrix< DDRMat > tJacobianRow;
            tJacobianRow = ( tResidualLSNormalBulkPert - tResidualLSNormalBulk ) / tPert;

            // check the value of the jacobian evaluated from the IWG with FD
            for( uint i = 0; i < tNDofsnPhi; i++ )
            {
                //std::cout<<tJacobianRow( i )<<std::endl;
                //std::cout<<tJacobiansLSNormalBulk( 1 )( i, k )<<std::endl;
                //std::cout<<"---------"<<std::endl;
                tCheckJacobianBulk = tCheckJacobianBulk
                                   && ( std::abs( tJacobianRow( i ) - tJacobiansLSNormalBulk( 1 )( i, k ) ) < tEpsilon );
            }
        }
        REQUIRE( tCheckJacobianBulk );

    }
}
