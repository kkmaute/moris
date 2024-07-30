/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_FS_Struc_Interface.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private public
// FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private

// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
// LINALG
#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

// This UT tests the pressure interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order
TEST_CASE( "IWG_FS_Struc_Interface", "[moris],[fem],[IWG_FS_Struc_Interface]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
        // create geometry type
        mtk::Geometry_Type tGeometryType;

        // create space coeff xHat
        Matrix< DDRMat > tXHat;

        // create evaluation point xi, tau
        Matrix< DDRMat > tParamPoint;

        // create list with number of coeffs
        Matrix< DDRMat > tNumCoeffs;

        // create the normal
        Matrix< DDRMat > tNormal;

        // dof type list
        Vector< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

        Vector< Vector< MSI::Dof_Type > > tVelDofTypes;
        Vector< Vector< MSI::Dof_Type > > tPDofTypes = { { MSI::Dof_Type::P } };

        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };

                // fill evaluation point xi, tau
                tParamPoint = { { 0.35 }, { -0.25 }, { 0.0 } };

                // number of coefficients
                tNumCoeffs = { { 4 }, { 9 }, { 16 } };

                // set the normal
                tNormal = { { 1.0 }, { 0.0 } };

                // set dof type list
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };

                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0, 0.0 },
                    { 1.0, 0.0, 0.0 },
                    { 1.0, 1.0, 0.0 },
                    { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 },
                    { 1.0, 0.0, 1.0 },
                    { 1.0, 1.0, 1.0 },
                    { 0.0, 1.0, 1.0 } };

                // fill evaluation point xi, tau
                tParamPoint = { { 0.35 }, { -0.25 }, { 0.75 }, { 0.0 } };

                // number of coefficients
                tNumCoeffs = { { 8 }, { 27 }, { 64 } };

                // set the normal
                tNormal = { { 1.0 }, { 0.0 }, { 0.0 } };

                // set dof type list
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };

                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // space and time geometry interpolators
        //------------------------------------------------------------------------------
        // create a space geometry interpolation rule
        mtk::Interpolation_Rule tGIRule( tGeometryType,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set the evaluation point
        tGI.set_space_time( tParamPoint );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder;

            // create random coefficients for leader FI
            arma::Mat< double > tLeaderMatrixVel;
            arma::Mat< double > tLeaderMatrixP;
            arma::Mat< double > tFollowerMatrixVel;
            arma::Mat< double > tFollowerMatrixP;

            // get number of dof
            int tNumDofVel = 0;
            int tNumDofP   = 0;

            // switch on interpolation order
            switch ( iInterpOrder )
            {
                case ( 1 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::LINEAR;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 0 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 0 );

                    // create random coefficients for leader FI
                    tLeaderMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tLeaderMatrixP.randu( tNumCoeffs( 0 ), 1 );

                    // create random coefficients for follower FI
                    tFollowerMatrixVel.randu( tNumCoeffs( 0 ), iSpaceDim );
                    tFollowerMatrixP.randu( tNumCoeffs( 0 ), 1 );
                    break;
                }
                case ( 2 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 1 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 1 );

                    // create random coefficients for leader FI
                    tLeaderMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tLeaderMatrixP.randu( tNumCoeffs( 1 ), 1 );

                    // create random coefficients for follower FI
                    tFollowerMatrixVel.randu( tNumCoeffs( 1 ), iSpaceDim );
                    tFollowerMatrixP.randu( tNumCoeffs( 1 ), 1 );
                    break;
                }
                case ( 3 ):
                {
                    // set interpolation type
                    tInterpolationOrder = mtk::Interpolation_Order::CUBIC;

                    // get number of dof
                    tNumDofVel = tNumCoeffs( 2 ) * iSpaceDim;
                    tNumDofP   = tNumCoeffs( 2 );

                    // create random coefficients for leader FI
                    tLeaderMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tLeaderMatrixP.randu( tNumCoeffs( 2 ), 1 );

                    // create random coefficients for follower FI
                    tFollowerMatrixVel.randu( tNumCoeffs( 2 ), iSpaceDim );
                    tFollowerMatrixP.randu( tNumCoeffs( 2 ), 1 );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "LINEAR, QUADRATIC or CUBIC only." );
                    break;
                }
            }

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::CONSTANT,
                    mtk::Interpolation_Order::CONSTANT );

            // fill random coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            tLeaderDOFHatVel = 10.0 * tLeaderMatrixVel;
            Matrix< DDRMat > tLeaderDOFHatP;
            tLeaderDOFHatP = 10.0 * tLeaderMatrixP;

            // fill random coefficients for follower FI
            Matrix< DDRMat > tFollowerDOFHatVel;
            tFollowerDOFHatVel = 10.0 * tFollowerMatrixVel;
            Matrix< DDRMat > tFollowerDOFHatP;
            tFollowerDOFHatP = 10.0 * tFollowerMatrixP;

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( 2 );

            // create the field interpolator
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );

            // set the coefficients uHat
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // set the evaluation point xi, tau
            tLeaderFIs( 0 )->set_space_time( tParamPoint );
            tLeaderFIs( 1 )->set_space_time( tParamPoint );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tFollowerFIs( 2 );

            // create the field interpolator
            tFollowerFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tFollowerFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );

            // set the coefficients uHat
            tFollowerFIs( 0 )->set_coeff( tFollowerDOFHatVel );
            tFollowerFIs( 1 )->set_coeff( tFollowerDOFHatP );

            // set the evaluation point xi, tau
            tFollowerFIs( 0 )->set_space_time( tParamPoint );
            tFollowerFIs( 1 )->set_space_time( tParamPoint );

            // define the IWGs
            fem::IWG_Factory tIWGFactory;

            std::shared_ptr< fem::IWG > tIWG =
                    tIWGFactory.create_IWG( fem::IWG_Type::FS_STRUC_INTERFACE );
            tIWG->set_residual_dof_type( tVelDofTypes );
            tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) }, mtk::Leader_Follower::LEADER );
            tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) }, mtk::Leader_Follower::FOLLOWER );

            // create and set the fem set for the IWG
            MSI::Equation_Set* tSet = new fem::Set();
            static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );

            // set pointer for IWG
            tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

            // set size for the set EqnObjDofTypeList
            tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

            // set size and populate the set dof type map
            tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

            // set size and populate the set leader and follower dof type map
            tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

            tIWG->mSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
            tIWG->mSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 4 );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofVel + tNumDofP, 2 * tNumDofVel + tNumDofP - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { 2 * tNumDofVel + tNumDofP, 2 * tNumDofVel + 2 * tNumDofP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVel - 1 },
                { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                { tNumDofVel + tNumDofP, 2 * tNumDofVel + tNumDofP - 1 },
                { 2 * tNumDofVel + tNumDofP, 2 * tNumDofVel + 2 * tNumDofP - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 4 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;

            // set IWG normal
            tIWG->set_normal( tNormal );

            // build global property type list
            tIWG->build_global_dof_dv_and_field_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes   = { tVelDofTypes, tPDofTypes };
            tIWG->mRequestedFollowerGlobalDofTypes = { tVelDofTypes, tPDofTypes };

            // create a field interpolator manager
            Vector< Vector< enum MSI::Dof_Type > > tDummy;
            Field_Interpolator_Manager             tLeaderFIManager( tDummy, tSet, mtk::Leader_Follower::LEADER );
            Field_Interpolator_Manager             tFollowerFIManager( tDummy, tSet, mtk::Leader_Follower::FOLLOWER );

            // populate the field interpolator manager
            tLeaderFIManager.mFI                       = tLeaderFIs;
            tLeaderFIManager.mIPGeometryInterpolator   = &tGI;
            tLeaderFIManager.mIGGeometryInterpolator   = &tGI;
            tFollowerFIManager.mFI                     = tFollowerFIs;
            tFollowerFIManager.mIPGeometryInterpolator = &tGI;
            tFollowerFIManager.mIGGeometryInterpolator = &tGI;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tLeaderFIManager );
            tIWG->set_field_interpolator_manager( &tFollowerFIManager, mtk::Leader_Follower::FOLLOWER );

            // reset residual and jacobian
            //------------------------------------------------------------------------------
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * ( tNumDofVel + tNumDofP ), 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * ( tNumDofVel + tNumDofP ), 2 * ( tNumDofVel + tNumDofP ), 0.0 );

            // check evaluation of the residual
            //------------------------------------------------------------------------------
            // evaluate the residual
            tIWG->compute_residual( 1.0 );

            // check evaluation of the jacobian  by FD
            //------------------------------------------------------------------------------
            // init the jacobian for IWG and FD evaluation
            Matrix< DDRMat > tJacobians;
            Matrix< DDRMat > tJacobiansFD;

            // check jacobian by FD
            bool tCheckJacobian = tIWG->check_jacobian_double(
                    tPerturbation,
                    tEpsilon,
                    1.0,
                    tJacobians,
                    tJacobiansFD );

            //            // print for debug
            //            print( tJacobians,"tJacobians");
            //            print( tJacobiansFD,"tJacobiansFD");

            // require check is true
            REQUIRE( tCheckJacobian );

            //            // print the treated case
            //            std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<std::endl;

            // clean up
            tLeaderFIs.clear();
            tFollowerFIs.clear();
        }
    }

} /* END_TEST_CASE */
