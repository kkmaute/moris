/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Linear_Bulk.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private public
// FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private
// MTK/src
#include "cl_MTK_Enums.hpp"
// LINALG/src
#include "op_equal_equal.hpp"
// FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

void IWG_Time_Continuity_Dof_Check(
        const real& aEndTime,
        const real& aLump )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    Vector< mtk::Interpolation_Order > tInterpolationOrders = {
        mtk::Interpolation_Order::LINEAR,
        mtk::Interpolation_Order::QUADRATIC,
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = tDispDofTypes;

    // init IWG
    std::shared_ptr< fem::Property > tPropWeightCurrent = std::make_shared< fem::Property >();
    tPropWeightCurrent->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropWeightPrevious = std::make_shared< fem::Property >();
    tPropWeightPrevious->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropInitialCondition = std::make_shared< fem::Property >();

    std::shared_ptr< fem::Property > tPropLump = std::make_shared< fem::Property >();
    tPropLump->set_parameters( { { { aLump } } } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );

    tIWG->set_property( tPropWeightCurrent, "WeightCurrent" );
    tIWG->set_property( tPropWeightPrevious, "WeightPrevious" );
    tIWG->set_property( tPropInitialCondition, "InitialCondition" );

    // omit lump parameter for aLump < 0
    if ( aLump >= 0 )
    {
        tIWG->set_property( tPropLump, "Lump" );
    }
    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                tPropInitialCondition->set_parameters( { { { 1.0 }, { 2.0 } } } );

                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                tPropInitialCondition->set_parameters( { { { 1.0 }, { 2.0 }, { 3.0 } } } );

                break;
            }
            default:
            {
                MORIS_ERROR( false, " QUAD or HEX only." );
                break;
            }
        }

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tGIInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            // create a space geometry interpolation rule
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tGIInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // create a space time geometry interpolator
            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            // create time coeff tHat
            Matrix< DDRMat > tTHat = { { 0.0 }, { aEndTime } };

            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );

            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            mtk::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            mtk::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof per type
            int tNumDofDisp = tNumCoeff * iSpaceDim;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
            Vector< Field_Interpolator* > tLeaderPreviousFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            tLeaderPreviousFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderPreviousFIs( 0 )->set_coeff( 0.5 * tLeaderDofHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp, tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager               tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager               tPreviousFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            tPreviousFIManager.mFI                     = tLeaderPreviousFIs;
            tPreviousFIManager.mIPGeometryInterpolator = &tGI;
            tPreviousFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager         = &tFIManager;
            tIWG->mSet->mLeaderPreviousFIManager = &tPreviousFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );
            tIWG->set_field_interpolator_manager_previous_time( &tPreviousFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mLeaderPreviousFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWG->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWG->mSet->mJacobian.fill( 0.0 );

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                // check jacobian by FD
                bool tCheckJacobian = tIWG->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                // print for debug
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << '\n';
                }

                // require check is true
                REQUIRE( tCheckJacobian );

                // check derivative with respect to previous dofs - only if not fist time step
                //------------------------------------------------------------------------------
                if ( aEndTime > 0 )
                {
                    tIWG->mSet->mJacobian.fill( 0.0 );

                    Matrix< DDRMat > tJacFD = tIWG->mSet->mJacobian;

                    tIWG->compute_jacobian_previous( 1.0 );

                    Matrix< DDRMat > tJacAna = tIWG->mSet->mJacobian;

                    Matrix< DDRMat > tCoefficients = tLeaderPreviousFIs( 0 )->get_coeff();

                    Matrix< DDRMat > tFDpertb = { { 1.0, -1.0 } };

                    uint tRowCtr = 0;

                    for ( uint iCol = 0; iCol < tCoefficients.n_cols(); ++iCol )
                    {
                        for ( uint iRow = 0; iRow < tCoefficients.n_rows(); ++iRow )
                        {
                            for ( uint iWgt = 0; iWgt < 2; ++iWgt )
                            {
                                tCoefficients( iRow, iCol ) += tFDpertb( iWgt ) * tPerturbation;

                                tLeaderPreviousFIs( 0 )->set_coeff( tCoefficients );

                                tLeaderPreviousFIs( 0 )->reset_eval_flags();

                                tIWG->reset_eval_flags();

                                tIWG->mSet->get_residual()( 0 ).fill( 0.0 );

                                tIWG->compute_residual( 1.0 );

                                tJacFD( { tRowCtr, tRowCtr }, { 0, tJacFD.n_rows() - 1 } ) +=    //
                                        tFDpertb( iWgt ) / ( 2.0 * tPerturbation ) * trans( tIWG->mSet->get_residual()( 0 ) );

                                tCoefficients( iRow, iCol ) -= tFDpertb( iWgt ) * tPerturbation;
                            }
                            tRowCtr++;
                        }
                    }

                    // check results
                    tCheckJacobian = true;

                    for ( uint iCol = 0; iCol < tJacFD.n_cols(); ++iCol )
                    {
                        for ( uint iRow = 0; iRow < tJacFD.n_rows(); ++iRow )
                        {
                            if ( std::abs( tJacFD( iRow, iCol ) - tJacAna( iRow, iCol ) ) > tEpsilon )
                            {
                                tCheckJacobian = false;

                                std::cout << "iiJacPrevious " << iRow << " - jjJacPrevious " << iCol << "\n"
                                          << std::flush;
                                std::cout << "aJacobianPrevious( iiJac, jjJac )   " << std::setprecision( 12 ) << tJacAna( iRow, iCol ) << "\n"
                                          << std::flush;
                                std::cout << "aJacobianFDPrevious( iiJac, jjJac ) " << std::setprecision( 12 ) << tJacFD( iRow, iCol ) << "\n"
                                          << std::flush;
                                std::cout << "Absolute difference " << std::abs( tJacFD( iRow, iCol ) - tJacAna( iRow, iCol ) ) << "\n"
                                          << std::flush;
                            }
                        }
                    }

                    // print for debug
                    if ( !tCheckJacobian )
                    {
                        std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << '\n';
                    }

                    // require check is true
                    REQUIRE( tCheckJacobian );
                }
            }

            // clean up
            tLeaderFIs.clear();
            tLeaderPreviousFIs.clear();
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Initial_Condition", "[IWG_Time_Continuity_Initial_Condition]" )
{
    // define end time for time slab; if zero: initial time step
    real tEndTime = 0.0;
    real tLump    = 0.0;

    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = 1.0;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = 0.7;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = -1.0;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );
}

TEST_CASE( "IWG_Time_Continuity_Dof_Previous_Time", "[IWG_Time_Continuity_Initial_Condition_Previous_Time]" )
{
    // define end time for time slab; if zero: initial time step
    real tEndTime = 1.0;
    real tLump    = 0.0;

    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = 1.0;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = 0.7;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );

    tLump = -1.0;
    IWG_Time_Continuity_Dof_Check( tEndTime, tLump );
}
/*END_TEST_CASE*/
