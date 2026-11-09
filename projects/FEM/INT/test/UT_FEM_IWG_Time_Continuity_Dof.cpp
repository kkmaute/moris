/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Time_Continuity_Dof.cpp
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
#include "fn_norm.hpp"

// FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

namespace
{
    struct Time_Continuity_Test_Options
    {
        bool mConstantTime             = false;
        bool mScalar                   = false;
        bool mUniformFields            = false;
        bool mCurrentWeightDependency  = false;
        bool mPreviousWeightDependency = false;
        real mCurrentTimeCoordinate    = -1.0;
        real mPreviousTimeCoordinate   = 1.0;
        real mIntegrationWeight        = 1.0;
        real mThickness                = 1.0;
        real mResidualWeight           = 1.0;
        uint mMaxSpaceOrder            = 3;
    };
}    // namespace

void IWG_Time_Continuity_Dof_Check(
        const real&                         aInterfaceTime,
        const real&                         aLump,
        const Time_Continuity_Test_Options& aOptions = {} )
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

    // Weight dependencies refer to the current solution even for WeightPrevious.
    if ( aOptions.mCurrentWeightDependency )
    {
        tPropWeightCurrent->set_dof_type_list( tDofTypes );
        tPropWeightCurrent->set_val_function( tUXFIValFunc_Elast );
        tPropWeightCurrent->set_dof_derivative_functions( { tUXFIDerFunc_Elast } );
    }

    if ( aOptions.mPreviousWeightDependency )
    {
        tPropWeightPrevious->set_dof_type_list( tDofTypes );
        tPropWeightPrevious->set_val_function( tUXFIValFunc_Elast );
        tPropWeightPrevious->set_dof_derivative_functions( { tUXFIDerFunc_Elast } );
    }

    auto tPropThickness = std::make_shared< fem::Property >();
    tPropThickness->set_parameters( { { { aOptions.mThickness } } } );

    auto tPropResidualWeight = std::make_shared< fem::Property >();
    tPropResidualWeight->set_parameters( { { { aOptions.mResidualWeight } } } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::TIME_CONTINUITY_DOF );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );

    tIWG->set_property( tPropWeightCurrent, "WeightCurrent" );
    tIWG->set_property( tPropWeightPrevious, "WeightPrevious" );
    tIWG->set_property( tPropInitialCondition, "InitialCondition" );

    if ( aOptions.mThickness != 1.0 )
    {
        tIWG->set_property( tPropThickness, "Thickness" );
    }

    if ( aOptions.mResidualWeight != 1.0 )
    {
        tIWG->set_property( tPropResidualWeight, "WeightResidual" );
    }

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

        const uint       tNumFields = aOptions.mScalar ? 1 : iSpaceDim;
        Matrix< DDRMat > tInitialValues( tNumFields, 1 );

        for ( uint iField = 0; iField < tNumFields; ++iField )
        {
            tInitialValues( iField ) = 1.0 + iField;
        }

        tPropInitialCondition->set_parameters( { tInitialValues } );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder <= aOptions.mMaxSpaceOrder; iInterpOrder++ )
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
            Geometry_Interpolator tGI( tGIRule );
            Geometry_Interpolator tPreviousGI( tGIRule );

            // create time coeff tHat
            Matrix< DDRMat > tTHat         = { { aInterfaceTime }, { aInterfaceTime + 1.0 } };
            Matrix< DDRMat > tPreviousTHat = { { aInterfaceTime - 1.0 }, { aInterfaceTime } };

            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );

            // set the coefficients xHat, tHat
            tGI.set_coeff( tXHat, tTHat );
            tPreviousGI.set_coeff( tXHat, tPreviousTHat );

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

            if ( aOptions.mConstantTime )
            {
                tNumCoeff /= 2;
            }

            // get number of dof per type
            int tNumDofDisp = tNumCoeff * tNumFields;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    aOptions.mConstantTime ? mtk::Interpolation_Type::CONSTANT : mtk::Interpolation_Type::LAGRANGE,
                    aOptions.mConstantTime ? mtk::Interpolation_Order::CONSTANT : mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );

            // The proxy supplies linear-time vector coefficients; select the needed block.
            Matrix< DDRMat > tSelectedCoefficients = tLeaderDofHatDisp(
                    { 0, tNumCoeff - 1 }, { 0, tNumFields - 1 } );
            tLeaderDofHatDisp = tSelectedCoefficients;

            if ( aOptions.mUniformFields )
            {
                for ( uint iField = 0; iField < tNumFields; ++iField )
                {
                    for ( uint iCoeff = 0; iCoeff < tNumCoeff; ++iCoeff )
                    {
                        tLeaderDofHatDisp( iCoeff, iField ) = 3.0 + iField;
                    }
                }
            }

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
            Vector< Field_Interpolator* > tLeaderPreviousFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( tNumFields, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );

            tLeaderPreviousFIs( 0 ) = new Field_Interpolator( tNumFields, tFIRule, &tPreviousGI, tDispDofTypes( 0 ) );
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
            tPreviousFIManager.mIPGeometryInterpolator = &tPreviousGI;
            tPreviousFIManager.mIGGeometryInterpolator = &tPreviousGI;

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
                INFO( "Dimension=" << iSpaceDim << " spatial order=" << iInterpOrder
                                   << " constant time=" << aOptions.mConstantTime << " scalar=" << aOptions.mScalar
                                   << " time=" << aInterfaceTime << " lump=" << aLump << " GP=" << iGP );
                tPropResidualWeight->set_parameters( { { { aOptions.mResidualWeight } } } );

                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // Default to matching endpoints; the temporal-block case uses interior coordinates.
                tParamPoint( iSpaceDim ) = aOptions.mCurrentTimeCoordinate;
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );
                tParamPoint( iSpaceDim ) = aOptions.mPreviousTimeCoordinate;
                tIWG->mSet->mLeaderPreviousFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWG->compute_residual( aOptions.mIntegrationWeight );

                // Uniform spatial fields give the same jump for every lumping fraction.
                // This reference uses partition of unity, independently of the lumping helper.
                if ( aOptions.mUniformFields )
                {
                    Matrix< DDRMat > tPreviousValue = aInterfaceTime > 0.0
                                                            ? Matrix< DDRMat >( tLeaderPreviousFIs( 0 )->val() )
                                                            : tInitialValues;
                    Matrix< DDRMat > tExpected      = aOptions.mIntegrationWeight * aOptions.mThickness * aOptions.mResidualWeight
                                               * tLeaderFIs( 0 )->N_trans()
                                               * ( tPropWeightCurrent->val()( 0 ) * tLeaderFIs( 0 )->val()
                                                       - tPropWeightPrevious->val()( 0 ) * tPreviousValue );

                    REQUIRE( norm( tIWG->mSet->mResidual( 0 ) - tExpected ) < tEpsilon );
                }

                // WeightResidual deliberately scales only the residual. Check that the
                // Jacobian is unchanged, then use unit residual weight for the FD check.
                if ( aOptions.mResidualWeight != 1.0 )
                {
                    tIWG->mSet->mJacobian.fill( 0.0 );
                    tIWG->compute_jacobian( aOptions.mIntegrationWeight );

                    Matrix< DDRMat > tWeightedJacobian = tIWG->mSet->mJacobian;
                    tPropResidualWeight->set_parameters( { { { 1.0 } } } );
                    tIWG->reset_eval_flags();
                    tIWG->mSet->mJacobian.fill( 0.0 );
                    tIWG->compute_jacobian( aOptions.mIntegrationWeight );

                    REQUIRE( norm( tWeightedJacobian - tIWG->mSet->mJacobian ) < tEpsilon );
                }

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
                        aOptions.mIntegrationWeight,
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
                if ( aInterfaceTime > 0 )
                {
                    // check_jacobian restores coefficients but leaves property caches
                    // from its final perturbation. Reevaluate at the restored state.
                    tIWG->reset_eval_flags();
                    tIWG->mSet->mJacobian.fill( 0.0 );

                    Matrix< DDRMat > tJacFD = tIWG->mSet->mJacobian;

                    tIWG->compute_jacobian_previous( aOptions.mIntegrationWeight );

                    Matrix< DDRMat > tJacAna = tIWG->mSet->mJacobian;

                    Matrix< DDRMat > tCoefficients = tLeaderPreviousFIs( 0 )->get_coeff();

                    Matrix< DDRMat > tFDpertb = { { 1.0, -1.0 } };

                    uint tColCtr = 0;

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

                                tIWG->compute_residual( aOptions.mIntegrationWeight );

                                tJacFD.get_column( tColCtr ) +=    //
                                        tFDpertb( iWgt ) / ( 2.0 * tPerturbation ) * tIWG->mSet->get_residual()( 0 );

                                tCoefficients( iRow, iCol ) -= tFDpertb( iWgt ) * tPerturbation;
                            }

                            tColCtr++;
                        }
                    }

                    // Restore the unperturbed previous state before the next point.
                    tLeaderPreviousFIs( 0 )->set_coeff( tCoefficients );
                    tIWG->reset_eval_flags();

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

                        print( tJacAna, "Jacobian Previous" );
                        print( tJacFD, "Jacobian Previous FD" );
                        print( tJacFD - tJacAna, "Jacobian Previous Diff" );
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

    // The stack-allocated managers have already been destroyed; Set must not delete them.
    tIWG->mSet->mLeaderFIManager         = nullptr;
    tIWG->mSet->mLeaderPreviousFIManager = nullptr;

    delete tSet;
}

TEST_CASE( "IWG_Time_Continuity_Dof_Initial_Condition", "[IWG_Time_Continuity_Initial_Condition],[IWG_Time_Continuity]" )
{
    // Vector jump at the initial interface; check current Jacobian for every lumping fraction.
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

TEST_CASE( "IWG_Time_Continuity_Dof_Previous_Time", "[IWG_Time_Continuity_Initial_Condition_Previous_Time],[IWG_Time_Continuity]" )
{
    // Consecutive vector-valued slabs; check current and previous Jacobians for every lumping fraction.
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

TEST_CASE( "IWG_Time_Continuity_Dof_Constant_Time", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // Scalar and vector jumps in 2D/3D with one temporal basis. Uniform fields
    // check the physical residual as well as both Jacobians under spatial lumping.
    Time_Continuity_Test_Options tOptions;
    tOptions.mConstantTime  = true;
    tOptions.mUniformFields = true;
    tOptions.mMaxSpaceOrder = 1;

    for ( bool tScalar : { false, true } )
    {
        tOptions.mScalar = tScalar;

        for ( real tTime : { 0.0, 1.0 } )
        {
            for ( real tLump : { -1.0, 0.0, 0.7, 1.0 } )
            {
                IWG_Time_Continuity_Dof_Check( tTime, tLump, tOptions );
            }
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Thickness_And_Residual_Weight", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // A nonunit thickness scales the physical jump and both Jacobians. Residual
    // weighting includes zero and fractional weights but must leave the Jacobian intact.
    Time_Continuity_Test_Options tOptions;
    tOptions.mUniformFields     = true;
    tOptions.mThickness         = 2.5;
    tOptions.mIntegrationWeight = 1.7;
    tOptions.mMaxSpaceOrder     = 1;

    for ( bool tConstantTime : { false, true } )
    {
        tOptions.mConstantTime = tConstantTime;

        for ( real tWeight : { 0.0, 0.4, 1.0 } )
        {
            tOptions.mResidualWeight = tWeight;

            for ( real tTime : { 0.0, 1.0 } )
            {
                for ( real tLump : { 0.0, 0.7, 1.0 } )
                {
                    IWG_Time_Continuity_Dof_Check( tTime, tLump, tOptions );
                }
            }
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Current_Weight_Dependency", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // WeightCurrent depends on the current field. FD checks exercise the product
    // rule at the initial interface and across consecutive slabs, including lumping.
    Time_Continuity_Test_Options tOptions;
    tOptions.mCurrentWeightDependency = true;
    tOptions.mMaxSpaceOrder           = 1;

    for ( bool tConstantTime : { false, true } )
    {
        tOptions.mConstantTime = tConstantTime;

        for ( real tTime : { 0.0, 1.0 } )
        {
            for ( real tLump : { 0.0, 0.7, 1.0 } )
            {
                IWG_Time_Continuity_Dof_Check( tTime, tLump, tOptions );
            }
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Previous_Weight_Dependency", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // WeightPrevious depends on the current field, not the previous coefficients.
    // Both FD Jacobians distinguish its product-rule term from the direct previous-state derivative.
    Time_Continuity_Test_Options tOptions;
    tOptions.mPreviousWeightDependency = true;
    tOptions.mMaxSpaceOrder            = 1;

    for ( bool tConstantTime : { false, true } )
    {
        tOptions.mConstantTime = tConstantTime;

        for ( real tTime : { 0.0, 1.0 } )
        {
            for ( real tLump : { 0.0, 0.7, 1.0 } )
            {
                DYNAMIC_SECTION( "constant time=" << tConstantTime << " interface=" << tTime << " lump=" << tLump )
                {
                    IWG_Time_Continuity_Dof_Check( tTime, tLump, tOptions );
                }
            }
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Scalar_Linear_Time", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // A scalar jump exercises the field/block indexing without vector components.
    // Check a physical residual reference and both Jacobians in 2D/3D.
    Time_Continuity_Test_Options tOptions;
    tOptions.mScalar        = true;
    tOptions.mUniformFields = true;
    tOptions.mMaxSpaceOrder = 1;

    for ( real tTime : { 0.0, 1.0 } )
    {
        for ( real tLump : { -1.0, 0.0, 0.7, 1.0 } )
        {
            IWG_Time_Continuity_Dof_Check( tTime, tLump, tOptions );
        }
    }
}

TEST_CASE( "IWG_Time_Continuity_Dof_Temporal_Blocks", "[IWG_Time_Continuity],[IWG_Time_Continuity_Extended]" )
{
    // Interior temporal coordinates keep both linear time bases active. This
    // exercises the off-diagonal temporal blocks in spatial-only lumping;
    // endpoint tests alone reduce those blocks to zero. Check both FD Jacobians.
    Time_Continuity_Test_Options tOptions;
    tOptions.mCurrentTimeCoordinate  = -0.3;
    tOptions.mPreviousTimeCoordinate = 0.4;
    tOptions.mMaxSpaceOrder          = 1;

    for ( real tLump : { 0.0, 0.7, 1.0 } )
    {
        IWG_Time_Continuity_Dof_Check( 1.0, tLump, tOptions );
    }
}

/*END_TEST_CASE*/
