/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_L2_Damage_Bulk.cpp
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
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
#include "cl_FEM_IWG_L2_Damage_Bulk.hpp"
#undef protected
#undef private
// LINALG/src
#include "op_equal_equal.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

void Test_IWG_L2_Damage_Bulk(
        fem::IWG_Type               aIWGType,
        Vector< Matrix< DDRMat > >& aParameters )
{
    // define an epsilon environment
    real tEpsilon = 5E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

    // create list of interpolation orders
    moris::Vector< mtk::Interpolation_Order > tInterpolationOrders = {
        mtk::Interpolation_Order::LINEAR,
        mtk::Interpolation_Order::QUADRATIC,
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Vector< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    moris::Vector< moris::Vector< MSI::Dof_Type > > tDisplDofTypes      = { { MSI::Dof_Type::UX } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tDamageDofTypes     = { { MSI::Dof_Type::L2 } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tDofTypes           = { tDisplDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ), tDamageDofTypes( 0 ) };

    // init IWG
    //------------------------------------------------------------------------------

    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 0 - no smoothing
    moris::Vector< Matrix< DDRMat > > tDamageParameters = { //
        { { 0.0 } },                                        //
        { { 2.0, 1.0e-3, 10.0 } },                          //
        { { 0.0 } }
    };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    std::shared_ptr< fem::Property > tPropWeight = std::make_shared< fem::Property >();
    tPropWeight->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropLump = std::make_shared< fem::Property >();
    tPropLump->set_parameters( { { { 0.5 } } } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();
    tCMLeader->set_parameters( tDamageParameters );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( aIWGType );
    tIWG->set_residual_dof_type( tDamageDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_property( tPropWeight, "Weight" );
    tIWG->set_property( tPropLump, "Lump" );
    tIWG->set_constitutive_model( tCMLeader, "ElasticDamage" );
    if ( aIWGType == fem::IWG_Type::L2_EQSTRAIN_BULK )
    {
        tIWG->set_parameters( aParameters );
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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::L2 ) )         = 3;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::L2 ) )         = 3;

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

                // fill space coeff xHat
                tXHat = { { 0.0, 0.0 },
                    { 1.0, 0.0 },
                    { 1.0, 1.0 },
                    { 0.0, 1.0 } };

                // set displacement dof types
                tDisplDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
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

                // set displacement dof types
                tDisplDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };
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

        // set space dimension to CM, SP
        tCMLeader->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {

            // set parameters for IWG
            tIWG->set_parameters( { { { 2.0 } }, { { static_cast< real >( iInterpOrder ) } } } );

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
            int tNumDofDispl      = tNumCoeff * iSpaceDim;
            int tNumDofNlEqStrain = tNumCoeff;
            int tNumDofHistory    = tNumCoeff;
            int tNumDofDamage	  = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatHistory;
            fill_phat_Elast( tLeaderDOFHatHistory, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatDamage;
            fill_phat_Elast( tLeaderDOFHatDamage, iSpaceDim, iInterpOrder );

            // create a Vector of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDisplDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatHistory );

            // create the field interpolator damage
            tLeaderFIs( 3 ) = new Field_Interpolator( 1, tFIRule, &tGI, tDamageDofTypes( 0 ) );
            tLeaderFIs( 3 )->set_coeff( tLeaderDOFHatDamage );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDispl - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory + tNumDofDamage - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDispl - 1 },
                { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 },
                { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 },
                { tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory + tNumDofDamage - 1}
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory + tNumDofDamage, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory + tNumDofDamage, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory + tNumDofDamage, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Vector< moris::Vector< enum gen::PDV_Type > >   tDummyDv;
            moris::Vector< moris::Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

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
                    std::cout << "Case: Geometry " << iSpaceDim
                              << " Order " << iInterpOrder
                              << " iGP " << iGP << '\n';
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}

TEST_CASE( "IWG_L2_EqStrain_Bulk_Order1", "[moris],[fem],[IWG_L2_EqStrain_Bulk_Order1]" )
{
    moris::Vector< Matrix< DDRMat > > tParameters = { { { 0.25 } }, { { 1.0 } } };

    Test_IWG_L2_Damage_Bulk( fem::IWG_Type::L2_EQSTRAIN_BULK, tParameters );
} /*END_TEST_CASE*/

TEST_CASE( "IWG_L2_EqStrain_Bulk_Order2", "[moris],[fem],[IWG_L2_EqStrain_Bulk_Order2]" )
{
    moris::Vector< Matrix< DDRMat > > tParameters = { { { 0.25 } }, { { 2.0 } } };

    Test_IWG_L2_Damage_Bulk( fem::IWG_Type::L2_EQSTRAIN_BULK, tParameters );
} /*END_TEST_CASE*/

TEST_CASE( "IWG_L2_EqStrain_Bulk_Order3", "[moris],[fem],[IWG_L2_EqStrain_Bulk_Order3]" )
{
    moris::Vector< Matrix< DDRMat > > tParameters = { { { 0.25 } }, { { 3.0 } } };

    Test_IWG_L2_Damage_Bulk( fem::IWG_Type::L2_EQSTRAIN_BULK, tParameters );
} /*END_TEST_CASE*/

TEST_CASE( "IWG_L2_History_Bulk", "[moris],[fem],[IWG_L2_History_Bulk]" )
{
    moris::Vector< Matrix< DDRMat > > tParameters;

    Test_IWG_L2_Damage_Bulk( fem::IWG_Type::L2_HISTORY_BULK, tParameters );
} /*END_TEST_CASE*/

TEST_CASE( "IWG_L2_Damage_Bulk", "[moris],[fem],[IWG_L2_Damage_Bulk]" )
{
    moris::Vector< Matrix< DDRMat > > tParameters;

    Test_IWG_L2_Damage_Bulk( fem::IWG_Type::L2_DAMAGE_BULK, tParameters );
} /*END_TEST_CASE*/

//---------------------------------------------------------------------------------------------

/*END_TEST_CASE*/
