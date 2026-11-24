/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Max_Damage.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private   public
//FEM/INT/src
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Design_Variable_Interface_Proxy.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

//FEM/MSI/src
#include "cl_MSI_Design_Variable_Interface.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IQI_Max_Damage","[moris],[fem],[IQI_Max_Damage]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Vector< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};
    Matrix< DDRMat > tNumHalfCoeffs = { { 4, 9, 16 }, { 8, 27, 64 } };


    // dof type list
    moris::Vector< moris::Vector< MSI::Dof_Type > > tDisplDofTypes      = { { MSI::Dof_Type::UX } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tNlEqStrainDofTypes = { { MSI::Dof_Type::NLEQSTRAIN } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tHistoryDofTypes    = { { MSI::Dof_Type::HISTORY } };
    moris::Vector< moris::Vector< MSI::Dof_Type > > tDofTypes           = { tDisplDofTypes( 0 ), tNlEqStrainDofTypes( 0 ), tHistoryDofTypes( 0 ) };


    // init IQI
    //------------------------------------------------------------------------------
    // create the properties
    // Damage setup - 0,2,0
    // Local equivalent strain - 0 - energy release rate
    // Damage law - 2 - smooth exponential
    // Smoothing law - 2 - ks smoothing
    moris::Vector< Matrix< DDRMat > > tDamageParameters = { //
        { { 2.0, 0.818 } },                                 //
        { { 1.0, 1.0e-4, 0.95, 100.0 } },                   //
        { { 2.0, 5.0 } }
    };

    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 200.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    std::shared_ptr< fem::Property > tPropReferenceValue = std::make_shared< fem::Property >();
    tPropReferenceValue->set_parameters( { { { 0.5 } } } );

    std::shared_ptr< fem::Property > tPropExponent = std::make_shared< fem::Property >();
    tPropExponent->set_parameters( { { { 2.0 } } } );

    std::shared_ptr< fem::Property > tPropShift = std::make_shared< fem::Property >();
    tPropShift->set_parameters( { { { 1.0 } } } );


    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeader =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO_DAMAGE);
    tCMLeader->set_model_type( fem::Model_Type::PLANE_STRESS );
    tCMLeader->set_dof_type_list( tDofTypes );
    tCMLeader->set_property( tPropEMod, "YoungsModulus" );
    tCMLeader->set_property( tPropNu, "PoissonRatio" );
    tCMLeader->set_local_properties();

    // define the IQI
    fem::IQI_Factory tIQIFactory;

    std::shared_ptr< fem::IQI > tIQI =
            tIQIFactory.create_IQI( fem::IQI_Type::MAX_DAMAGE);
    tIQI->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIQI->set_property( tPropReferenceValue, "ReferenceValue" );
    tIQI->set_property( tPropExponent, "Exponent" );
    tIQI->set_property( tPropShift, "Shift" );
    tIQI->set_constitutive_model( tCMLeader, "ElasticDamage" );
    tIQI->set_name("Max_Damage");

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIQI->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIQI->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;


    // set size and populate the set leader dof type map
    tIQI->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::NLEQSTRAIN ) ) = 1;
    tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::HISTORY ) )    = 2;


    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
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

                // set displacement dof types
                tDisplDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

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

        // set space dimension to CM, SP
        tCMLeader->set_space_dim( iSpaceDim );
        tCMLeader->set_parameters( tDamageParameters );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};
            Matrix< DDRMat > tXHat;
            fill_xhat( tXHat, iSpaceDim, iInterpOrder );

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
            int tNumDofDispl      = tNumCoeff * iSpaceDim;
            int tNumDofNlEqStrain = tNumCoeff;
            int tNumDofHistory    = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // get number of coefficients
            uint tNumHalfCoeffCurrent = tNumHalfCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            uint tNumCoeffCurrent     = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // create a cell of field interpolators current
            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator for displacement
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDisplDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatDispl;
            fill_uhat_Elast( tLeaderDOFHatDispl, iSpaceDim, iInterpOrder );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatDispl );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tNlEqStrainDofTypes( 0 ) );
            // fill coefficients for master FI
            Matrix< DDRMat > tLeaderDOFHatNlEqStrain;
            fill_phat_Elast( tLeaderDOFHatNlEqStrain, iSpaceDim, iInterpOrder );
            tLeaderDOFHatNlEqStrain = 1.0e-3 * tLeaderDOFHatNlEqStrain;
            tLeaderDOFHatNlEqStrain( { tNumHalfCoeffCurrent, tNumCoeffCurrent - 1 }, { 0, 0 } ) =
                    2.0 * tLeaderDOFHatNlEqStrain( { 0, tNumHalfCoeffCurrent - 1 }, { 0, 0 } );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // create the field interpolator for nonlocal equivalent strain
            tLeaderFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tHistoryDofTypes( 0 ) );
            // fill coefficients for master FI
            tLeaderFIs( 2 )->set_coeff( tLeaderDOFHatNlEqStrain );

            // set size and fill the set residual assembly map
            tIQI->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIQI->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDispl - 1 } };
            tIQI->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDispl, tNumDofDispl + tNumDofNlEqStrain - 1 } };
            tIQI->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofDispl + tNumDofNlEqStrain, tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory - 1 } };

            // set size and init the set residual
            tIQI->mSet->mResidual.resize( 1 );
            tIQI->mSet->mResidual( 0 ).set_size( tNumDofDispl + tNumDofNlEqStrain + tNumDofHistory, 1, 0.0 );

            // fill requested IQI map
            moris::map< std::string, moris_index > tRequestedIQINamesAssemblyMap;
            tRequestedIQINamesAssemblyMap[ "Max_Damage" ] = 0;
            tIQI->mSet->mRequestedIQINamesAssemblyMap = tRequestedIQINamesAssemblyMap;

            // set size and init the set mQI
            tIQI->mSet->mQI.resize( 1 );
            tIQI->mSet->mQI( 0 ).set_size( 1, 1, 0.0 );

            // build global dof type list
            tIQI->get_global_dof_type_list();

            // populate the requested leader dof type
            tIQI->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            moris::Vector< moris::Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIQI->mSet->mLeaderFIManager = &tFIManager;

            // set IQI field interpolator manager
            tIQI->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IQI evaluation flags
                tIQI->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIQI->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // check evaluation of the quantity of interest
                //------------------------------------------------------------------------------
                // reset residual
                tIQI->mSet->mQI( 0 ).fill( 0.0 );

                // compute residual
                tIQI->compute_QI( 1.0 );

                // check evaluation of the quantity of interest derivative
                //------------------------------------------------------------------------------

                // reset residual
                tIQI->mSet->mResidual( 0 ).fill( 0.0 );

                Matrix< DDRMat > tdQIdu;
                Matrix< DDRMat > tdQIduFD;
                bool tCheckdQIdu = tIQI->check_dQIdu_FD(
                        1.0,
                        tPerturbation,
                        tEpsilon,
                        tdQIdu,
                        tdQIduFD,
                        true);

                // print( tdQIdu, "tdQIdu" );
                // print( tdQIduFD, "tdQIduFD" );

                // print for debug
                if( !tCheckdQIdu )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << '\n';
                }

                REQUIRE( tCheckdQIdu );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}/*END_TEST_CASE*/


