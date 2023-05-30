/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Dof.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private
// FEM/INT/src
#include "cl_FEM_Enums.hpp"    //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Design_Variable_Interface_Proxy.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

// FEM/MSI/src
#include "cl_MSI_Design_Variable_Interface.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
// LINALG/src
#include "op_equal_equal.hpp"

using namespace moris;
using namespace fem;

void
UT_FEM_IQI_Dof( uint aTestConfiguration )
{
    // define an epsilon environment
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
        mtk::Interpolation_Order::LINEAR,
        mtk::Interpolation_Order::QUADRATIC,
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes   = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tVelDofTypes( 0 ), tPDofTypes( 0 ) };

    // define the IQI
    fem::IQI_Factory tIQIFactory;

    std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::DOF );

    tIQI->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIQI->set_quantity_dof_type( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } } );
    tIQI->set_output_type_index( 1 );    // get vy component
    tIQI->set_name( "DofSet" );

    switch ( aTestConfiguration )
    {
        // no parameters are set - IQI is just vy value
        case 0:
        {
            break;
        }
        // IQI: 1st spatial derivative of vy in y-direction
        case 1:
        {
            tIQI->set_parameters( { { { 1, 1 } } } );
            break;
        }
        // IQI: 1st time derivative of vy
        case 2:
        {
            tIQI->set_parameters( { { { 0, 0 } }, { { 1 } } } );
            break;
        }
    }

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIQI->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIQI->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

    // set size and populate the set leader dof type map
    tIQI->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIQI->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } };
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
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

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
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofP   = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            ;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatP;
            fill_phat( tLeaderDOFHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator pressure
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatP );

            // set size and fill the set residual assembly map
            tIQI->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIQI->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIQI->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };

            // set size and init the set residual and jacobian
            tIQI->mSet->mResidual.resize( 1 );
            tIQI->mSet->mResidual( 0 ).set_size(
                    tNumDofVel + tNumDofP,
                    1,
                    0.0 );

            // fill requested IQI map
            moris::map< std::string, moris_index > tRequestedIQINamesAssemblyMap;
            tRequestedIQINamesAssemblyMap[ "DofSet" ] = 0;
            tIQI->mSet->mRequestedIQINamesAssemblyMap = tRequestedIQINamesAssemblyMap;

            // set size and init the set mQI
            tIQI->mSet->mQI.resize( 1 );
            tIQI->mSet->mQI( 0 ).set_size( 1, 1, 0.0 );

            // build global dof type list
            tIQI->get_global_dof_type_list();

            // populate the requested leader dof type
            tIQI->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;

            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIQI->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIQI->set_field_interpolator_manager( &tFIManager );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
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

                // reset jacobian
                tIQI->mSet->mResidual( 0 ).fill( 0.0 );

                Matrix< DDRMat > tdQIdu;
                Matrix< DDRMat > tdQIduFD;
                bool             tCheckdQIdu = tIQI->check_dQIdu_FD(
                        1.0,
                        tPerturbation,
                        tEpsilon,
                        tdQIdu,
                        tdQIduFD,
                        true );

                // print for debug
                if ( !tCheckdQIdu )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                }

                REQUIRE( tCheckdQIdu );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}

TEST_CASE( "IQI_Dof_Value", "[moris],[fem],[IQI_Dof_Value]" )
{
    // test Dof value
    UT_FEM_IQI_Dof( 0 );
}

TEST_CASE( "IQI_Dof_Spatial_Derivative", "[moris],[fem],[IQI_Dof_Spatial_Derivative]" )
{
    // test Dof spatial derivative
    UT_FEM_IQI_Dof( 1 );
}

TEST_CASE( "IQI_Dof_Time_Derivative", "[moris],[fem],[IQI_Dof_Time_Derivative]" )
{
    // test Dof temporal derivative
    UT_FEM_IQI_Dof( 2 );
}
/*END_TEST_CASE*/

