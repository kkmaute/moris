/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_SP_Crosswind.cpp
 *
 */

#include "catch.hpp"

#define protected public
#define private   public
//FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private

//LINALG/src
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "SP_Crosswind", "[SP_Crosswind]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // set geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

    // create list of interpolation orders
    moris::Cell< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC,
            mtk::Interpolation_Order::CUBIC };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< MSI::Dof_Type > tTargetDofTypes = { MSI::Dof_Type::VX };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes( 0 ), tTargetDofTypes };

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPCrosswind =
            tSPFactory.create_SP( fem::Stabilization_Type::CROSSWIND );
    //moris::Cell< std::string > tCellStringDof = {"Velocity","Target"};
    //tSPCrosswind->set_dof_type_list( tDofTypes, tCellStringDof );
    tSPCrosswind->set_parameters( {{{ 0.7 }}} );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPCrosswind->set_cluster( tCluster );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tSPCrosswind->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tSPCrosswind->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tSPCrosswind->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPCrosswind->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tSPCrosswind->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set leader dof type map
    tSPCrosswind->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPCrosswind->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tSPCrosswind->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set space dim
        tSPCrosswind->set_space_dim( iSpaceDim );

        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tSPCrosswind->set_normal( tNormal );

        // switch on space dimension
        switch( iSpaceDim )
        {
            case 2 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0 },
                        { 1.0, 0.0 },
                        { 1.0, 1.0 },
                        { 0.0, 1.0 }};

                // set velocity dof types
                tVelDofTypes = { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } };
                break;
            }
            case 3 :
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;

                // fill space coeff xHat
                tXHat = {{ 0.0, 0.0, 0.0 },
                        { 1.0, 0.0, 0.0 },
                        { 1.0, 1.0, 0.0 },
                        { 0.0, 1.0, 0.0 },
                        { 0.0, 0.0, 1.0 },
                        { 1.0, 0.0, 1.0 },
                        { 1.0, 1.0, 1.0 },
                        { 0.0, 1.0, 1.0 }};

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
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
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
                    mtk::Integration_Order::BAR_2 );

            // create an integrator
            mtk::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatVel;
            fill_uhat( tLeaderDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_phat( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatVel );

            // create the field interpolator temperature
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTargetDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDOFHatTEMP );

            // build global dof type list
            tSPCrosswind->get_global_dof_type_list();

            // set order
            tSPCrosswind->set_interpolation_order( iInterpOrder );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tSPCrosswind->mSet->mLeaderFIManager = &tFIManager;

            // set SP field interpolator manager
            tSPCrosswind->set_field_interpolator_manager( &tFIManager );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset SP evaluation flags
                tSPCrosswind->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tSPCrosswind->mSet->mLeaderFIManager->set_space_time( tParamPoint );

                // evaluate the SP for debug
                //print(tSPCrosswind->val(),"tSPCrosswind");

                // populate the requested leader dof type for SP
                moris::Cell< moris::Cell< MSI::Dof_Type > > tLeaderDofTypes =
                        tSPCrosswind->get_global_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tLeaderDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tLeaderDofTypes( iRequestedDof );

                    // print derivative dof for debug
                    // std::cout<<"Derivative dof type treated "<<static_cast<uint>(tDofDerivative(0))<<std::endl;

                    // evaluate dspdu
                    Matrix< DDRMat > tdspdu = tSPCrosswind->dSPdLeaderDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdspduFD;
                    tSPCrosswind->eval_dSPdLeaderDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                    // print for debug
                    // print( tdspdu, "tdspdu" );
                    // print( tdspduFD, "tdspduFD" );

                    // check that analytical and FD match
                    bool tCheckCrosswind = fem::check( tdspdu, tdspduFD, tEpsilon );
                    REQUIRE( tCheckCrosswind );
                }
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}

