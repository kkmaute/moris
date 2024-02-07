/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Diffusion_Robin_Nitsche.cpp
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
// FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Diffusion_UT.cpp"
// LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

inline void
tGeoValFunction_UTIWGDIFFDIR(
        moris::Matrix< moris::DDRMat >&                aPropMatrix,
        Vector< moris::Matrix< moris::DDRMat > >& aParameters,
        moris::fem::Field_Interpolator_Manager*        aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
}

using namespace moris;
using namespace fem;

inline void
UT_FEM_IWG_Diffusion_Robin_Nitsche_Core( enum fem::IWG_Type tIWGType )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // create space coeff xHat
    Matrix< DDRMat > tXHat;

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
    Vector< Vector< MSI::Dof_Type > > tTempDofTypes = { { MSI::Dof_Type::TEMP } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property >();
    tPropLeaderConductivity->set_parameters( { { { 2.0 } } } );
    tPropLeaderConductivity->set_val_function( tConstValFunc_Diff );
    //         tPropLeaderConductivity->set_dof_type_list( { tTempDofTypes } );
    //         tPropLeaderConductivity->set_val_function( tTEMPFIValFunc_Diff );
    //         tPropLeaderConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc_Diff } );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_parameters( { { { 1.0 } } } );
    tPropLeaderDirichlet->set_val_function( tConstValFunc_Diff );
    //         tPropLeaderDirichlet->set_dof_type_list( { tTempDofTypes } );
    //         tPropLeaderDirichlet->set_val_function( tTEMPFIValFunc_Diff );
    //         tPropLeaderDirichlet->set_dof_derivative_functions( { tTEMPFIDerFunc_Diff } );

    std::shared_ptr< fem::Property > tPropLeaderNemannPen = std::make_shared< fem::Property >();
    tPropLeaderNemannPen->set_parameters( { { { 0.1 } } } );
    tPropLeaderNemannPen->set_val_function( tConstValFunc_Diff );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiffLinIso->set_dof_type_list( { tTempDofTypes } );

    tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
    tCMLeaderDiffLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory                                 tSPFactory;
    std::shared_ptr< fem::Stabilization_Parameter > tSPRobinNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::ROBIN_NITSCHE );
    tSPRobinNitsche->set_parameters( { { { 10.0 } } } );
    tSPRobinNitsche->set_property( tPropLeaderNemannPen, "NeumannPenalty", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPRobinNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( tIWGType );
    tIWG->set_residual_dof_type( tTempDofTypes );
    tIWG->set_dof_type_list( { tTempDofTypes } );
    tIWG->set_stabilization_parameter( tSPRobinNitsche, "RobinNitsche" );
    tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderNemannPen, "NeumannPenalty" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

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
        tCMLeaderDiffLinIso->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
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
            int tNumDofTEMP = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FI
            Matrix< DDRMat > tLeaderDOFHatTEMP;
            fill_that( tLeaderDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Vector< Field_Interpolator* > tLeaderFIs( tTempDofTypes.size() );

            // create the field interpolator temperature
            tLeaderFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDOFHatTEMP );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tTempDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofTEMP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = { { 0, tNumDofTEMP - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( tTempDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofTEMP, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofTEMP, tNumDofTEMP, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tTempDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >        tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tTempDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mLeaderFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integartion points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tParamPoint );

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

                // This check does not apply anymore!
                if ( tIWGType == fem::IWG_Type::SPATIALDIFF_ROBIN_SYMMETRIC_NITSCHE )
                {
                    // real tRelError = norm( tJacobian - trans( tJacobian ) ) / norm( tJacobian );
                    // REQUIRE( tRelError < 1e-12 );
                }

                // print for debug
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << "iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}

TEST_CASE( "IWG_Diff_Robin_Symmetric", "[moris],[fem],[IWG_Diff_Robin_Symmetric]" )
{
    UT_FEM_IWG_Diffusion_Robin_Nitsche_Core( fem::IWG_Type::SPATIALDIFF_ROBIN_SYMMETRIC_NITSCHE );
}

TEST_CASE( "IWG_Diff_Robin_Unsymmetric", "[moris],[fem],[IWG_Diff_Robin_Unsymmetric]" )
{
    UT_FEM_IWG_Diffusion_Robin_Nitsche_Core( fem::IWG_Type::SPATIALDIFF_ROBIN_UNSYMMETRIC_NITSCHE );
}

TEST_CASE( "IWG_Diff_Robin_Geo_Prop", "[moris],[fem],[IWG_Diff_Robin_Geo_Prop]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropLeaderConductivity = std::make_shared< fem::Property >();
    tPropLeaderConductivity->set_parameters( { { { 1.0 } } } );
    tPropLeaderConductivity->set_val_function( tGeoValFunction_UTIWGDIFFDIR );

    std::shared_ptr< fem::Property > tPropLeaderDirichlet = std::make_shared< fem::Property >();
    tPropLeaderDirichlet->set_parameters( { { { 1.0 } } } );
    tPropLeaderDirichlet->set_val_function( tGeoValFunction_UTIWGDIFFDIR );

    std::shared_ptr< fem::Property > tPropLeaderNemannPen = std::make_shared< fem::Property >();
    tPropLeaderNemannPen->set_parameters( { { { 0.01 } } } );
    tPropLeaderNemannPen->set_val_function( tConstValFunc_Diff );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderDiffLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMLeaderDiffLinIso->set_dof_type_list( { { MSI::Dof_Type::TEMP } } );
    tCMLeaderDiffLinIso->set_property( tPropLeaderConductivity, "Conductivity" );
    tCMLeaderDiffLinIso->set_space_dim( 3 );
    tCMLeaderDiffLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory                                 tSPFactory;
    std::shared_ptr< fem::Stabilization_Parameter > tSPRobinNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::ROBIN_NITSCHE );
    tSPRobinNitsche->set_parameters( { { { 10.0 } } } );
    tSPRobinNitsche->set_property( tPropLeaderNemannPen, "NeumannPenalty", mtk::Leader_Follower::LEADER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPRobinNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_ROBIN_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( { { MSI::Dof_Type::TEMP } } );
    tIWG->set_dof_type_list( { { MSI::Dof_Type::TEMP } }, mtk::Leader_Follower::LEADER );
    tIWG->set_stabilization_parameter( tSPRobinNitsche, "RobinNitsche" );
    tIWG->set_constitutive_model( tCMLeaderDiffLinIso, "Diffusion" );
    tIWG->set_property( tPropLeaderDirichlet, "Dirichlet" );
    tIWG->set_property( tPropLeaderNemannPen, "NeumannPenalty" );

    // set the normal
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tNormal = { { 1.0 }, { 0.0 }, { 0.0 } };
    tIWG->set_normal( tNormal );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = { { 0.35 }, { -0.25 }, { 0.75 }, { 0.0 } };

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = { { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 } };

    // create time coeff tHat
    Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI.set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule( mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::CONSTANT,
            mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 1 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Vector< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    // set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::SIDESET );

    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    tIWG->mSet->mResDofAssemblyMap.resize( 1 );
    tIWG->mSet->mJacDofAssemblyMap.resize( 1 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 7 } };
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 7 } };

    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 8, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 8, 8, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    tIWG->mRequestedLeaderGlobalDofTypes = { { MSI::Dof_Type::TEMP } };

    Vector< Vector< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager                       tFIManager( tDummy, tSet );

    // set interpolators to the manager
    tFIManager.mFI                     = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager( &tFIManager );

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    tIWG->compute_residual( 1.0 );

    // check evaluation of the jacobian  by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tJacobians;
    Matrix< DDRMat > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
            tEpsilon,
            1.0,
            tJacobians,
            tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

} /* END_TEST_CASE */
