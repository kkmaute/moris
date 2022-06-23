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
// FEM//INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
// LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Struc_Nonlinear_Interface_Unsymmetric_Saint_Venant_Kirchhoff", "[moris],[fem],[IWG_Struc_Nonlinear_Interface_Unsymmetric_Saint_Venant_Kirchhoff]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property >();
    tPropMasterEMod->set_parameters( { { { 10.0 } } } );
    tPropMasterEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( { { { 0.3 } } } );
    tPropMasterNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveEMod = std::make_shared< fem::Property >();
    tPropSlaveEMod->set_parameters( { { { 20.0 } } } );
    tPropSlaveEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveNu = std::make_shared< fem::Property >();
    tPropSlaveNu->set_parameters( { { { 0.3 } } } );
    tPropSlaveNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMMasterStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMSlaveStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
    tSPNitscheInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_SE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveStrucLinIso, "ElastLinIso", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMMasterStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMMasterStrucLinIso->set_space_dim( iSpaceDim );
        tCMSlaveStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMSlaveStrucLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
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

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatDisp;
            fill_uhat_Elast( tMasterDOFHatDisp, iSpaceDim, iInterpOrder );

            tMasterDOFHatDisp = tMasterDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatDisp );

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatDisp;
            fill_uhat_Elast( tSlaveDOFHatDisp, iSpaceDim, iInterpOrder );

            tSlaveDOFHatDisp = tSlaveDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, 2 * tNumDofDisp - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tMasterFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tSlaveFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI                     = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI                      = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator  = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager  = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mSlaveFIManager->set_space_time( tParamPoint );

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


                //                // print for debug
                //                if ( !tCheckJacobian )
                //                {
                //                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                //                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
} /* END_TEST_CASE */


TEST_CASE( "IWG_Struc_Nonlinear_Interface_Symmetric_Saint_Venant_Kirchhoff", "[moris],[fem],[IWG_Struc_Nonlinear_Interface_Symmetric_Saint_Venant_Kirchhoff]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-5;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property >();
    tPropMasterEMod->set_parameters( { { { 10.0 } } } );
    tPropMasterEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( { { { 0.3 } } } );
    tPropMasterNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveEMod = std::make_shared< fem::Property >();
    tPropSlaveEMod->set_parameters( { { { 20.0 } } } );
    tPropSlaveEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveNu = std::make_shared< fem::Property >();
    tPropSlaveNu->set_parameters( { { { 0.3 } } } );
    tPropSlaveNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMMasterStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF );
    tCMSlaveStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
    tSPNitscheInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_SE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveStrucLinIso, "ElastLinIso", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMMasterStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMMasterStrucLinIso->set_space_dim( iSpaceDim );
        tCMSlaveStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMSlaveStrucLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
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

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatDisp;
            fill_uhat_Elast( tMasterDOFHatDisp, iSpaceDim, iInterpOrder );

            tMasterDOFHatDisp = tMasterDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatDisp );

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatDisp;
            fill_uhat_Elast( tSlaveDOFHatDisp, iSpaceDim, iInterpOrder );

            tSlaveDOFHatDisp = tSlaveDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, 2 * tNumDofDisp - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tMasterFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tSlaveFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI                     = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI                      = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator  = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager  = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mSlaveFIManager->set_space_time( tParamPoint );

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


                //                // print for debug
                //                if ( !tCheckJacobian )
                //                {
                //                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                //                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
} /* END_TEST_CASE */


TEST_CASE( "IWG_Struc_Nonlinear_Interface_Unsymmetric_Neo_Hookean", "[moris],[fem],[IWG_Struc_Nonlinear_Interface_Unsymmetric_Neo_Hookean]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property >();
    tPropMasterEMod->set_parameters( { { { 10.0 } } } );
    tPropMasterEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( { { { 0.3 } } } );
    tPropMasterNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveEMod = std::make_shared< fem::Property >();
    tPropSlaveEMod->set_parameters( { { { 20.0 } } } );
    tPropSlaveEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveNu = std::make_shared< fem::Property >();
    tPropSlaveNu->set_parameters( { { { 0.3 } } } );
    tPropSlaveNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN );
    tCMMasterStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN );
    tCMSlaveStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
    tSPNitscheInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_SE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveStrucLinIso, "ElastLinIso", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMMasterStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMMasterStrucLinIso->set_space_dim( iSpaceDim );
        tCMSlaveStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMSlaveStrucLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
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

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatDisp;
            fill_uhat_Elast( tMasterDOFHatDisp, iSpaceDim, iInterpOrder );

            tMasterDOFHatDisp = tMasterDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatDisp );

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatDisp;
            fill_uhat_Elast( tSlaveDOFHatDisp, iSpaceDim, iInterpOrder );

            tSlaveDOFHatDisp = tSlaveDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, 2 * tNumDofDisp - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tMasterFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tSlaveFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI                     = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI                      = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator  = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager  = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mSlaveFIManager->set_space_time( tParamPoint );

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


                //                // print for debug
                //                if ( !tCheckJacobian )
                //                {
                //                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                //                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
} /* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Nonlinear_Interface_Symmetric_Neo_Hookean", "[moris],[fem],[IWG_Struc_Nonlinear_Interface_Symmetric_Neo_Hookean]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 1E-5;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property >();
    tPropMasterEMod->set_parameters( { { { 10.0 } } } );
    tPropMasterEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( { { { 0.3 } } } );
    tPropMasterNu->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveEMod = std::make_shared< fem::Property >();
    tPropSlaveEMod->set_parameters( { { { 20.0 } } } );
    tPropSlaveEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropSlaveNu = std::make_shared< fem::Property >();
    tPropSlaveNu->set_parameters( { { { 0.3 } } } );
    tPropSlaveNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN );
    tCMMasterStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_NON_LIN_ISO_NEO_HOOKEAN );
    tCMSlaveStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { { { 100.0 } } } );
    tSPNitscheInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPNitscheInterface->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_SE );
    tIWG->set_residual_dof_type( tDispDofTypes );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDispDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveStrucLinIso, "ElastLinIso", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set space dim for CM
        tCMMasterStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMMasterStrucLinIso->set_space_dim( iSpaceDim );
        tCMSlaveStrucLinIso->set_model_type( fem::Model_Type::PLANE_STRAIN );
        tCMSlaveStrucLinIso->set_space_dim( iSpaceDim );
        tSPNitscheInterface->set_space_dim( iSpaceDim );

        // set geometry inputs
        //------------------------------------------------------------------------------
        // switch on space dimension
        switch ( iSpaceDim )
        {
            case 2:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::QUAD;
                break;
            }
            case 3:
            {
                // set geometry type
                tGeometryType = mtk::Geometry_Type::HEX;
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

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatDisp;
            fill_uhat_Elast( tMasterDOFHatDisp, iSpaceDim, iInterpOrder );

            tMasterDOFHatDisp = tMasterDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatDisp );

            // fill random coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatDisp;
            fill_uhat_Elast( tSlaveDOFHatDisp, iSpaceDim, iInterpOrder );

            tSlaveDOFHatDisp = tSlaveDOFHatDisp * 0.01;

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDispDofTypes.size() );

            // create the field interpolator displacement
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatDisp );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDispDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, 2 * tNumDofDisp - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, 2 * tNumDofDisp - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( 2 * tNumDofDisp, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( 2 * tNumDofDisp, 2 * tNumDofDisp, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDispDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDispDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tMasterFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager                         tSlaveFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI                     = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI                      = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator  = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator  = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager  = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tIWG->mSet->mSlaveFIManager->set_space_time( tParamPoint );

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


                //                // print for debug
                //                if ( !tCheckJacobian )
                //                {
                //                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                //                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
} /* END_TEST_CASE */
