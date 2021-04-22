#include <string>
#include <catch.hpp>
#include "assert.hpp"

#define protected public
#define private   public
 //FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private

//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"
//LINALG
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"

using namespace moris;
using namespace fem;

// This UT tests the velocity interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order

TEST_CASE( "IWG_Incompressible_NS_Velocity_Interface_Symmetric_Nitsche",
        "[IWG_Incompressible_NS_Velocity_Interface_Symmetric_Nitsche]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the master properties
    std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
    tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
    tPropMasterViscosity->set_val_function( tConstValFunc );
    //tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropMasterViscosity->set_val_function( tVXFIValFunc );
    //tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
    tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
    tPropMasterDensity->set_val_function( tConstValFunc );
    //tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropMasterViscosity->set_val_function( tVXFIValFunc );
    //tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the slave properties
    std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
    tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveViscosity->set_val_function( tConstValFunc );
    //tPropSlaveViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropSlaveViscosity->set_val_function( tVXFIValFunc );
    //tPropSlaveViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropSlaveDensity = std::make_shared< fem::Property >();
    tPropSlaveDensity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveDensity->set_val_function( tConstValFunc );
    //tPropSlaveDensity->set_dof_type_list( { tVelDofTypes } );
    //tPropSlaveDensity->set_val_function( tVXFIValFunc );
    //tPropSlaveDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMMasterTurbulence->set_property( tPropMasterViscosity, "Viscosity" );
    tCMMasterTurbulence->set_property( tPropMasterDensity, "Density" );
    tCMMasterTurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMSlaveTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMSlaveTurbulence->set_property( tPropSlaveViscosity, "Viscosity" );
    tCMSlaveTurbulence->set_property( tPropSlaveDensity, "Density" );
    tCMSlaveTurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVelDofTypes );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterTurbulence, "IncompressibleFluid", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveTurbulence, "IncompressibleFluid", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) )       = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )        = 1;
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
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

        // set space dimension to CM, SP
        tCMMasterTurbulence->set_space_dim( iSpaceDim );
        tCMSlaveTurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

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
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatVis;
            fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // create the field interpolator viscosity
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );

            // fill coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatVel;
            fill_uhat( tSlaveDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatP;
            fill_phat( tSlaveDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatVis;
            fill_phat( tSlaveDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVel );

            // create the field interpolator pressure
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tSlaveFIs( 1 )->set_coeff( tSlaveDOFHatP );

            // create the field interpolator pressure
            tSlaveFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tSlaveFIs( 2 )->set_coeff( tSlaveDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tMasterFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tSlaveFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Velocity_Interface_Unsymmetric_Nitsche",
        "[IWG_Incompressible_NS_Velocity_Interface_Unsymmetric_Nitsche]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the master properties
    std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
    tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
    tPropMasterViscosity->set_val_function( tConstValFunc );
    //tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropMasterViscosity->set_val_function( tVXFIValFunc );
    //tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
    tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
    tPropMasterDensity->set_val_function( tConstValFunc );
    //tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropMasterViscosity->set_val_function( tVXFIValFunc );
    //tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the slave properties
    std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
    tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveViscosity->set_val_function( tConstValFunc );
    //tPropSlaveViscosity->set_dof_type_list( { tVelDofTypes } );
    //tPropSlaveViscosity->set_val_function( tVXFIValFunc );
    //tPropSlaveViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropSlaveDensity = std::make_shared< fem::Property >();
    tPropSlaveDensity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveDensity->set_val_function( tConstValFunc );
    //tPropSlaveDensity->set_dof_type_list( { tVelDofTypes } );
    //tPropSlaveDensity->set_val_function( tVXFIValFunc );
    //tPropSlaveDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMMasterTurbulence->set_property( tPropMasterViscosity, "Viscosity" );
    tCMMasterTurbulence->set_property( tPropMasterDensity, "Density" );
    tCMMasterTurbulence->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMSlaveTurbulence->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMSlaveTurbulence->set_property( tPropSlaveViscosity, "Viscosity" );
    tCMSlaveTurbulence->set_property( tPropSlaveDensity, "Density" );
    tCMSlaveTurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tVelDofTypes );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes }, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterTurbulence, "IncompressibleFluid", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveTurbulence, "IncompressibleFluid", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) )       = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )        = 1;
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
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

        // set space dimension to CM, SP
        tCMMasterTurbulence->set_space_dim( iSpaceDim );
        tCMSlaveTurbulence->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

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
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatVis;
            fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // create the field interpolator viscosity
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );

            // fill coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatVel;
            fill_uhat( tSlaveDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatP;
            fill_phat( tSlaveDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatVis;
            fill_phat( tSlaveDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVel );

            // create the field interpolator pressure
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tSlaveFIs( 1 )->set_coeff( tSlaveDOFHatP );

            // create the field interpolator pressure
            tSlaveFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tSlaveFIs( 2 )->set_coeff( tSlaveDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tMasterFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tSlaveFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
}/*END_TEST_CASE*/
// This UT tests the pressure interface IWG for incompressible NS
// for QUAD, HEX geometry type
// for LINEAR, QUADRATIC and CUBIC interpolation order

TEST_CASE( "IWG_Incompressible_NS_Pressure_Interface_Symmetric_Nitsche",
        "[IWG_Incompressible_NS_Pressure_Interface_Symmetric_Nitsche]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the master properties
    std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
    tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
    tPropMasterViscosity->set_val_function( tConstValFunc );
//    tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropMasterViscosity->set_val_function( tVXFIValFunc );
//    tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
    tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
    tPropMasterDensity->set_val_function( tConstValFunc );

//    tPropMasterDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropMasterDensity->set_val_function( tVXFIValFunc );
//    tPropMasterDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the slave properties
    std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
    tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveViscosity->set_val_function( tConstValFunc );
//    tPropSlaveViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropSlaveViscosity->set_val_function( tVXFIValFunc );
//    tPropSlaveViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropSlaveDensity = std::make_shared< fem::Property >();
    tPropSlaveDensity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveDensity->set_val_function( tConstValFunc );
//    tPropSlaveDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropSlaveDensity->set_val_function( tVXFIValFunc );
//    tPropSlaveDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMMasterFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMMasterFluid->set_property( tPropMasterViscosity, "Viscosity" );
    tCMMasterFluid->set_property( tPropMasterDensity, "Density" );
    tCMMasterFluid->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMSlaveFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMSlaveFluid->set_property( tPropSlaveViscosity, "Viscosity" );
    tCMSlaveFluid->set_property( tPropSlaveDensity, "Density" );
    tCMSlaveFluid->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_SYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterFluid, "IncompressibleFluid", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveFluid, "IncompressibleFluid", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
         Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
         tNormal = tNormal / norm( tNormal );
         tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
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

        // set space dimension to CM, SP
        tCMMasterFluid->set_space_dim( iSpaceDim );
        tCMSlaveFluid->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatVis;
            fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // create the field interpolator viscosity
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );

            // fill coefficients for slave FI
            Matrix< DDRMat > tSlaveDOFHatVel;
            fill_uhat( tSlaveDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatP;
            fill_phat( tSlaveDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatVis;
            fill_phat( tSlaveDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVel );

            // create the field interpolator pressure
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tSlaveFIs( 1 )->set_coeff( tSlaveDOFHatP );

            // create the field interpolator pressure
            tSlaveFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tSlaveFIs( 2 )->set_coeff( tSlaveDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tMasterFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tSlaveFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Pressure_Interface_Unsymmetric_Nitsche",
        "[IWG_Incompressible_NS_Pressure_Interface_Unsymmetric_Nitsche]" )
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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes    = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the master properties
    std::shared_ptr< fem::Property > tPropMasterViscosity = std::make_shared< fem::Property >();
    tPropMasterViscosity->set_parameters( { {{ 1.0 }} } );
    tPropMasterViscosity->set_val_function( tConstValFunc );
//    tPropMasterViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropMasterViscosity->set_val_function( tVXFIValFunc );
//    tPropMasterViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropMasterDensity = std::make_shared< fem::Property >();
    tPropMasterDensity->set_parameters( { {{ 1.0 }} } );
    tPropMasterDensity->set_val_function( tConstValFunc );

//    tPropMasterDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropMasterDensity->set_val_function( tVXFIValFunc );
//    tPropMasterDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // create the slave properties
    std::shared_ptr< fem::Property > tPropSlaveViscosity = std::make_shared< fem::Property >();
    tPropSlaveViscosity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveViscosity->set_val_function( tConstValFunc );
//    tPropSlaveViscosity->set_dof_type_list( { tVelDofTypes } );
//    tPropSlaveViscosity->set_val_function( tVXFIValFunc );
//    tPropSlaveViscosity->set_dof_derivative_functions( { tVXFIDerFunc } );

    std::shared_ptr< fem::Property > tPropSlaveDensity = std::make_shared< fem::Property >();
    tPropSlaveDensity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveDensity->set_val_function( tConstValFunc );
//    tPropSlaveDensity->set_dof_type_list( { tVelDofTypes } );
//    tPropSlaveDensity->set_val_function( tVXFIValFunc );
//    tPropSlaveDensity->set_dof_derivative_functions( { tVXFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMMasterFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMMasterFluid->set_property( tPropMasterViscosity, "Viscosity" );
    tCMMasterFluid->set_property( tPropMasterDensity, "Density" );
    tCMMasterFluid->set_local_properties();

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMSlaveFluid->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes } );
    tCMSlaveFluid->set_property( tPropSlaveViscosity, "Viscosity" );
    tCMSlaveFluid->set_property( tPropSlaveDensity, "Density" );
    tCMSlaveFluid->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitsche->set_parameters( { {{ 1.0 }} } );
    tSPNitsche->set_property( tPropMasterViscosity, "Material", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropSlaveViscosity, "Material", mtk::Master_Slave::SLAVE );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_UNSYMMETRIC_NITSCHE );
    tIWG->set_residual_dof_type( tPDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_model( tCMMasterFluid, "IncompressibleFluid", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveFluid, "IncompressibleFluid", mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitscheInterface" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::DOUBLE_SIDESET );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set slave dof type map
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >( MSI::Dof_Type::P ) )  = 1;
    tIWG->mSet->mSlaveDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
         Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
         tNormal = tNormal / norm( tNormal );
         tIWG->set_normal( tNormal );

        // set geometry inputs
        //------------------------------------------------------------------------------
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

        // set space dimension to CM, SP
        tCMMasterFluid->set_space_dim( iSpaceDim );
        tCMSlaveFluid->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );

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
            int tNumDofVel  = tNumCoeff * iSpaceDim;
            int tNumDofP    = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule ( tGeometryType,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         tInterpolationOrder,
                                         mtk::Interpolation_Type::LAGRANGE,
                                         mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatVis;
            fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // create the field interpolator viscosity
            tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 2 )->set_coeff( tMasterDOFHatVis );

            // fill coefficients for master FI
            Matrix< DDRMat > tSlaveDOFHatVel;;
            fill_uhat( tSlaveDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatP;
            fill_phat( tSlaveDOFHatP, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tSlaveDOFHatVis;
            fill_phat( tSlaveDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tSlaveFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tSlaveFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tSlaveFIs( 0 )->set_coeff( tSlaveDOFHatVel );

            // create the field interpolator pressure
            tSlaveFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tSlaveFIs( 1 )->set_coeff( tSlaveDOFHatP );

            // create the field interpolator pressure
            tSlaveFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tSlaveFIs( 2 )->set_coeff( tSlaveDOFHatVis );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofP + tNumDofVel - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 3 ) = { { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 4 ) = { { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 5 ) = { { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                    { 0, tNumDofVel-1 },
                    { tNumDofVel, tNumDofVel + tNumDofP - 1 },
                    { tNumDofP + tNumDofVel, tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { tNumDofVel + tNumDofP + tNumDofVis, 2 * tNumDofVel + tNumDofP + tNumDofVis - 1 },
                    { 2 * tNumDofVel + tNumDofP + tNumDofVis, 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis - 1 },
                    { 2 * ( tNumDofVel + tNumDofP ) + tNumDofVis, 2 * ( tNumDofVel + tNumDofP + tNumDofVis ) - 1 } };
            tIWG->mSet->mJacDofAssemblyMap.resize( 2 * tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 3 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 4 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 5 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    1,
                    0.0 );
            tIWG->mSet->mJacobian.set_size(
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    2 * ( tNumDofVel + tNumDofP + tNumDofVis ),
                    0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested master dof type
            tIWG->mRequestedMasterGlobalDofTypes = tDofTypes;
            tIWG->mRequestedSlaveGlobalDofTypes  = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tMasterFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
            Field_Interpolator_Manager tSlaveFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tMasterFIManager.mFI = tMasterFIs;
            tMasterFIManager.mIPGeometryInterpolator = &tGI;
            tMasterFIManager.mIGGeometryInterpolator = &tGI;
            tSlaveFIManager.mFI = tSlaveFIs;
            tSlaveFIManager.mIPGeometryInterpolator = &tGI;
            tSlaveFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWG->mSet->mMasterFIManager = &tMasterFIManager;
            tIWG->mSet->mSlaveFIManager = &tSlaveFIManager;

            // set IWG field interpolator manager
            tIWG->set_field_interpolator_manager( &tMasterFIManager, mtk::Master_Slave::MASTER );
            tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

            // loop iver integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
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

                // print for debug
                if( !tCheckJacobian )
                {
                    std::cout<<"Case: Geometry "<<iSpaceDim<<" Order "<<iInterpOrder<<"iGP "<<iGP<<std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
            tSlaveFIs.clear();
        }
    }
}/*END_TEST_CASE*/
