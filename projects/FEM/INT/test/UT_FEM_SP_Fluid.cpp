
#include "catch.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#undef protected
#undef private

// LINALG/src
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
// FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "SP_Fluid", "[SP_Fluid]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // material parameter
    real tFluidDensity   = 1.0;
    real tFluidViscosity = 1.0;
    real tGammaNitsche   = 1.0;
    real tGammaGPmu      = 1.0;
    real tGammaGPu       = 1.0;
    real tGammaGPp       = 1.0;

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
    moris::Cell< MSI::Dof_Type > tVisDofTypes = { MSI::Dof_Type::VISCOSITY };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes   = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tVelDofTypes( 0 ), tPDofTypes( 0 ), tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropFluidDensity = std::make_shared< fem::Property >();
    tPropFluidDensity->set_parameters( { { { tFluidDensity } } } );
    tPropFluidDensity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropFluidViscosity = std::make_shared< fem::Property >();
    tPropFluidViscosity->set_parameters( { { { tFluidViscosity } } } );
    tPropFluidViscosity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
    tPropWallDistance->set_parameters( { { { 1.0 } } } );
    tPropWallDistance->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMSATurbulence->set_dof_type_list( tDofTypes );
    tCMSATurbulence->set_property( tPropWallDistance, "WallDistance" );
    tCMSATurbulence->set_property( tPropFluidViscosity, "KinViscosity" );
    tCMSATurbulence->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow =
            tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
    tSPIncFlow->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } }, mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_parameters( { { { 36.0 } } } );

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPNitsche->set_parameters( { { { tGammaNitsche } }, { { 1.0 } } } );

    std::shared_ptr< fem::Stabilization_Parameter > tSPViscousGhost =
            tSPFactory.create_SP( fem::Stabilization_Type::VISCOUS_GHOST );
    tSPViscousGhost->set_parameters( { { { tGammaGPmu } } } );
    tSPViscousGhost->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );

    std::shared_ptr< fem::Stabilization_Parameter > tSPConvectiveGhost =
            tSPFactory.create_SP( fem::Stabilization_Type::CONVECTIVE_GHOST );
    tSPConvectiveGhost->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
    tSPConvectiveGhost->set_parameters( { { { tGammaGPu } } } );
    tSPConvectiveGhost->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );

    std::shared_ptr< fem::Stabilization_Parameter > tSPPressureGhost =
            tSPFactory.create_SP( fem::Stabilization_Type::PRESSURE_GHOST );
    tSPPressureGhost->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY } }, mtk::Master_Slave::MASTER );
    tSPPressureGhost->set_parameters( { { { tGammaGPp } }, { { 1.0 } } } );
    tSPPressureGhost->set_property( tPropFluidViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPPressureGhost->set_property( tPropFluidDensity, "Density", mtk::Master_Slave::MASTER );

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPGSA =
            tSPFactory.create_SP( fem::Stabilization_Type::SUPG_SPALART_ALLMARAS_TURBULENCE );
    tSPSUPGSA->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::VISCOSITY } }, mtk::Master_Slave::MASTER );
    tSPSUPGSA->set_constitutive_model( tCMSATurbulence, "SpalartAllmarasTurbulence" );
    tSPSUPGSA->set_parameters( { { { 3.0 } }, { { 1.0 } } } );

    // create a dummy fem cluster and set it to SP
    fem::Cluster* tCluster = new fem::Cluster();
    tSPIncFlow->set_cluster( tCluster );
    tSPNitsche->set_cluster( tCluster );
    tSPViscousGhost->set_cluster( tCluster );
    tSPConvectiveGhost->set_cluster( tCluster );
    tSPPressureGhost->set_cluster( tCluster );
    tSPSUPGSA->set_cluster( tCluster );

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    tSPIncFlow->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tSPNitsche->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tSPViscousGhost->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tSPConvectiveGhost->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tSPPressureGhost->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tSPSUPGSA->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tSPIncFlow->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tSPIncFlow->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPIncFlow->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tSPIncFlow->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tSPIncFlow->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // set size and populate the set master dof type map
    tSPIncFlow->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPIncFlow->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tSPIncFlow->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;
    tSPIncFlow->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 2;

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set space dim
        tSPIncFlow->set_space_dim( iSpaceDim );
        tSPNitsche->set_space_dim( iSpaceDim );
        tSPViscousGhost->set_space_dim( iSpaceDim );
        tSPConvectiveGhost->set_space_dim( iSpaceDim );
        tSPPressureGhost->set_space_dim( iSpaceDim );
        tSPSUPGSA->set_space_dim( iSpaceDim );
        tCMSATurbulence->set_space_dim( iSpaceDim );

        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tSPIncFlow->set_normal( tNormal );
        tSPNitsche->set_normal( tNormal );
        tSPViscousGhost->set_normal( tNormal );
        tSPConvectiveGhost->set_normal( tNormal );
        tSPPressureGhost->set_normal( tNormal );
        tSPSUPGSA->set_normal( tNormal );

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

                // set velocity dof types
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
        Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

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

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // loop over different state variable configurations
            for ( uint iCu = 0; iCu < 4; ++iCu )
            {
                for ( uint iCp = 0; iCp < 3; ++iCp )
                {
                    // fill coefficients for master FI
                    Matrix< DDRMat > tMasterDOFHatVel;
                    fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tMasterDOFHatP;
                    fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
                    Matrix< DDRMat > tMasterDOFHatVis;
                    fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

                    // create a cell of field interpolators
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

                    // build global dof type list
                    tSPIncFlow->get_global_dof_type_list();
                    tSPNitsche->get_global_dof_type_list();
                    tSPViscousGhost->get_global_dof_type_list();
                    tSPConvectiveGhost->get_global_dof_type_list();
                    tSPPressureGhost->get_global_dof_type_list();
                    tSPSUPGSA->get_global_dof_type_list();

                    // set order
                    tSPIncFlow->set_interpolation_order( iInterpOrder );
                    tSPNitsche->set_interpolation_order( iInterpOrder );
                    tSPViscousGhost->set_interpolation_order( iInterpOrder );
                    tSPConvectiveGhost->set_interpolation_order( iInterpOrder );
                    tSPPressureGhost->set_interpolation_order( iInterpOrder );
                    tSPSUPGSA->set_interpolation_order( iInterpOrder );

                    // create a field interpolator manager
                    moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
                    moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
                    Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

                    // populate the field interpolator manager
                    tFIManager.mFI                     = tMasterFIs;
                    tFIManager.mIPGeometryInterpolator = &tGI;
                    tFIManager.mIGGeometryInterpolator = &tGI;

                    // set the interpolator manager to the set
                    tSPIncFlow->mSet->mMasterFIManager = &tFIManager;

                    // set SP field interpolator manager
                    tSPIncFlow->set_field_interpolator_manager( &tFIManager );
                    tSPNitsche->set_field_interpolator_manager( &tFIManager );
                    tSPViscousGhost->set_field_interpolator_manager( &tFIManager );
                    tSPConvectiveGhost->set_field_interpolator_manager( &tFIManager );
                    tSPPressureGhost->set_field_interpolator_manager( &tFIManager );
                    tSPSUPGSA->set_field_interpolator_manager( &tFIManager );

                    // loop over integration points
                    uint tNumGPs = tIntegPoints.n_cols();
                    for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
                    {
                        // reset SP evaluation flags
                        tSPIncFlow->reset_eval_flags();
                        tSPNitsche->reset_eval_flags();
                        tSPViscousGhost->reset_eval_flags();
                        tSPConvectiveGhost->reset_eval_flags();
                        tSPPressureGhost->reset_eval_flags();
                        tSPSUPGSA->reset_eval_flags();

                        // create evaluation point xi, tau
                        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                        // set integration point
                        tSPIncFlow->mSet->mMasterFIManager->set_space_time( tParamPoint );

                        // populate the requested master dof type for SP
                        moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                                tSPIncFlow->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPIncFlow->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPIncFlow->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPSUPGPSPG = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPSUPGPSPG );
                        }

                        // populate the requested master dof type for SP
                        tMasterDofTypes = tSPNitsche->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPNitsche->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPNitsche->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPNitsche = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPNitsche );
                        }

                        // populate the requested master dof type for SP
                        tMasterDofTypes = tSPViscousGhost->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPViscousGhost->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPViscousGhost->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPViscousGhost = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPViscousGhost );
                        }

                        // populate the requested master dof type for SP
                        tMasterDofTypes = tSPConvectiveGhost->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPConvectiveGhost->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPConvectiveGhost->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPConvectiveGhost = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPConvectiveGhost );
                        }

                        // populate the requested master dof type for SP
                        tMasterDofTypes = tSPPressureGhost->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPPressureGhost->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPPressureGhost->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPPressureGhost = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPPressureGhost );
                        }

                        // populate the requested master dof type for SP
                        tMasterDofTypes = tSPSUPGSA->get_global_dof_type_list();

                        // loop over requested dof type
                        for ( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                        {
                            // derivative dof type
                            Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                            // evaluate dspdu
                            Matrix< DDRMat > tdspdu = tSPSUPGSA->dSPdMasterDOF( tDofDerivative );

                            // evaluate dfluxdu by FD
                            Matrix< DDRMat > tdspduFD;
                            tSPSUPGSA->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                            // check that analytical and FD match
                            bool tCheckSPSUPGSA = fem::check( tdspdu, tdspduFD, tEpsilon );
                            REQUIRE( tCheckSPSUPGSA );
                        }
                    }

                    // clean up
                    tMasterFIs.clear();
                }
            }
        }
    }
} /*END_TEST_CASE*/
