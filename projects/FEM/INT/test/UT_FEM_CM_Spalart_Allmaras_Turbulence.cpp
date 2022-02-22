
#include "catch.hpp"

#define protected public
#define private   public
//FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

//LINALG/src
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
//FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_MTK_Integrator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Spalart_Allmaras_Turbulence", "[CM_Spalart_Allmaras_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-3;

    // define a perturbation relative size
    real tPerturbation = 1.0E-6;

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
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tVisDofTypes  = { { MSI::Dof_Type::VISCOSITY } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes( 0 ), tVisDofTypes( 0 ) };

    // create viscosity property
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // create wall distance property
    std::shared_ptr< fem::Property > tPropWallDistance = std::make_shared< fem::Property >();
    tPropWallDistance->set_parameters( { {{ 1.0 }} } );
    tPropWallDistance->set_val_function( tConstValFunc );

    // define constitutive model
    fem::CM_Factory tCMFactory;

    // CM for Spalart-Allmaras turbulence
    std::shared_ptr< fem::Constitutive_Model > tCMMasterSATurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );
    tCMMasterSATurbulence->set_dof_type_list( tDofTypes );
    tCMMasterSATurbulence->set_property( tPropWallDistance, "WallDistance" );
    tCMMasterSATurbulence->set_property( tPropViscosity,    "KinViscosity" );
    tCMMasterSATurbulence->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMMasterSATurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMMasterSATurbulence->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMMasterSATurbulence->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterSATurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterSATurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set master dof type map
    tCMMasterSATurbulence->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterSATurbulence->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterSATurbulence->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // build global dof type list
    tCMMasterSATurbulence->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( iSpaceDim, 1, 10.0 );

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
        mtk::Interpolation_Rule tGIRule(
                tGeometryType,
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

        // set space dimensions for property, CM and SP
        tCMMasterSATurbulence->set_space_dim( iSpaceDim );

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
            mtk::Interpolation_Rule tFIRule (
                    tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatPVis;
            fill_phat( tMasterDOFHatPVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator viscosity
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatPVis );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMMasterSATurbulence->mSet->mMasterFIManager = &tFIManager;

            // set the field interpolator manager
            tCMMasterSATurbulence->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMMasterSATurbulence->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMMasterSATurbulence->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes =
                        tCMMasterSATurbulence->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                        tCMMasterSATurbulence->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                    // production coeff
                    //------------------------------------------------------------------------------
                    // evaluate dproductioncoeffdu
                    Matrix< DDRMat > tdproductioncoeffdu = tCMMasterSATurbulence->dproductioncoeffdu( tDofDerivative );


                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMMasterSATurbulence->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMMasterSATurbulence->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTractionFluid = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTractionFluid );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tMasterDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tMasterDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMMasterSATurbulence->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMMasterSATurbulence->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal);

                        // check that analytical and FD match
                        bool tCheckTestTractionFluid = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionFluid );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMMasterSATurbulence->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMMasterSATurbulence->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxFluid = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxFluid );

                }
            }
            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/

