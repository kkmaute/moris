
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

TEST_CASE( "CM_Diff_Lin_Iso_Turbulence", "[CM_Diff_Lin_Iso_Turbulence]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-5;

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
    moris::Cell< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };
    moris::Cell< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tTempDofTypes, tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropTurbPrandtl = std::make_shared< fem::Property >();
    tPropTurbPrandtl->set_val_function( tConstValFunc );
    tPropTurbPrandtl->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_val_function( tConstValFunc );
    tPropConductivity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_val_function( tConstValFunc );
    tPropHeatCapacity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbDiff =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO_TURBULENCE );
    tCMMasterTurbDiff->set_dof_type_list( { tTempDofTypes, tVisDofTypes } );
    tCMMasterTurbDiff->set_property( tPropKinViscosity, "KinematicViscosity" );
    tCMMasterTurbDiff->set_property( tPropTurbPrandtl,  "TurbulentPrandtl" );
    tCMMasterTurbDiff->set_property( tPropDensity,      "Density" );
    tCMMasterTurbDiff->set_property( tPropConductivity, "Conductivity" );
    tCMMasterTurbDiff->set_property( tPropHeatCapacity, "HeatCapacity" );

    tCMMasterTurbDiff->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMMasterTurbDiff->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMMasterTurbDiff->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMMasterTurbDiff->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterTurbDiff->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 0;
    tCMMasterTurbDiff->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set master dof type map
    tCMMasterTurbDiff->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterTurbDiff->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )        = 0;
    tCMMasterTurbDiff->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) )         = 1;

    // build global dof type list
    tCMMasterTurbDiff->get_global_dof_type_list();

    // standard diffusion for comparison
    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiff =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiff->set_dof_type_list( { tTempDofTypes } );
    tCMMasterDiff->set_property( tPropDensity,      "Density" );
    tCMMasterDiff->set_property( tPropConductivity, "Conductivity" );
    tCMMasterDiff->set_property( tPropHeatCapacity, "HeatCapacity" );

    tCMMasterDiff->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set * tSet2 = new fem::Set();
    tCMMasterDiff->set_set_pointer( static_cast< fem::Set* >( tSet2 ) );

    // set size for the set EqnObjDofTypeList
    tCMMasterDiff->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMMasterDiff->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterDiff->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )      = 0;

    // set size and populate the set master dof type map
    tCMMasterDiff->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterDiff->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )        = 0;

    // build global dof type list
    tCMMasterDiff->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // create the jump
        Matrix< DDRMat > tJump( 1, 1, 10.0 );

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

                // set viscosity space derivative parameters
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{ 1.0 },{2.0}} } );
                tPropConductivity->set_parameters( { {{ 2.0 }}, {{ 3.0 },{4.0}} } );
                tPropTurbPrandtl->set_parameters( {  {{ 0.7 }}, {{ 5.0 },{6.0}} } );
                tPropHeatCapacity->set_parameters( { {{ 3.0 }}, {{ 7.0 },{8.0}} } );

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

                // set viscosity space derivative parameters
                tPropKinViscosity->set_parameters( { {{ 1.0 }}, {{1.0},{2.0},{3.0}} } );
                tPropConductivity->set_parameters( { {{ 2.0 }}, {{4.0},{5.0},{6.0}} } );
                tPropTurbPrandtl->set_parameters(  { {{ 0.7 }}, {{7.0},{8.0},{9.0}} } );
                tPropHeatCapacity->set_parameters( { {{ 3.0 }}, {{10.0},{11.0},{12.0}} } );

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
        tCMMasterTurbDiff->set_space_dim( iSpaceDim );
        tCMMasterDiff->set_space_dim( iSpaceDim );

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
            Matrix< DDRMat > tMasterDOFHatTemp;
            fill_phat( tMasterDOFHatTemp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatViscosity;
            fill_phat( tMasterDOFHatViscosity, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatTemp );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatViscosity );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMMasterTurbDiff->mSet->mMasterFIManager = &tFIManager;
            tCMMasterDiff->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMMasterTurbDiff->set_field_interpolator_manager( &tFIManager );
            tCMMasterDiff->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMMasterTurbDiff->reset_eval_flags();
                tCMMasterDiff->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMMasterTurbDiff->mSet->mMasterFIManager->set_space_time( tParamPoint );
                tCMMasterDiff->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes =
                        tCMMasterTurbDiff->get_global_dof_type_list();
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes2 =
                                       tCMMasterDiff->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                        tCMMasterTurbDiff->get_dof_type_list();
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes2 =
                                       tCMMasterDiff->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    Matrix< DDRMat > tflux  = tCMMasterTurbDiff->flux();
                    Matrix< DDRMat > tflux2 = tCMMasterDiff->flux();
                    //print(tflux,"tflux");
                    //print(tflux2,"tflux2");

                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMMasterTurbDiff->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMMasterTurbDiff->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxFluid = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxFluid );

                    // strain
                    //------------------------------------------------------------------------------
                    Matrix< DDRMat > tstrain  = tCMMasterTurbDiff->strain();
                    Matrix< DDRMat > tstrain2 = tCMMasterDiff->strain();
                    //print(tstrain,"tstrain");
                    //print(tstrain2,"tstrain2");

                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMMasterTurbDiff->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMMasterTurbDiff->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrain = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrain );

                    // traction
                    //------------------------------------------------------------------------------
                    Matrix< DDRMat > ttraction  = tCMMasterTurbDiff->traction( tNormal );
                    Matrix< DDRMat > ttraction2 = tCMMasterDiff->traction( tNormal );
                    //print(ttraction,"ttraction");
                    //print(ttraction2,"ttraction2");

                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMMasterTurbDiff->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMMasterTurbDiff->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTraction = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTraction );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tMasterDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tMasterDofTypes( iTestDof );

                        if( tDofTest( 0 ) == MSI::Dof_Type::TEMP )
                        {
                            Matrix< DDRMat > ttesttraction  = tCMMasterTurbDiff->testTraction( tNormal, tDofTest );
                            Matrix< DDRMat > ttesttraction2 = tCMMasterDiff->testTraction( tNormal, tDofTest );
                            //print(ttesttraction,"ttesttraction");
                            //print(ttesttraction2,"ttesttraction2");
                        }

//                        // evaluate dtesttractiondu
//                        Matrix< DDRMat > tdtesttractiondu = tCMMasterTurbDiff->dTestTractiondDOF(
//                                tDofDerivative,
//                                tNormal,
//                                tJump,
//                                tDofTest );
//
//                        // evaluate dtractiondu by FD
//                        Matrix< DDRMat > tdtesttractionduFD;
//                        tCMMasterTurbDiff->eval_dtesttractiondu_FD(
//                                tDofDerivative,
//                                tDofTest,
//                                tdtesttractionduFD,
//                                tPerturbation,
//                                tNormal,
//                                tJump );
//
//                        // check that analytical and FD match
//                        bool tCheckTestTraction = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
//                        REQUIRE( tCheckTestTraction );
                    }

                    // test traction 2
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tMasterDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tMasterDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu2 = tCMMasterTurbDiff->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD2;
                        tCMMasterTurbDiff->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD2,
                                tPerturbation,
                                tNormal);

                        // check that analytical and FD match
                        bool tCheckTestTraction2 = fem::check( tdtesttractiondu2, tdtesttractionduFD2, tEpsilon );
                        REQUIRE( tCheckTestTraction2 );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    Matrix< DDRMat > tdivflux  = tCMMasterTurbDiff->divflux();
                    Matrix< DDRMat > tdivflux2 = tCMMasterDiff->divflux();
                    //print(tdivflux,"tdivflux");
                    //print(tdivflux2,"tdivflux2");

                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMMasterTurbDiff->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMMasterTurbDiff->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFlux = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFlux );

                    // div strain
                    //------------------------------------------------------------------------------
                    Matrix< DDRMat > tdivstrain  = tCMMasterTurbDiff->divstrain();
                    Matrix< DDRMat > tdivstrain2 = tCMMasterDiff->divstrain();
                    //print(tdivstrain,"tdivstrain");
                    //print(tdivstrain2,"tdivstrain2");

                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMMasterTurbDiff->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMMasterTurbDiff->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrain = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrain );
                }
            }
            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/
