
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
#include "cl_FEM_Integrator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "fn_FEM_Check.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

TEST_CASE( "CM_Fluid", "[CM_Fluid]" )
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
    moris::Cell< fem::Integration_Order > tIntegrationOrders = {
            fem::Integration_Order::QUAD_2x2,
            fem::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tVelDofTypes  = { MSI::Dof_Type::VX };
    moris::Cell< MSI::Dof_Type > tPDofTypes    = { MSI::Dof_Type::P };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes, tPDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMMasterFluid->set_dof_type_list( { tVelDofTypes, tPDofTypes } );
    tCMMasterFluid->set_property( tPropViscosity, "Viscosity" );
    tCMMasterFluid->set_property( tPropDensity, "Density" );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMMasterFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMMasterFluid->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMMasterFluid->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // set size and populate the set master dof type map
    tCMMasterFluid->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterFluid->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterFluid->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )         = 1;

    // build global dof type list
    tCMMasterFluid->get_global_dof_type_list();

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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

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
        Interpolation_Rule tGIRule(
                tGeometryType,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimensions for property, CM and SP
        tCMMasterFluid->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            fem::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            fem::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    Integration_Type::GAUSS,
                    fem::Integration_Order::BAR_2 );

            // create an integrator
            fem::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofP   = tNumCoeff;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            Interpolation_Rule tFIRule (
                    tGeometryType,
                    Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMMasterFluid->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMMasterFluid->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMMasterFluid->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMMasterFluid->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes =
                        tCMMasterFluid->get_global_dof_type_list();

                // populate the test master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                        tCMMasterFluid->get_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMMasterFluid->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMMasterFluid->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxFluid = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxFluid );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMMasterFluid->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMMasterFluid->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrainFluid = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrainFluid );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMMasterFluid->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMMasterFluid->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

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
                        Matrix< DDRMat > tdtesttractiondu = tCMMasterFluid->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMMasterFluid->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        bool tCheckTestTractionFluid = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionFluid );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMMasterFluid->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMMasterFluid->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxFluid = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxFluid );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMMasterFluid->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMMasterFluid->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrainFluid = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrainFluid );
                }
            }
            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/


TEST_CASE( "CM_Fluid_Turbulence", "[CM_Fluid_Turbulence]" )
{
    // define an epsilon environment
    // FIXME
    real tEpsilon = 2E-5;

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
    moris::Cell< fem::Integration_Order > tIntegrationOrders = {
            fem::Integration_Order::QUAD_2x2,
            fem::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = {{ 8, 18, 32 },{ 16, 54, 128 }};

    // dof type list
    moris::Cell< MSI::Dof_Type > tVelDofTypes  = { MSI::Dof_Type::VX };
    moris::Cell< MSI::Dof_Type > tVisDofTypes  = { MSI::Dof_Type::VISCOSITY };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = { tVelDofTypes, tVisDofTypes };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunc );
    //tPropViscosity->set_dof_type_list( { tPDofTypes } );
    //tPropViscosity->set_val_function( tPFIValFunc );
    //tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //tPropDensity->set_dof_type_list( { tPDofTypes } );
    //tPropDensity->set_val_function( tPFIValFunc );
    //tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterTurbulence =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_TURBULENCE );
    tCMMasterTurbulence->set_dof_type_list( { tVelDofTypes, tVisDofTypes } );
    tCMMasterTurbulence->set_property( tPropViscosity, "Viscosity" );
    tCMMasterTurbulence->set_property( tPropDensity, "Density" );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tCMMasterTurbulence->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tCMMasterTurbulence->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMMasterTurbulence->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterTurbulence->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // set size and populate the set master dof type map
    tCMMasterTurbulence->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMMasterTurbulence->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )        = 0;
    tCMMasterTurbulence->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VISCOSITY ) ) = 1;

    // build global dof type list
    tCMMasterTurbulence->get_global_dof_type_list();

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create normal
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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

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
                tVelDofTypes = { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ };

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
        Interpolation_Rule tGIRule( tGeometryType,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        // create a space time geometry interpolator
        Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

        // create time coeff tHat
        Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

        // set the coefficients xHat, tHat
        tGI.set_coeff( tXHat, tTHat );

        // set space dimensions for property, CM and SP
        tCMMasterTurbulence->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for( uint iInterpOrder = 1; iInterpOrder < 4; iInterpOrder++ )
        {
            // integration points
            //------------------------------------------------------------------------------
            // get an integration order
            fem::Integration_Order tIntegrationOrder = tIntegrationOrders( iSpaceDim - 2 );

            // create an integration rule
            fem::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    Integration_Type::GAUSS,
                    tIntegrationOrder,
                    mtk::Geometry_Type::LINE,
                    Integration_Type::GAUSS,
                    fem::Integration_Order::BAR_2 );

            // create an integrator
            fem::Integrator tIntegrator( tIntegrationRule );

            // get integration points
            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // number of dof for interpolation order
            uint tNumCoeff = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );

            // get number of dof
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofVis = tNumCoeff;

            //create a space time interpolation rule
            Interpolation_Rule tFIRule (
                    tGeometryType,
                    Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill random coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatVis;
            fill_phat( tMasterDOFHatVis, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator viscosity
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tVisDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatVis );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tCMMasterTurbulence->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tCMMasterTurbulence->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset IWG evaluation flags
                tCMMasterTurbulence->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tCMMasterTurbulence->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // populate the requested master dof type for CM
                moris::Cell< moris::Cell< MSI::Dof_Type > > tRequestedMasterGlobalDofTypes =
                        tCMMasterTurbulence->get_global_dof_type_list();

                // populate the test master dof type for CM
                //                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                //                        tCMMasterTurbulence->get_dof_type_list();
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes = { tVelDofTypes };

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tRequestedMasterGlobalDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tRequestedMasterGlobalDofTypes( iRequestedDof );

                    // flux
                    //------------------------------------------------------------------------------
                    // evaluate dfluxdu
                    Matrix< DDRMat > tdfluxdu = tCMMasterTurbulence->dFluxdDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdfluxduFD;
                    tCMMasterTurbulence->eval_dFluxdDOF_FD( tDofDerivative, tdfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckFluxTurbulence = fem::check( tdfluxdu, tdfluxduFD, tEpsilon );
                    REQUIRE( tCheckFluxTurbulence );

                    // strain
                    //------------------------------------------------------------------------------
                    // evaluate dstraindu
                    Matrix< DDRMat > tdstraindu = tCMMasterTurbulence->dStraindDOF( tDofDerivative );

                    // evaluate dstraindu by FD
                    Matrix< DDRMat > tdstrainduFD;
                    tCMMasterTurbulence->eval_dStraindDOF_FD( tDofDerivative, tdstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckStrainTurbulence = fem::check( tdstraindu, tdstrainduFD, tEpsilon );
                    REQUIRE( tCheckStrainTurbulence );

                    // traction
                    //------------------------------------------------------------------------------
                    // evaluate dtractiondu
                    Matrix< DDRMat > tdtractiondu = tCMMasterTurbulence->dTractiondDOF( tDofDerivative, tNormal );

                    // evaluate dtractiondu by FD
                    Matrix< DDRMat > tdtractionduFD;
                    tCMMasterTurbulence->eval_dtractiondu_FD( tDofDerivative, tdtractionduFD, tPerturbation, tNormal );

                    // check that analytical and FD match
                    bool tCheckTractionTurbulence = fem::check( tdtractiondu, tdtractionduFD, tEpsilon );
                    REQUIRE( tCheckTractionTurbulence );

                    // test traction
                    //------------------------------------------------------------------------------
                    // loop over test dof type
                    for( uint iTestDof = 0; iTestDof < tMasterDofTypes.size(); iTestDof++ )
                    {
                        // get the test dof type
                        Cell< MSI::Dof_Type > tDofTest = tMasterDofTypes( iTestDof );

                        // evaluate dtesttractiondu
                        Matrix< DDRMat > tdtesttractiondu = tCMMasterTurbulence->dTestTractiondDOF(
                                tDofDerivative,
                                tNormal,
                                tJump,
                                tDofTest );

                        // evaluate dtractiondu by FD
                        Matrix< DDRMat > tdtesttractionduFD;
                        tCMMasterTurbulence->eval_dtesttractiondu_FD(
                                tDofDerivative,
                                tDofTest,
                                tdtesttractionduFD,
                                tPerturbation,
                                tNormal,
                                tJump );

                        // check that analytical and FD match
                        bool tCheckTestTractionTurbulence = fem::check( tdtesttractiondu, tdtesttractionduFD, tEpsilon );
                        REQUIRE( tCheckTestTractionTurbulence );
                    }

                    // div flux
                    //------------------------------------------------------------------------------
                    // evaluate ddivfluxdu
                    Matrix< DDRMat > tddivfluxdu = tCMMasterTurbulence->ddivfluxdu( tDofDerivative );

                    // evaluate ddivfluxdu by FD
                    Matrix< DDRMat > tddivfluxduFD;
                    tCMMasterTurbulence->eval_ddivfluxdu_FD( tDofDerivative, tddivfluxduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivFluxTurbulence = fem::check( tddivfluxdu, tddivfluxduFD, tEpsilon );
                    REQUIRE( tCheckDivFluxTurbulence );

                    // div strain
                    //------------------------------------------------------------------------------
                    // evaluate ddivstraindu
                    Matrix< DDRMat > tddivstraindu = tCMMasterTurbulence->ddivstraindu( tDofDerivative );

                    // evaluate ddivstraindu by FD
                    Matrix< DDRMat > tddivstrainduFD;
                    tCMMasterTurbulence->eval_ddivstraindu_FD( tDofDerivative, tddivstrainduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckDivStrainTurbulence = fem::check( tddivstraindu, tddivstrainduFD, tEpsilon );
                    REQUIRE( tCheckDivStrainTurbulence );
                }
            }
            // clean up
            tMasterFIs.clear();
        }
    }
}/*END_TEST_CASE*/


