
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

void UT_FEM_SP_YZBETA_Advection_Core( real aBeta )
{
    // define an epsilon environment
    real tEpsilon = 1E-7;

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
    moris::Cell< MSI::Dof_Type > tTempDofTypes = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes  = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes     = { tVelDofTypes( 0 ), tTempDofTypes };

    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );
    //    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropHeatCapacity->set_val_function( tTEMPFIValFunc );
    //    tPropHeatCapacity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_val_function( tConstValFunc );
    //    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropDensity->set_val_function( tTEMPFIValFunc );
    //    tPropDensity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_val_function( tConstValFunc );
    //    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    //    tPropConductivity->set_val_function( tTEMPFIValFunc );
    //    tPropConductivity->set_dof_derivative_functions( { tTEMPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropLoad = std::make_shared< fem::Property >();
    tPropLoad->set_parameters( { {{ 1.5 }} } );
    tPropLoad->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropAdvectionSource = std::make_shared< fem::Property >();
    tPropAdvectionSource->set_parameters( { {{ 2.0 * 1.5 }} } );
    tPropAdvectionSource->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropBeta = std::make_shared< fem::Property >();
    tPropBeta->set_parameters( { {{ aBeta }} } );
    tPropBeta->set_val_function( tConstValFunc );

    std::shared_ptr< fem::Property > tPropRefState = std::make_shared< fem::Property >();
    tPropRefState->set_parameters( { {{ 0.001 }} } );
    tPropRefState->set_val_function( tConstValFunc );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion =
            tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( { tTempDofTypes } );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity,      "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "HeatCapacity" );
    tCMDiffusion->set_local_properties();

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPYZBeta =
            tSPFactory.create_SP( fem::Stabilization_Type::YZBETA_ADVECTION );
    tSPYZBeta->set_dof_type_list( { tTempDofTypes } );
    tSPYZBeta->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tSPYZBeta->set_property( tPropBeta,      "Beta",           mtk::Master_Slave::MASTER );
    tSPYZBeta->set_property( tPropRefState,  "ReferenceState", mtk::Master_Slave::MASTER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPYZBeta->set_cluster( tCluster );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tSPYZBeta->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tCMDiffusion->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tSPYZBeta->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tSPYZBeta->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPYZBeta->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tSPYZBeta->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set master dof type map
    tSPYZBeta->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSPYZBeta->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )   = 0;
    tSPYZBeta->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // loop on the space dimension
    for( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // set space dim
        tSPYZBeta->set_space_dim( iSpaceDim );
        tCMDiffusion->set_space_dim( iSpaceDim );

        // create and set normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );
        tSPYZBeta->set_normal( tNormal );

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

            // fill coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tMasterDOFHatTEMP;
            fill_phat( tMasterDOFHatTEMP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( tMasterDOFHatVel );

            // create the field interpolator temperature
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDofTypes );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatTEMP );

            // build global dof type list
            tSPYZBeta->get_global_dof_type_list();

            // set order
            tSPYZBeta->set_interpolation_order( iInterpOrder );

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tSPYZBeta->mSet->mMasterFIManager = &tFIManager;

            // set SP field interpolator manager
            tSPYZBeta->set_field_interpolator_manager( &tFIManager );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
            {
                // reset SP evaluation flags
                tSPYZBeta->reset_eval_flags();
                tCMDiffusion->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tSPYZBeta->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // populate the requested master dof type for SP
                moris::Cell< moris::Cell< MSI::Dof_Type > > tMasterDofTypes =
                        tSPYZBeta->get_global_dof_type_list();

                // loop over requested dof type
                for( uint iRequestedDof = 0; iRequestedDof < tMasterDofTypes.size(); iRequestedDof++ )
                {
                    // derivative dof type
                    Cell< MSI::Dof_Type > tDofDerivative = tMasterDofTypes( iRequestedDof );

                    // evaluate dspdu
                    Matrix< DDRMat > tdspdu = tSPYZBeta->dSPdMasterDOF( tDofDerivative );

                    // evaluate dfluxdu by FD
                    Matrix< DDRMat > tdspduFD;
                    tSPYZBeta->eval_dSPdMasterDOF_FD( tDofDerivative, tdspduFD, tPerturbation );

                    // check that analytical and FD match
                    bool tCheckYZBeta = fem::check( tdspdu, tdspduFD, tEpsilon );
                    REQUIRE( tCheckYZBeta );
                }
            }

            // clean up
            tMasterFIs.clear();
        }
    }
}

TEST_CASE( "SP_YZBeta_Advection_Beq1", "[SP_YZBeta_Advection_Beq1]" )
{
    UT_FEM_SP_YZBETA_Advection_Core( 1.0 );
}

TEST_CASE( "SP_YZBeta_Advection_Beq2", "[SP_YZBeta_Advection_Beq2]" )
{
    UT_FEM_SP_YZBETA_Advection_Core( 2.0 );
}
