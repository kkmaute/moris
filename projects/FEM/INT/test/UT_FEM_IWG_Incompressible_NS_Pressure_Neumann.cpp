#include <string>
#include <catch.hpp>
#include <memory>
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
// LINALG/src
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
// FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_MTK_Integrator.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"

using namespace moris;
using namespace fem;

void
UT_FEM_IWG_Incompressible_NS_Pressure_Neumann(
        bool aUsePressure,
        bool aUseTotalPressure,
        bool aUseBackflowPrevention )
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
        mtk::Interpolation_Order::CUBIC
    };

    // create list of integration orders
    moris::Cell< mtk::Integration_Order > tIntegrationOrders = {
        mtk::Integration_Order::QUAD_2x2,
        mtk::Integration_Order::HEX_2x2x2
    };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18, 32 }, { 16, 54, 128 } };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tVelDofTypes = { { MSI::Dof_Type::VX } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tPDofTypes   = { { MSI::Dof_Type::P } };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes    = { tVelDofTypes( 0 ), tPDofTypes( 0 ) };

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { { { 0.3 } } } );
    tPropViscosity->set_val_function( tConstValFunc );
    // tPropViscosity->set_dof_type_list( { tPDofTypes } );
    // tPropViscosity->set_val_function( tPFIValFunc );
    // tPropViscosity->set_dof_derivative_functions( { tPFIDerFunc } );

    std::shared_ptr< fem::Property > tPropKinViscosity = std::make_shared< fem::Property >();
    tPropKinViscosity->set_val_function( tConstValFunc );
    tPropKinViscosity->set_space_der_functions( { tVISCOSITYFISpaceDerFunc } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { { { 0.6 } } } );
    tPropDensity->set_val_function( tConstValFunc );
    // tPropDensity->set_dof_type_list( { tPDofTypes } );
    // tPropDensity->set_val_function( tPFIValFunc );
    // tPropDensity->set_dof_derivative_functions( { tPFIDerFunc } );


    std::shared_ptr< fem::Property > tPropPressure           = nullptr;
    std::shared_ptr< fem::Property > tPropTotalPressure      = nullptr;
    std::shared_ptr< fem::Property > tPropBackflowPrevention = nullptr;

    if ( aUsePressure )
    {
        tPropPressure = std::make_shared< fem::Property >();
        tPropPressure->set_parameters( { { { 3.33 } } } );
        tPropPressure->set_val_function( tConstValFunc );
    }

    if ( aUseTotalPressure )
    {
        tPropTotalPressure = std::make_shared< fem::Property >();
        tPropTotalPressure->set_parameters( { { { 3.33 } } } );
        tPropTotalPressure->set_val_function( tConstValFunc );
    }

    if ( aUseBackflowPrevention )
    {
        tPropBackflowPrevention = std::make_shared< fem::Property >();
        tPropBackflowPrevention->set_parameters( { { { 1.0 } } } );
        tPropBackflowPrevention->set_val_function( tConstValFunc );
    }

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWGNeumann =
            tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_IMPOSED_PRESSURE );
    tIWGNeumann->set_residual_dof_type( tVelDofTypes );
    tIWGNeumann->set_dof_type_list( { tVelDofTypes( 0 ), tPDofTypes( 0 ) }, mtk::Master_Slave::MASTER );
    tIWGNeumann->set_property( tPropDensity, "Density" );
    tIWGNeumann->set_property( tPropPressure, "Pressure" );
    tIWGNeumann->set_property( tPropTotalPressure, "TotalPressure" );
    tIWGNeumann->set_property( tPropBackflowPrevention, "BackFlowPrevention" );

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::SIDESET );
    tIWGNeumann->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWGNeumann->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWGNeumann->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGNeumann->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGNeumann->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

    // set size and populate the set master dof type map
    tIWGNeumann->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGNeumann->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGNeumann->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

    // build global dof type list
    tIWGNeumann->get_global_dof_type_list();

    // loop on the space dimension
    for ( uint iSpaceDim = 2; iSpaceDim < 4; iSpaceDim++ )
    {
        // create the normal
        Matrix< DDRMat > tNormal( iSpaceDim, 1, 0.5 );
        tNormal = tNormal / norm( tNormal );

        // velocity
        Matrix< DDRMat > tVelocity( iSpaceDim, 1, 10.0 );

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

                // set viscosity property parameters
                tPropKinViscosity->set_parameters( { { { 1.0 } }, { { 0.0 }, { 0.0 } } } );
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

                // set viscosity property parameters
                tPropKinViscosity->set_parameters( { { { 1.0 } }, { { 0.0 }, { 0.0 }, { 0.0 } } } );

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

        // set the normal
        tIWGNeumann->set_normal( tNormal );

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

            // get number of dof
            int tNumDofVel = tNumCoeff * iSpaceDim;
            int tNumDofP   = tNumCoeff;

            // create a space time interpolation rule
            mtk::Interpolation_Rule tFIRule(
                    tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill random coefficients for master FI
            Matrix< DDRMat > tMasterDOFHatVel;
            fill_uhat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );

            Matrix< DDRMat > tMasterDOFHatP;
            fill_phat( tMasterDOFHatP, iSpaceDim, iInterpOrder );

            // create a cell of field interpolators for IWG
            Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

            // create the field interpolator velocity ( use negative of default values to trigger backflow )
            tMasterFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelDofTypes( 0 ) );
            tMasterFIs( 0 )->set_coeff( -1.0 * tMasterDOFHatVel );

            // create the field interpolator pressure
            tMasterFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPDofTypes( 0 ) );
            tMasterFIs( 1 )->set_coeff( tMasterDOFHatP );

            // set size and fill the set residual assembly map
            tIWGNeumann->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWGNeumann->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofVel - 1 } };
            tIWGNeumann->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofVel, tNumDofVel + tNumDofP - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofVel - 1 },
                { tNumDofVel, tNumDofVel + tNumDofP - 1 }
            };

            tIWGNeumann->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWGNeumann->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWGNeumann->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWGNeumann->mSet->mResidual.resize( 1 );
            tIWGNeumann->mSet->mResidual( 0 ).set_size( tNumDofVel + tNumDofP, 1, 0.0 );
            tIWGNeumann->mSet->mJacobian.set_size( tNumDofVel + tNumDofP, tNumDofVel + tNumDofP, 0.0 );

            // populate the requested master dof type
            tIWGNeumann->mRequestedMasterGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            moris::Cell< moris::Cell< enum PDV_Type > >        tDummyDv;
            moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager                         tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tMasterFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set
            tIWGNeumann->mSet->mMasterFIManager = &tFIManager;

            // set IWG field interpolator manager
            tIWGNeumann->set_field_interpolator_manager( &tFIManager );

            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWGNeumann->reset_eval_flags();

                // create evaluation point xi, tau
                Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );

                // set integration point
                tIWGNeumann->mSet->mMasterFIManager->set_space_time( tParamPoint );

                // check evaluation of the residual for IWG
                //------------------------------------------------------------------------------
                // reset residual
                tIWGNeumann->mSet->mResidual( 0 ).fill( 0.0 );

                // compute residual
                tIWGNeumann->compute_residual( 1.0 );

                // check evaluation of the jacobian by FD
                //------------------------------------------------------------------------------
                // reset jacobian
                tIWGNeumann->mSet->mJacobian.fill( 0.0 );

                // init the jacobian for IWG and FD evaluation
                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                // check jacobian by FD
                bool tCheckJacobian = tIWGNeumann->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                // check that jacobian is symmetric for symmetric version using upper diagonal block
                uint             tNumRowBlock = tJacobian.n_rows() / iSpaceDim;
                Matrix< DDRMat > tBlock       = tJacobian( { 0, tNumRowBlock - 1 }, { 0, tNumRowBlock - 1 } );

                real tRelError = norm( tBlock - trans( tBlock ) ) / ( norm( tBlock ) + MORIS_REAL_EPS );
                REQUIRE( tRelError < 1e-12 );

                // print for debug
                if ( !tCheckJacobian )
                {
                    std::cout << "Case: Geometry " << iSpaceDim << " Order " << iInterpOrder << " iGP " << iGP << std::endl;
                }

                // require check is true
                REQUIRE( tCheckJacobian );
            }

            // clean up
            tMasterFIs.clear();
        }
    }
}


TEST_CASE( "IWG_Incompressible_NS_Pressure_Neumann",
        "[IWG_Incompressible_NS_Pressure_Neumann]" )
{
    UT_FEM_IWG_Incompressible_NS_Pressure_Neumann( true, false, false );
}


TEST_CASE( "IWG_Incompressible_NS_Pressure_Neumann_TotalPressure",
        "[IWG_Incompressible_NS_Pressure_Neumann_TotalPressure]" )
{
    UT_FEM_IWG_Incompressible_NS_Pressure_Neumann( false, true, false );
}

TEST_CASE( "IWG_Incompressible_NS_Pressure_Neumann_BackflowPrevention",
        "[IWG_Incompressible_NS_Pressure_Neumann_BackflowPrevention]" )
{
    UT_FEM_IWG_Incompressible_NS_Pressure_Neumann( false, false, true );
}
/*END_TEST_CASE*/
