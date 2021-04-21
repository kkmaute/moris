#include <string>
#include <catch.hpp>
#include <memory>
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
//LINALG/src
#include "op_equal_equal.hpp"
//MTK/src
#include "cl_MTK_Enums.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_UT.cpp"

#include "fn_trans.hpp"
#include "fn_FEM_Check.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Compressible_NS_Bulk_Perfect_Gas_Pressure_Primitive_Analytical",
        "[IWG_Compressible_NS_Bulk_Perfect_Gas_Pressure_Primitive_Analytical]" )
{
    // define an epsilon environment
    real tEpsilon = 5.0E-3;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // dof type list
    moris::Cell< MSI::Dof_Type > tPressureDof = { MSI::Dof_Type::P };
    moris::Cell< MSI::Dof_Type > tVelocityDof = { MSI::Dof_Type::VX };
    moris::Cell< MSI::Dof_Type > tTempDof     = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes         = { tPressureDof, tVelocityDof, tTempDof };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tResidualDofTypes = { { MSI::Dof_Type::P, MSI::Dof_Type::VX, MSI::Dof_Type::TEMP } };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties

    // dummy factor by which strong form of the residual gets multiplied and added to weak form
    real tDummyFactor = 1.3;

    // dynamic viscosity
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.716e-5 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 0.718e3 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 287.058}} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 24.35e-3 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    // dummy for SP
    std::shared_ptr< fem::Property > tPropDummyForSP = std::make_shared< fem::Property >();
    tPropDummyForSP->set_parameters( { {{ 1.0 }} } );
    tPropDummyForSP->set_val_function( tConstValFunc );

    // define material model and assign properties
    fem::MM_Factory tMMFactory;

    std::shared_ptr< fem::Material_Model > tMMFluid =
            tMMFactory.create_CM( fem::Material_Type::PERFECT_GAS );
    tMMFluid->set_dof_type_list( {tPressureDof, tTempDof } );
    tMMFluid->set_property( tPropHeatCapacity, "IsochoricHeatCapacity" );
    tMMFluid->set_property( tPropGasConstant,  "SpecificGasConstant" );    

    // define constitutive model and assign properties
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterFluid =
            tCMFactory.create_CM( fem::Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN );
    tCMMasterFluid->set_dof_type_list( {tPressureDof, tVelocityDof, tTempDof } );
    tCMMasterFluid->set_property( tPropViscosity,    "DynamicViscosity" );
    tCMMasterFluid->set_property( tPropConductivity, "ThermalConductivity" );
    tCMMasterFluid->set_material_model( tMMFluid, "ThermodynamicMaterialModel" );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;
    std::shared_ptr< fem::Stabilization_Parameter > tSP =
            tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSP->set_parameters( { {{ tDummyFactor }} });
    tSP->set_property( tPropDummyForSP, "Material", mtk::Master_Slave::MASTER );

    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSP->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_BULK );

    tIWG->set_residual_dof_type( tResidualDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_property( tPropViscosity,    "DynamicViscosity" );
    tIWG->set_property( tPropConductivity, "ThermalConductivity" );
    tIWG->set_material_model( tMMFluid, "FluidMM" );
    tIWG->set_constitutive_model( tCMMasterFluid, "FluidCM" );

    // FIXME: generic SP for testing strong form only
    tIWG->set_stabilization_parameter( tSP, "GenericSP" );

    //------------------------------------------------------------------------------
    // set a fem set pointer

    MSI::Equation_Set * tSet = new fem::Set();
    static_cast<fem::Set*>(tSet)->set_set_type( fem::Element_Type::BULK );
    tMMFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCMMasterFluid->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tMMFluid->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tMMFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )   = 0;
    tMMFluid->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 1;

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )     = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) )    = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) )  = 2;

    // set number of spatial dimensions
    uint iSpaceDim = 2;

    // set geometry type
    tGeometryType = mtk::Geometry_Type::QUAD;

    // set velocity dof types
    tVelocityDof = { MSI::Dof_Type::VX, MSI::Dof_Type::VY };

    // set space dimension to CM
    tMMFluid->set_space_dim( iSpaceDim );
    tCMMasterFluid->set_space_dim( iSpaceDim );

    // set interpolation order
    uint iInterpOrder = 2;

    // create an interpolation order
    mtk::Interpolation_Order tGIInterpolationOrder = mtk::Interpolation_Order::LINEAR;

    //------------------------------------------------------------------------------
    // space and time geometry interpolators
    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGIRule( tGeometryType,
            mtk::Interpolation_Type::LAGRANGE,
            tGIInterpolationOrder,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 7.9 }, { 9.3 }};

    Matrix< DDRMat > tXHat = { 
            { -0.4, -3.9 }, 
            {  2.4, -3.9 }, 
            {  2.4, -1.3 }, 
            { -0.4, -1.3 } };

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    //------------------------------------------------------------------------------
    // integration points
    // get an integration order
    mtk::Integration_Order tIntegrationOrder = mtk::Integration_Order::QUAD_2x2;

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

    //------------------------------------------------------------------------------
    // field interpolators
    // create an interpolation order
    mtk::Interpolation_Order tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

    // number of dof for interpolation order
    uint tNumCoeff = 18;

    // get number of dof per type
    int tNumDofRho  = tNumCoeff;
    int tNumDofVel  = tNumCoeff * iSpaceDim;
    int tNumDofTemp = tNumCoeff;
    int tTotalNumDof = tNumDofRho + tNumDofVel + tNumDofTemp;

    //create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule (
            tGeometryType,
            mtk::Interpolation_Type::LAGRANGE,
            tInterpolationOrder,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // fill coefficients for master FI
    Matrix< DDRMat > tMasterDOFHatRho;
    fill_RhoHat( tMasterDOFHatRho, iSpaceDim, iInterpOrder );
    Matrix< DDRMat > tMasterDOFHatVel;
    fill_UHat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
    Matrix< DDRMat > tMasterDOFHatTemp;
    fill_TempHat( tMasterDOFHatTemp, iSpaceDim, iInterpOrder );

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

    // create the field interpolator density
    tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDof );
    tMasterFIs( 0 )->set_coeff( tMasterDOFHatRho );

    // create the field interpolator velocity
    tMasterFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof );
    tMasterFIs( 1 )->set_coeff( tMasterDOFHatVel );

    // create the field interpolator pressure
    tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof );
    tMasterFIs( 2 )->set_coeff( tMasterDOFHatTemp );

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofRho - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofRho, tNumDofRho + tNumDofVel - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };

    // set size and fill the set jacobian assembly map
    Matrix< DDSMat > tJacAssembly = {
            { 0, tNumDofRho - 1 },
            { tNumDofRho, tNumDofRho + tNumDofVel - 1 },
            { tNumDofRho + tNumDofVel, tTotalNumDof - 1 } };
    tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;
    tIWG->mSet->mJacDofAssemblyMap( 2 ) = tJacAssembly;

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size(
            tTotalNumDof,
            1,
            0.0 );
    tIWG->mSet->mJacobian.set_size(
            tTotalNumDof,
            tTotalNumDof,
            0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = tDofTypes;

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
    moris::Cell< moris::Cell< enum mtk::Field_Type > > tDummyField;
    Field_Interpolator_Manager tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tMasterFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager( &tFIManager );

    // set the interpolator manager to the set
    tCMMasterFluid->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tCMMasterFluid->set_field_interpolator_manager( &tFIManager );

    // init the jacobian for IWG
    Matrix<DDRMat> tResidual;
    Matrix<DDRMat> tJacobian;

    // loop over integration points
    // uint tNumGPs = tIntegPoints.n_cols();
    //for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
    for( uint iGP = 0; iGP < 1; iGP ++ )
    {
        // output for debugging
        // std::cout << "-------------------------------------------------------------------\n" << std::flush;
        // std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

        // reset CM evaluation flags
        tCMMasterFluid->reset_eval_flags();

        // reset IWG evaluation flags
        tIWG->reset_eval_flags();

        // create evaluation point xi, tau
        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );       

        // set integration point
        tCMMasterFluid->mSet->mMasterFIManager->set_space_time( tParamPoint );
        tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );

        // check evaluation of the residual for IWG
        //------------------------------------------------------------------------------
        // reset residual & jacobian
        tIWG->mSet->mResidual( 0 ).fill( 0.0 );      
        tIWG->mSet->mJacobian.fill( 0.0 );    

        // compute residual & jacobian
        tIWG->compute_residual( 1.0 );
        tIWG->compute_jacobian( 1.0 );

        // init the jacobian for IWG
        if ( iGP == 0 )
        {
            tResidual = tIWG->mSet->mResidual( 0 );
            tJacobian = tIWG->mSet->mJacobian;
        }
        else
        {
            tResidual += tIWG->mSet->mResidual( 0 );
            tJacobian += tIWG->mSet->mJacobian;
        }

        // debug: print matrices
        //print( trans( tResidual ), "tResidual" );
        //print( tJacobian, "tJacobian" );
    }

    // Analytical residuals from Matlab code
    Matrix< DDRMat > tResidualWeakAnalytical = { {
            -0.00580390,0.001555152,-0.00041670,0.001555152,-0.00849751,0.002276901,0.002276901,-0.00849751,-0.01244122,
            -0.00155515,0.000416701,-0.00011165,0.000416701,-0.00227690,0.000610093,0.000610093,-0.00227690,-0.00333361,
            -1.28897612,0.345438679,-0.09256364,0.345393661,-1.88762917,0.505784425,0.505808550,-1.88729316,-2.76382898,
            -0.34538011,0.092560015,-0.02480235,0.092547952,-0.50578871,0.135524528,0.135530992,-0.50569867,-0.74056574,
            -0.51160558,0.137096883,-0.03671580,0.137012677,-0.74913646,0.200580516,0.200625643,-0.74850793,-1.09602916,
            -0.13708430,0.036734999,-0.00983797,0.036712436,-0.20073051,0.053745387,0.053757479,-0.20056209,-0.29368012,
            0.755137341,-0.17853015,0.068735091,-0.28033113,0.927890368,-0.41737168,-0.36281670,1.687743210,2.210845821,
            0.202338440,-0.04783701,0.018417512,-0.07511450,0.248627474,-0.11183440,-0.09721644,0.452229430,0.592394352 } };

    Matrix< DDRMat > tResidualStrongAnalytical = { {   // = N' * R_strong
            -0.00580390,0.001555152,-0.00041670,0.001555152,-0.00849751,0.002276901,0.002276901,-0.00849751,-0.01244122,
            -0.00155515,0.000416701,-0.00011165,0.000416701,-0.00227690,0.000610093,0.000610093,-0.00227690,-0.00333361,
            -1.28922905,0.345447883,-0.09256248,0.345447883,-1.88756233,0.505770804,0.505770804,-1.88756233,-2.76358306,
            -0.34544788,0.092562481,-0.02480204,0.092562481,-0.50577080,0.135520878,0.135520878,-0.50577080,-0.74049985,
            -0.51139386,0.137027573,-0.03671642,0.137027573,-0.74873258,0.200622291,0.200622291,-0.74873258,-1.09622058,
            -0.13702757,0.036716427,-0.00983813,0.036716427,-0.20062229,0.053756581,0.053756581,-0.20062229,-0.29373142,
            1.066347801,-0.28572703,0.076560327,-0.28572703,1.561241538,-0.41833340,-0.41833340,1.561241538,2.285816258,
            0.285727032,-0.07656032,0.020514277,-0.07656032,0.418333409,-0.11209209,-0.11209209,0.418333409,0.612482620 } };

    Matrix< DDRMat > tResidualAnalytical = trans( tResidualWeakAnalytical ) + tDummyFactor * trans( tResidualStrongAnalytical );

    // check jacobian against analytical solution
    bool tCheckResidual = fem::check( tResidual, tResidualAnalytical, tEpsilon, true );
    REQUIRE( tCheckResidual );

    // clean up
    tMasterFIs.clear();

}/*END_TEST_CASE*/