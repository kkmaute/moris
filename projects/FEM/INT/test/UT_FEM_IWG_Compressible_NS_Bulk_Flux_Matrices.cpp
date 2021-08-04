#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IWG_Compressible_NS_Bulk.hpp"
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
//#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_1D.cpp"
#include "FEM_Test_Proxy/cl_FEM_Fields_for_NS_Compressible_UT.cpp"
#include "FEM_Test_Proxy/cl_FEM_Flux_Matrix_Reference_Values.cpp"

#include "fn_trans.hpp"
#include "fn_FEM_Check.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"
#include "fn_FEM_Check.hpp"
#include "fn_sqrtmat.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Compressible_NS_Bulk_Flux_Matrices",
        "[IWG_Compressible_NS_Bulk_Flux_Matrices]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-12;

    // define absolute tolerance accepted as numerical error
    real tAbsTol = 1.0E-14;

    // init geometry inputs
    //------------------------------------------------------------------------------
    // create geometry type
    mtk::Geometry_Type tGeometryType = mtk::Geometry_Type::UNDEFINED;

    // dof type list
    moris::Cell< MSI::Dof_Type > tPressureDof = { MSI::Dof_Type::P };
    moris::Cell< MSI::Dof_Type > tVelocityDof = { MSI::Dof_Type::VX };
    moris::Cell< MSI::Dof_Type > tTempDof     = { MSI::Dof_Type::TEMP };

    moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes         = { tPressureDof, tVelocityDof, tTempDof };
    moris::Cell< moris::Cell< MSI::Dof_Type > > tResidualDofTypes = tDofTypes;

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties

    // dynamic viscosity
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.2 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 5.7 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 2.4 }} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 0.8 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    // Body heat load
    std::shared_ptr< fem::Property > tPropHeatLoad = std::make_shared< fem::Property >();
    tPropHeatLoad->set_parameters( { {{ 1.3 }} } );
    tPropHeatLoad->set_val_function( tConstValFunc );

    // dummy for SP
    std::shared_ptr< fem::Property > tPropDummyForSP = std::make_shared< fem::Property >();
    tPropDummyForSP->set_parameters( { {{ 1.0 }} } );
    tPropDummyForSP->set_val_function( tConstValFunc );

    // define material model and assign properties
    fem::MM_Factory tMMFactory;

    std::shared_ptr< fem::Material_Model > tMMFluid =
            tMMFactory.create_MM( fem::Material_Type::PERFECT_GAS );
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
    tSP->set_parameters( { {{ 1.8 }} });
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
    tIWG->set_property( tPropHeatLoad, "BodyHeatLoad" ); 
    tIWG->set_material_model( tMMFluid, "FluidMM" );
    tIWG->set_constitutive_model( tCMMasterFluid, "FluidCM" );

    // FIXME: generic SP for testing strong form only
    tIWG->set_stabilization_parameter( tSP, "GLS" );

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
    mtk::Interpolation_Order tGIInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

    //------------------------------------------------------------------------------
    // prepare element information

    // initialize matrices to be filled
    Matrix< DDRMat >tXHat;
    Matrix< DDRMat >tTHat;
    Matrix< DDRMat >tMasterDOFHatP;
    Matrix< DDRMat >tMasterDOFHatVel;
    Matrix< DDRMat >tMasterDOFHatTemp;

    // fill in data for element size
    fill_data_rectangle_element( tXHat, tTHat );

    // fill DoF values
    fill_smooth_PHat( tMasterDOFHatP, iSpaceDim, iInterpOrder );
    tMasterDOFHatP = trans( tMasterDOFHatP );
    fill_smooth_UHat( tMasterDOFHatVel, iSpaceDim, iInterpOrder );
    fill_smooth_TempHat( tMasterDOFHatTemp, iSpaceDim, iInterpOrder );
    tMasterDOFHatTemp = trans( tMasterDOFHatTemp ); 

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

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    //------------------------------------------------------------------------------
    // integration points
    // get an integration order
    mtk::Integration_Order tIntegrationOrder = mtk::Integration_Order::QUAD_3x3;

    // create an integration rule
    mtk::Integration_Rule tIntegrationRule(
            tGeometryType,
            mtk::Integration_Type::GAUSS,
            tIntegrationOrder,
            mtk::Geometry_Type::LINE,
            mtk::Integration_Type::GAUSS,
            mtk::Integration_Order::BAR_3 );

    // create an integrator
    mtk::Integrator tIntegrator( tIntegrationRule );

    // get integration points
    Matrix< DDRMat > tIntegPoints;
    tIntegrator.get_points( tIntegPoints );

    // get integration weights
    Matrix< DDRMat > tIntegWeights;
    tIntegrator.get_weights( tIntegWeights );


    //------------------------------------------------------------------------------
    // field interpolators
    // create an interpolation order
    mtk::Interpolation_Order tInterpolationOrder = mtk::Interpolation_Order::QUADRATIC;

    // number of dof for interpolation order
    uint tNumCoeff = 18;

    // get number of dof per type
    int tNumDofP  = tNumCoeff;
    int tNumDofVel  = tNumCoeff * iSpaceDim;
    int tNumDofTemp = tNumCoeff;
    int tTotalNumDof = tNumDofP + tNumDofVel + tNumDofTemp;

    //create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule (
            tGeometryType,
            mtk::Interpolation_Type::LAGRANGE,
            tInterpolationOrder,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tMasterFIs( tDofTypes.size() );

    // create the field interpolator density
    tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, tPressureDof );
    tMasterFIs( 0 )->set_coeff( tMasterDOFHatP );

    // create the field interpolator velocity
    tMasterFIs( 1 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tVelocityDof );
    tMasterFIs( 1 )->set_coeff( tMasterDOFHatVel );

    // create the field interpolator pressure
    tMasterFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI, tTempDof );
    tMasterFIs( 2 )->set_coeff( tMasterDOFHatTemp );

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofP - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofP, tNumDofP + tNumDofVel - 1 } };
    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { tNumDofP + tNumDofVel, tTotalNumDof - 1 } };

    // set size and fill the set jacobian assembly map
    Matrix< DDSMat > tJacAssembly = {
            { 0, tNumDofP - 1 },
            { tNumDofP, tNumDofP + tNumDofVel - 1 },
            { tNumDofP + tNumDofVel, tTotalNumDof - 1 } };
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

    // reset IWG evaluation flags
    tIWG->reset_eval_flags();

    // create evaluation point xi, tau
    Matrix< DDRMat > tParamPoint = {
            {-0.67},
            {+0.22},
            {+0.87}};      

    // set integration point
    tCMMasterFluid->mSet->mMasterFIManager->set_space_time( tParamPoint );
    tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );

    // for debug
    print( tIWG->mSet->mMasterFIManager->get_IP_geometry_interpolator()->valx(), "x-pos" );
    print( tIWG->mSet->mMasterFIManager->get_IP_geometry_interpolator()->valt(), "t-pos" ); 

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // reset residual & jacobian
    tIWG->mSet->mResidual( 0 ).fill( 0.0 );      
    tIWG->mSet->mJacobian.fill( 0.0 );    

    // compute residual & jacobian
    tIWG->compute_residual( 1.0 );
    tIWG->compute_jacobian( 1.0 );

    // get the child
    fem::IWG_Compressible_NS_Bulk * tChildIWG = dynamic_cast< fem::IWG_Compressible_NS_Bulk * > ( tIWG.get() );

    //------------------------------------------------------------------------------
    // check flux matrices

    // state vars and derivs
    std::cout << "Checking State Variable Vector and its Derivatives ... \n" << std::flush;
    Matrix< DDRMat > tY     = tChildIWG->Y();
    Matrix< DDRMat > tY_ref = get_reference_Y();
    REQUIRE( fem::check( tY, tY_ref, tEpsilon, true, true, tAbsTol ) );

    Matrix< DDRMat > tdYdt     = tChildIWG->dYdt();
    Matrix< DDRMat > tdYdt_ref = get_reference_dYdt();
    REQUIRE( fem::check( tdYdt, tdYdt_ref, tEpsilon, true, true, tAbsTol ) );

    Matrix< DDRMat > tdYdx     = tChildIWG->dYdx( 0 );
    Matrix< DDRMat > tdYdx_ref = get_reference_dYdx()( 0 );
    REQUIRE( fem::check( tdYdx, tdYdx_ref, tEpsilon, true, true, tAbsTol ) );

    Matrix< DDRMat > tdYdy     = tChildIWG->dYdx( 1 );
    Matrix< DDRMat > tdYdy_ref = get_reference_dYdx()( 1 );
    REQUIRE( fem::check( tdYdy, tdYdy_ref, tEpsilon, true, true, tAbsTol ) );
    
    // Force term
    std::cout << "Checking C-Operator ... \n" << std::flush;
    Matrix< DDRMat > tC     = tChildIWG->C();
    Matrix< DDRMat > tC_ref = get_reference_C();
    REQUIRE( fem::check( tC, tC_ref, tEpsilon, true, true, tAbsTol ) );

    // A - Flux Matrices
    std::cout << "Checking A-Matrices ... \n" << std::flush;
    Cell< Matrix< DDRMat > > tA     = { tChildIWG->A( 0 ), tChildIWG->A( 1 ), tChildIWG->A( 2 ) };
    Cell< Matrix< DDRMat > > tA_ref = get_reference_A();
    REQUIRE( fem::check( tA( 0 ), tA_ref( 0 ), tEpsilon, true, true, tAbsTol ) );
    REQUIRE( fem::check( tA( 1 ), tA_ref( 1 ), tEpsilon, true, true, tAbsTol ) );
    REQUIRE( fem::check( tA( 2 ), tA_ref( 2 ), tEpsilon, true, true, tAbsTol ) );

    // K - Flux Matrices
    std::cout << "Checking K-Matrices ... \n" << std::flush;
    Cell< Matrix< DDRMat > > tK     = { tChildIWG->K( 0, 0 ), tChildIWG->K( 0, 1 ), tChildIWG->K( 1, 0 ), tChildIWG->K( 1, 1 ) };
    Cell< Matrix< DDRMat > > tK_ref = get_reference_K();
    REQUIRE( fem::check( tK( 0 ), tK_ref( 0 ), tEpsilon, true, true, tAbsTol ) );
    REQUIRE( fem::check( tK( 1 ), tK_ref( 1 ), tEpsilon, true, true, tAbsTol ) );
    REQUIRE( fem::check( tK( 2 ), tK_ref( 2 ), tEpsilon, true, true, tAbsTol ) );
    REQUIRE( fem::check( tK( 3 ), tK_ref( 3 ), tEpsilon, true, true, tAbsTol ) );

    // For Stabilization
    std::cout << "Checking LY (strong form of the residual) ... \n" << std::flush;
    Matrix< DDRMat > tLY     = tChildIWG->LY();
    Matrix< DDRMat > tLY_ref = get_reference_LY();
    REQUIRE( fem::check( tLY, tLY_ref, tEpsilon, true, true, tAbsTol ) );

    std::cout << "Checking inv(A0) ... \n" << std::flush;
    Matrix< DDRMat > tA0inv     = tChildIWG->A0inv();
    Matrix< DDRMat > tA0inv_ref = get_reference_A0inv();
    REQUIRE( fem::check( tA0inv, tA0inv_ref, tEpsilon, true, true, tAbsTol ) );

    std::cout << "Checking G ... \n" << std::flush;
    Matrix< DDRMat > tG     = tChildIWG->G();
    Matrix< DDRMat > tG_ref = get_reference_G();
    REQUIRE( fem::check( tG, tG_ref, tEpsilon, true, true, tAbsTol ) );

    std::cout << "Checking M ... \n" << std::flush;
    Matrix< DDRMat > tM     = tChildIWG->M();
    Matrix< DDRMat > tM_ref = get_reference_M();
    REQUIRE( fem::check( tM, tM_ref, tEpsilon, true, true, tAbsTol ) );

    // // initialize storage matrices
    // Matrix<DDRMat> tW = tChildIWG->W();
    // Cell< Matrix< DDRMat > > tdWdx = { tChildIWG->dWdt(), tChildIWG->dWdx( 0 ), tChildIWG->dWdx( 1 ) };
    // Matrix< DDRMat > tC = tChildIWG->C();

    // Matrix< DDRMat > tA0inv = tChildIWG->A0inv();
    // Matrix< DDRMat > tM = tChildIWG->M();
    // Matrix< DDRMat > tTau = tChildIWG->Tau();
    // Matrix< DDRMat > tLW = tChildIWG->LW();
    // Matrix< DDRMat > tGLSTestFunc = tChildIWG->GLSTestFunc();
    // Cell< Matrix< DDRMat > > tdMdY = { tChildIWG->dMdY( 0 ), tChildIWG->dMdY( 1 ), tChildIWG->dMdY( 2 ), tChildIWG->dMdY( 3 ) };
    // //Cell< Matrix< DDRMat > > tdTaudY = { tChildIWG->dTaudY( 0 ), tChildIWG->dTaudY( 1 ), tChildIWG->dTaudY( 2 ), tChildIWG->dTaudY( 3 ) };

    //------------------------------------------------------------------------------

    // clean up
    tMasterFIs.clear();

}/*END_TEST_CASE*/