#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_IWG_Compressible_NS_Dirichlet_Nitsche.hpp"
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

TEST_CASE( "IWG_Compressible_NS_Dirichlet_Analytical",
        "[IWG_Compressible_NS_Dirichlet_Analytical]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-10;

    // define absolute tolerance accepted as numerical error
    real tAbsTol = 1.0E-12;

    // weight for the GLS term, use 0 for turning GLS off, use 1 for weak form + regular GLS formulation 
    real tUpwindFactor = 0.0;

    // penalty factors
    real tPpenatly = 1.0;
    real tUXpenatly = 1.0;
    real tUYpenatly = 1.0;
    real tTpenatly = 1.0;

    // prescribed values
    real tPpresc = 151.2;
    real tUXpresc =  0.0;
    real tUYpresc =  0.0;
    real tTpresc = 121.7;

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

    //------------------------------------------------------------------------------
    // prescribed values

    // prescribe pressure
    std::shared_ptr< fem::Property > tPropPrescPres = std::make_shared< fem::Property >();
    tPropPrescPres->set_parameters( { {{ tPpresc }} } );
    tPropPrescPres->set_val_function( tConstValFunc );

    // prescribed velocity
    std::shared_ptr< fem::Property > tPropPrescVel = std::make_shared< fem::Property >();
    tPropPrescVel->set_parameters( { {{ tUXpresc },{ tUYpresc }} } );
    tPropPrescVel->set_val_function( tConstValFunc );

    // prescribed temperature
    std::shared_ptr< fem::Property > tPropPrescTemp = std::make_shared< fem::Property >();
    tPropPrescTemp->set_parameters( { {{ tTpresc }} } );
    tPropPrescTemp->set_val_function( tConstValFunc );

    // upwinding weight
    std::shared_ptr< fem::Property > tPropUpwind = std::make_shared< fem::Property >();
    tPropUpwind->set_parameters( { {{ tUpwindFactor }} } );
    tPropUpwind->set_val_function( tConstValFunc );


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

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche =
            tSPFactory.create_SP( fem::Stabilization_Type::COMPRESSIBLE_DIRICHLET_NITSCHE );
    tSPNitsche->set_property( tPropViscosity, "DynamicViscosity", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropConductivity, "ThermalConductivity", mtk::Master_Slave::MASTER );


    // create a dummy fem cluster and set it to SP
    fem::Cluster * tCluster = new fem::Cluster();
    tSPNitsche->set_cluster( tCluster );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE );

    tIWG->set_residual_dof_type( tResidualDofTypes );
    tIWG->set_dof_type_list( tDofTypes, mtk::Master_Slave::MASTER );
    tIWG->set_property( tPropViscosity,    "DynamicViscosity" );
    tIWG->set_property( tPropConductivity, "ThermalConductivity" );
    tIWG->set_material_model( tMMFluid, "FluidMM" );
    tIWG->set_constitutive_model( tCMMasterFluid, "FluidCM" );
    tIWG->set_stabilization_parameter( tSPNitsche, "NitschePenaltyParameter" );

    tIWG->set_property( tPropPrescPres, "PrescribedDof1" );
    tIWG->set_property( tPropPrescVel,  "PrescribedVelocity" );
    tIWG->set_property( tPropPrescTemp, "PrescribedDof3" );
    tIWG->set_property( tPropUpwind,    "PressureUpwind" );

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
    tSPNitsche->set_space_dim( iSpaceDim );

    // create normal for Nitsche
    Matrix< DDRMat > tNormal = {{3.0/5.0},{4.0/5.0}};

    // assign normal to IWG
    tIWG->set_normal( tNormal );

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
    fem::IWG_Compressible_NS_Dirichlet_Nitsche * tChildIWG = dynamic_cast< fem::IWG_Compressible_NS_Bulk * > ( tIWG.get() );

    // init the jacobian for IWG
    Matrix< DDRMat > tResidual = tIWG->mSet->get_residual()( 0 );
    Matrix< DDRMat > tJacobian = tIWG->mSet->get_jacobian();

    //------------------------------------------------------------------------------
    // get reference solution
    Matrix< DDRMat > tRefResidual = get_reference_ResVelNitsche() + get_reference_ResTmpNitsche() + tUpwindFactor * get_reference_ResUpwind();
    Matrix< DDRMat > tRefJacobian = get_reference_JacVelNitsche() + get_reference_JacTmpNitsche() + tUpwindFactor * get_reference_JacUpwind();

    //------------------------------------------------------------------------------
    // check residual against 2D analytical solution
    bool tCheckResidual = fem::check( tResidual, tRefResidual, tEpsilon, true, true, tAbsTol );
    std::cout << "Residual check passed? : " << tCheckResidual << " \n" << std::flush;
    REQUIRE( tCheckResidual );

    // check Jacobian against 2D analytical solution
    bool tCheckJacobian = fem::check( tJacobian, tRefJacobian, tEpsilon, true, true, tAbsTol );
    std::cout << "Jacobian check passed? : " << tCheckJacobian << " \n" << std::flush;
    REQUIRE( tCheckJacobian );

    //------------------------------------------------------------------------------

    // clean up
    tMasterFIs.clear();

}/*END_TEST_CASE*/