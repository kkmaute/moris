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
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Compressible_1D.cpp"

#include "fn_trans.hpp"
#include "fn_FEM_Check.hpp"

// debug - output to hdf5
#include "paths.hpp"
#include "HDF5_Tools.hpp"
#include "fn_FEM_Check.hpp"
#include "fn_sqrtmat.hpp"

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Compressible_NS_Bulk_Perfect_Gas_Pressure_Primitive_Analytical",
        "[IWG_Compressible_NS_Bulk_Perfect_Gas_Pressure_Primitive_Analytical]" )
{
    // define an epsilon environment
    real tEpsilon = 1.0E-7;

    // define absolute tolerance accepted as numerical error
    real tAbsTol = 1.0E-12;

    // weight for the GLS term, use 0 for turning GLS off, use 1 for weak form + regular GLS formulation 
    real tGLSWeightFactor = 0.0;

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
    tPropViscosity->set_parameters( { {{ 1.716e-5 }} } );
    tPropViscosity->set_val_function( tConstValFunc );

    // isochoric heat capacity
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 0.718e3 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunc );

    // specific gas constant
    std::shared_ptr< fem::Property > tPropGasConstant = std::make_shared< fem::Property >();
    tPropGasConstant->set_parameters( { {{ 287.058 }} } );
    tPropGasConstant->set_val_function( tConstValFunc );

    // thermal conductivity
    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 24.35e-3 }} } );
    tPropConductivity->set_val_function( tConstValFunc );

    // Heat Loading
    std::shared_ptr< fem::Property > tPropHeatLoad = std::make_shared< fem::Property >();
    tPropHeatLoad->set_parameters( { {{ 1.0e5 }} } ); // maximum volumetric heat load in center of channel in W/m^3
    tPropHeatLoad->set_val_function( Func_Heat_Load_Distribution );

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
    tSP->set_parameters( { {{ tGLSWeightFactor }} });
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
    //uint iInterpOrder = 2;

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

    // fill in data
    fill_data( tXHat, tTHat, tMasterDOFHatP, tMasterDOFHatVel, tMasterDOFHatTemp );

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

    // init the jacobian for IWG
    Matrix<DDRMat> tResidual;
    Matrix<DDRMat> tJacobian;

    // get number of integration points
    uint tNumGPs = tIntegPoints.n_cols();
    std::cout << "Number of Gauss points: " << tNumGPs << " \n" << std::flush;

    // initialize matrix storing individual Gausspoint residuals
    Matrix< DDRMat > tGPResiduals( 18, tNumGPs, 0.0 );

    // initialize cell storing individual Gausspoint residuals
    moris::Cell< Matrix< DDRMat > > tGPJacobians( tNumGPs );

    // loop over integration points
    //for( uint iGP = 0; iGP < tNumGPs - tNumGPs + 1; iGP ++ ) // uncomment for debug - only use first GP
    for( uint iGP = 0; iGP < tNumGPs; iGP ++ )
    {
        // output for debugging
        std::cout << "-------------------------------------------------------------------\n" << std::flush;
        std::cout << "Looping over Gauss points. Current GP-#: " << iGP << "\n\n" << std::flush;

        // reset CM evaluation flags
        tCMMasterFluid->reset_eval_flags();

        // reset IWG evaluation flags
        tIWG->reset_eval_flags();

        // create evaluation point xi, tau
        Matrix< DDRMat > tParamPoint = tIntegPoints.get_column( iGP );
        // Matrix< DDRMat > tParamPoint = {
        //         {-0.67},
        //         {+0.22},
        //         {+0.87}};      

        // set integration point
        tCMMasterFluid->mSet->mMasterFIManager->set_space_time( tParamPoint );
        tIWG->mSet->mMasterFIManager->set_space_time( tParamPoint );

        // for debug
        // print( tIWG->mSet->mMasterFIManager->get_IP_geometry_interpolator()->valx(), "x-pos" );
        // print( tIWG->mSet->mMasterFIManager->get_IP_geometry_interpolator()->valt(), "t-pos" ); 

        // check evaluation of the residual for IWG
        //------------------------------------------------------------------------------
        // reset residual & jacobian
        tIWG->mSet->mResidual( 0 ).fill( 0.0 );      
        tIWG->mSet->mJacobian.fill( 0.0 );    

        // compute detJ of integration domain
        real tDetJ = tIWG->mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->det_J();

        // compute residual & jacobian
        tIWG->compute_residual( tIntegWeights( iGP ) * tDetJ );
        tIWG->compute_jacobian( tIntegWeights( iGP ) * tDetJ );

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

        // temporarily store Gauss point residual and jacobian
        Matrix< DDRMat > tResidualGP = tIWG->mSet->mResidual( 0 );
        Matrix< DDRMat > tJacobianGP = tIWG->mSet->mJacobian;

        // convert the 2D Gauss point residual to 1D
        Matrix< DDRMat > tResidualGP1D = fem::convert_comp_flow_residual_2D_to_1D_quadratic( tResidualGP );
        Matrix< DDRMat > tJacobianGP1D = fem::convert_comp_flow_jacobian_2D_to_1D_quadratic( tJacobianGP );

        // store GP residual and jacobian to list of all residuals and jacobians
        tGPResiduals( { 0, 17 }, { iGP } ) = tResidualGP1D.matrix_data();
        tGPJacobians( iGP ) = tJacobianGP1D;

    } // end Gauss loop

    //------------------------------------------------------------------------------

    // convert Gauss sum of the the 2D residual and jabobian to 1D
    Matrix< DDRMat > tResidual1D = fem::convert_comp_flow_residual_2D_to_1D_quadratic( tResidual );
    Matrix< DDRMat > tJacobian1D = fem::convert_comp_flow_jacobian_2D_to_1D_quadratic( tJacobian );

    // Analytical residual for the weak form from Matlab code
    Matrix< DDRMat > tResWeakFormAna1D = { 
                { -2.719421205234789e-02 },
                {  3.254348159384793e-02 },
                {  7.661289950608770e-02 },
                { -2.228817689751322e-02 },
                {  3.887942571473649e-02 },
                {  1.091741631355034e-01 },
                {  1.719141550024611e+02 },
                { -1.399967373612037e+02 },
                {  6.381144789420286e+01 },
                {  2.306186859049692e+02 },
                { -1.819219563676959e+02 },
                {  9.736792124830566e+01 },
                { -7.569324106788088e+02 },
                {  1.314743811630449e+03 },
                {  6.418214259941517e+03 },
                { -5.270194544693870e+02 },
                {  1.486446479850517e+03 },
                {  7.643294989136290e+03 } };

    // Analytical Jacobian for the weak form from Matlab code
    Matrix< DDRMat > tJacWeakFormAna1D = {
            { -5.970568286e-07,  1.397150775e-07, -3.039911886e-07,  4.147316760e-07, -1.439087247e-07,  1.817415458e-07, -3.008716944e-03, -1.265553707e-03,  3.895610866e-03, -1.028950874e-03, -6.305310878e-04,  2.022074088e-03,  1.217241697e-04, -2.667272093e-05,  6.299189495e-05, -5.963350162e-05,  3.251173121e-05, -1.858890040e-05 },
            {  1.428850211e-07, -5.188488160e-07, -2.708721752e-07, -1.405014718e-07,  7.381560071e-07,  3.769142495e-07,  1.062960666e-03,  5.343974025e-03, -4.507014728e-03,  5.481278756e-04,  2.559265445e-03, -2.657735479e-03, -2.747991670e-05,  9.141400835e-05,  4.762829528e-05,  3.157538049e-05, -2.006351186e-04, -1.069718740e-04 },
            { -3.055132566e-07, -2.500237740e-07, -1.563761427e-06,  1.798582628e-07,  3.991267281e-07,  2.167282894e-06, -4.061575968e-03,  5.483351995e-03,  7.066116427e-04, -1.888946975e-03,  2.854221214e-03,  6.613003385e-04,  6.335502946e-05,  4.263533085e-05,  2.107882006e-04, -1.807259106e-05, -1.127150812e-04, -4.836794749e-04 },
            { -4.637417653e-07,  1.374933429e-07, -2.163060365e-07,  3.299827088e-07, -1.571308775e-07,  1.196817608e-07, -1.028950874e-03, -6.305310878e-04,  2.022074088e-03, -1.329462828e-03, -1.227101362e-03,  4.159950698e-03,  7.580719134e-05, -2.851324537e-05,  3.102157719e-05, -3.558811220e-05,  4.289527560e-05,  2.933109058e-06 },
            {  1.409005958e-07, -6.345870933e-07, -3.506464389e-07, -1.461406514e-07,  9.313837277e-07,  4.546744167e-07,  5.481278756e-04,  2.559265445e-03, -2.657735479e-03,  1.189443996e-03,  4.766328852e-03, -6.275460337e-03, -2.944959608e-05,  1.521240011e-04,  8.751755466e-05,  3.964617121e-05, -3.079156972e-04, -1.591239879e-04 },
            { -2.181893195e-07, -3.284339603e-07, -1.924034849e-06,  1.128946959e-07,  5.258352424e-07,  2.757556933e-06, -1.888946975e-03,  2.854221214e-03,  6.613003385e-04, -3.666849356e-03,  5.932988461e-03,  2.186613754e-03,  3.153788653e-05,  8.177434746e-05,  3.679359497e-04,  4.991724062e-06, -1.789666324e-04, -7.999437580e-04 },
            { -5.032111561e-03, -1.677389357e-03,  6.709444660e-03, -2.516070000e-03, -8.387324088e-04,  3.354688708e-03, -4.843304494e-02,  1.116808808e-02, -2.467160121e-02,  3.316996367e-02, -1.227850951e-02,  1.415535237e-02,  1.418213679e-06,  6.027158245e-06,  4.729673595e-06,  1.766915552e-06,  1.425514700e-05,  1.029002922e-05 },
            {  1.677345123e-03,  5.032292989e-03, -6.709359186e-03,  8.386361078e-04,  2.516462316e-03, -3.354500846e-03,  1.171268006e-02, -4.062446927e-02, -2.273343600e-02, -1.168553478e-02,  6.530958953e-02,  3.202392123e-02,  5.727188423e-06, -4.958458793e-05, -2.747237521e-05,  1.374929922e-05, -1.157269206e-04, -6.500410823e-05 },
            { -6.709489076e-03,  6.709580539e-03,  2.484191177e-07, -3.354778749e-03,  3.354976092e-03,  5.484468764e-07, -2.498776375e-02, -1.913211379e-02, -1.290493193e-01,  1.376086580e-02,  3.591511671e-02,  2.035229627e-01,  4.841360189e-06, -2.911604524e-05, -5.409821442e-05,  1.048185479e-05, -6.778556072e-05, -1.344781449e-04 },
            { -2.516038966e-03, -8.386629392e-04,  3.354755710e-03, -5.032132337e-03, -1.677481341e-03,  6.709374752e-03, -3.703593426e-02,  1.124506160e-02, -1.696619883e-02,  2.615800238e-02, -1.407598905e-02,  8.906484493e-03, -1.328992109e-06, -5.382120559e-06, -4.865420929e-06,  7.658343402e-07,  3.756004333e-05,  2.381393646e-05 },
            {  8.387055774e-04,  2.515875524e-03, -3.354845728e-03,  1.677259817e-03,  5.033048422e-03, -6.708951698e-03,  1.183803633e-02, -5.267308148e-02, -3.094881202e-02, -1.213992346e-02,  8.710841574e-02,  3.792369985e-02, -5.887968343e-06,  4.446712752e-05,  2.870041219e-05,  3.530134667e-05, -3.008999047e-04, -1.637443762e-04 },
            { -3.354711747e-03,  3.354631211e-03, -2.436585004e-07, -6.709562079e-03,  6.710025120e-03,  1.285739421e-06, -1.736068540e-02, -2.705761655e-02, -1.730154174e-01,  7.472468683e-03,  5.055915489e-02,  2.758975531e-01, -4.673595360e-06,  2.591895970e-05,  5.293837444e-05,  2.468554779e-05, -1.762006680e-04, -3.699762347e-04 },
            { -1.651307936e-01,  4.098202879e-02, -8.273629481e-02,  1.683673934e-01, -4.239252616e-02,  8.401275017e-02, -1.189975074e+03, -4.872906410e+02,  2.044311715e+03, -5.577611041e+02, -2.454286564e+02,  1.038091550e+03,  5.715947041e-04,  8.308114226e-05, -6.527669600e-04,  2.851985340e-04,  4.555908322e-05, -3.241876127e-04 },
            {  4.204211114e-02, -1.617096769e-01, -8.681604444e-02, -4.133246551e-02,  1.717886185e-01,  7.993320298e-02,  4.808291359e+02,  1.312686692e+03, -2.095965769e+03,  2.407243779e+02,  6.470163033e+02, -1.075641717e+03,  8.301771790e-05,  5.608945315e-04, -6.595198623e-04,  4.541531965e-05,  2.489118287e-04, -3.475161681e-04 },
            { -8.339228781e-02, -7.915193267e-02, -6.641931716e-01,  8.335684440e-02,  8.759716914e-02,  6.697997701e-01, -1.897495136e+03,  1.974995211e+03,  5.168429712e+01, -9.441345815e+02,  1.000501861e+03,  3.761837742e+01, -6.527548661e-04, -6.598598916e-04,  1.299860780e-03, -3.241592882e-04, -3.482876121e-04,  6.282226402e-04 },
            { -1.651308350e-01,  4.098205112e-02, -8.273635063e-02,  1.716039167e-01, -4.380301511e-02,  8.528907453e-02, -5.577587135e+02, -2.454225844e+02,  1.038097192e+03, -1.041066503e+03, -4.944190071e+02,  2.108059863e+03,  2.861251020e-04,  3.964643893e-05, -3.275116635e-04,  5.694237006e-04,  9.880044737e-05, -6.442734713e-04 },
            {  4.204211177e-02, -1.617097712e-01, -8.681600536e-02, -4.062285399e-02,  1.818676550e-01,  7.305060005e-02,  2.407304499e+02,  6.469653364e+02, -1.075671647e+03,  4.820735666e+02,  1.275335828e+03, -2.206627357e+03,  3.950267536e-05,  2.952864151e-04, -3.205415401e-04,  9.801811936e-05,  4.374002590e-04, -7.275951310e-04 },
            { -8.339225640e-02, -7.915203920e-02, -6.641932555e-01,  8.332144424e-02,  9.604235759e-02,  6.754064211e-01, -9.441289397e+02,  1.000471931e+03,  3.754419860e+01, -1.879037923e+03,  2.026987178e+03,  9.873055590e+01, -3.274833391e-04, -3.213129840e-04,  6.604497260e-04, -6.441162938e-04, -7.317955102e-04,  1.212959588e-03 }
            };

    // Analytical residual for the GLS term from Matlab code
    Matrix< DDRMat > tResStabAna1D = 0.0 * tResWeakFormAna1D;

    // Analytical Jacobian for the GLS term from Matlab code
    Matrix< DDRMat > tJacStabAna1D = 0.0 * tJacWeakFormAna1D;

    // check residual against 1D analytical solution
    Matrix< DDRMat > tResAna1D = tResWeakFormAna1D + tGLSWeightFactor * tResStabAna1D;
    bool tCheckResidual = fem::check( tResidual1D, tResAna1D, tEpsilon, true, true, tAbsTol );
    std::cout << "Residual check passed? : " << tCheckResidual << " \n" << std::flush;
    REQUIRE( tCheckResidual );

    // check Jacobian against 1D analytical solution
    Matrix< DDRMat > tJacAna1D = tJacWeakFormAna1D + tGLSWeightFactor * tJacStabAna1D;
    bool tCheckJacobian = fem::check( tJacobian1D, tJacAna1D, tEpsilon, true, true, tAbsTol );
    std::cout << "Jacobian check passed? : " << tCheckJacobian << " \n" << std::flush;
    REQUIRE( tCheckJacobian );

    //------------------------------------------------------------------------------

    // clean up
    tMasterFIs.clear();

}/*END_TEST_CASE*/