#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                //FEM//INT/src

#include "op_equal_equal.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                   //FEM//INT//src

void tConstValFunction_NSPBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tTEMPFIValFunction_NSPBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tTEMPFIDerFunction_NSPBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

void tPFIValFunction_NSPBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::P )->val();
}

void tPFIDerFunction_NSPBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::P )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Pressure_Bulk_2D", "[IWG_Incompressible_NS_Pressure_Bulk_2D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropViscosity->set_val_function( tPFIValFunction_NSPBULK );
    tPropViscosity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 2.0 }} } );
    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropDensity->set_val_function( tPFIValFunction_NSPBULK );
    tPropDensity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
    tPropGravity->set_parameters( { {{ 10.0 }, { 10.0 }} } );
    tPropGravity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropGravity->set_val_function( tPFIValFunction_NSPBULK );
    tPropGravity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropThermalExp = std::make_shared< fem::Property >();
    tPropThermalExp->set_parameters( { {{ 23.0 }} } );
    tPropThermalExp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropThermalExp->set_val_function( tTEMPFIValFunction_NSPBULK );
    tPropThermalExp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropRefTemp = std::make_shared< fem::Property >();
    tPropRefTemp->set_parameters( { {{ 15.0 }} } );
    tPropRefTemp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropRefTemp->set_val_function( tTEMPFIValFunction_NSPBULK );
    tPropRefTemp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSPBULK } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMMasterIncFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
    tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
    tCMMasterIncFluid->set_property( tPropDensity, "Density" );
    tCMMasterIncFluid->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
    tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_parameters( { {{ 36.0 }} } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::P } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
    tIWG->set_property( tPropDensity, "Density" );
    tIWG->set_property( tPropGravity, "Gravity" );
    tIWG->set_property( tPropThermalExp, "ThermalExpansion" );
    tIWG->set_property( tPropRefTemp, "ReferenceTemp" );
    tIWG->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.1 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::QUAD,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                              { 1.0, 0.0 },
                              { 1.0, 1.0 },
                              { 0.0, 1.0}};

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }, {1.0 }};

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI.set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tFIRule ( mtk::Geometry_Type::QUAD,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 2 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    arma::Mat< double > tMatrixP;
    tMatrixP.randu( 8, 1 );
    Matrix< DDRMat > tDOFHatP;
    tDOFHatP.matrix_data() = 10.0 * tMatrixP;

    arma::Mat< double > tMatrixT;
    tMatrixT.randu( 8, 1 );
    Matrix< DDRMat > tDOFHatT;
    tDOFHatT.matrix_data() = 10.0 * tMatrixT;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 3 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI,{ { MSI::Dof_Type::VX, MSI::Dof_Type::VY } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ { MSI::Dof_Type::P } } );
    tFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ { MSI::Dof_Type::TEMP } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatP );
    tFIs( 2 )->set_coeff( tDOFHatT );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );
    tFIs( 2 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 2;

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( 3 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 15 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 16, 23 } };
    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { 24, 31 } };

    // set size and fill the set jacobian assembly map
    tIWG->mSet->mJacDofAssemblyMap.resize( 3 );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 15 }, { 16, 23 }, { 24, 31 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 15 }, { 16, 23 }, { 24, 31 } };
    tIWG->mSet->mJacDofAssemblyMap( 2 ) = { { 0, 15 }, { 16, 23 }, { 24, 31 } };

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 32, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 32, 32, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P }, { MSI::Dof_Type::TEMP }};

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum PDV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager(&tFIManager);

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    tIWG->compute_residual( 1.0 );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tJacobians;
    Matrix< DDRMat > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                1.0,
                                                tJacobians,
                                                tJacobiansFD );

//    print(tJacobians({ 0, 7 }, { 0, 15 }),   "tJacobiansPV" );
//    print(tJacobiansFD({ 0, 7 }, { 0, 15 }), "tJacobiansFDPV" );
//    print(tJacobians({ 0, 7 }, { 16, 23 }),   "tJacobiansPP" );
//    print(tJacobiansFD({ 0, 7 }, { 16, 23 }), "tJacobiansFDPP" );
//    print(tJacobians({ 0, 7 }, { 24, 31 }),   "tJacobiansPT" );
//    print(tJacobiansFD({ 0, 7 }, { 24, 31 }), "tJacobiansFDPT" );

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Pressure_Bulk_3D", "[IWG_Incompressible_NS_Pressure_Bulk_3D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropViscosity->set_val_function( tPFIValFunction_NSPBULK );
    tPropViscosity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 2.0 }} } );
    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropDensity->set_val_function( tPFIValFunction_NSPBULK );
    tPropDensity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
    tPropGravity->set_parameters( { {{ 10.0 }, { 10.0 }, { 10.0 }} } );
    tPropGravity->set_dof_type_list( {{ MSI::Dof_Type::P }} );
    tPropGravity->set_val_function( tPFIValFunction_NSPBULK );
    tPropGravity->set_dof_derivative_functions( { tPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropThermalExp = std::make_shared< fem::Property >();
    tPropThermalExp->set_parameters( { {{ 23.0 }} } );
    tPropThermalExp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropThermalExp->set_val_function( tTEMPFIValFunction_NSPBULK );
    tPropThermalExp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSPBULK } );

    std::shared_ptr< fem::Property > tPropRefTemp = std::make_shared< fem::Property >();
    tPropRefTemp->set_parameters( { {{ 15.0 }} } );
    tPropRefTemp->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropRefTemp->set_val_function( tTEMPFIValFunction_NSPBULK );
    tPropRefTemp->set_dof_derivative_functions( { tTEMPFIDerFunction_NSPBULK } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMMasterIncFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }} );
    tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
    tCMMasterIncFluid->set_property( tPropDensity, "Density" );
    tCMMasterIncFluid->set_space_dim( 3 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPIncFlow = tSPFactory.create_SP( fem::Stabilization_Type::INCOMPRESSIBLE_FLOW );
    tSPIncFlow->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPIncFlow->set_parameters( { {{ 36.0 }} } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::P } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
    tIWG->set_property( tPropDensity, "Density" );
    tIWG->set_property( tPropGravity, "Gravity" );
    tIWG->set_property( tPropThermalExp, "ThermalExpansion" );
    tIWG->set_property( tPropRefTemp, "ReferenceTemp" );
    tIWG->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.75 }, { 0.1 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
                              { 1.0, 0.0, 0.0 },
                              { 1.0, 1.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 },
                              { 1.0, 0.0, 1.0 },
                              { 1.0, 1.0, 1.0 },
                              { 0.0, 1.0, 1.0 }};

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }, { 1.0 }};

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    // set the evaluation point
    tGI.set_space_time( tParamPoint );

    // field interpolators
    //------------------------------------------------------------------------------
    //create a space time interpolation rule
    Interpolation_Rule tFIRule ( mtk::Geometry_Type::HEX,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR,
                                 Interpolation_Type::LAGRANGE,
                                 mtk::Interpolation_Order::LINEAR );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 16, 3 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    arma::Mat< double > tMatrixP;
    tMatrixP.randu( 16, 1 );
    Matrix< DDRMat > tDOFHatP;
    tDOFHatP.matrix_data() = 10.0 * tMatrixP;

    arma::Mat< double > tMatrixT;
    tMatrixT.randu( 16, 1 );
    Matrix< DDRMat > tDOFHatT;
    tDOFHatT.matrix_data() = 10.0 * tMatrixT;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 3 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 3, tFIRule, &tGI,{ { MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ { MSI::Dof_Type::P } } );
    tFIs( 2 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ { MSI::Dof_Type::TEMP } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatP );
    tFIs( 2 )->set_coeff( tDOFHatT );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );
    tFIs( 2 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 2;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 2;

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( 3 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 47 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 48, 63 } };
    tIWG->mSet->mResDofAssemblyMap( 2 ) = { { 64, 79 } };

    // set size and fill the set jacobian assembly map
    tIWG->mSet->mJacDofAssemblyMap.resize( 3 );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 47 }, { 48, 63 }, { 64, 79 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 47 }, { 48, 63 }, { 64, 79 } };
    tIWG->mSet->mJacDofAssemblyMap( 2 ) = { { 0, 47 }, { 48, 63 }, { 64, 79 } };

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 80, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 80, 80, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P }, { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum PDV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager(&tFIManager);

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    tIWG->compute_residual( 1.0 );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tJacobians;
    Matrix< DDRMat > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                1.0,
                                                tJacobians,
                                                tJacobiansFD );

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/
