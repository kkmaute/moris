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
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"

void tConstValFunction_UTVelocityDirichlet
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Dirichlet_Nitsche_2D", "[IWG_Incompressible_NS_Dirichlet_Nitsche_2D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
    tPropVelocity->set_parameters( { {{ 10.0 }, { 20.0 } } } );
    tPropVelocity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMMasterIncFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }} );
    tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
    tCMMasterIncFluid->set_property( tPropDensity, "Density" );
    tCMMasterIncFluid->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }}, mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWGVelocity = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWGVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
    tIWGVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tIWGVelocity->set_property( tPropVelocity, "Dirichlet" );
    tIWGVelocity->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
    tIWGVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

    std::shared_ptr< fem::IWG > tIWGPressure = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWGPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
    tIWGPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tIWGPressure->set_property( tPropVelocity, "Dirichlet" );
    tIWGPressure->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.1 }};

    // set the normal
    Matrix< DDRMat > tNormal = {{ 1.0 }, { 0.0 }};
    tIWGVelocity->set_normal( tNormal );
    tIWGPressure->set_normal( tNormal );

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
    Matrix< DDRMat > tTHat = {{ 0.0}, {1.0 }};

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

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI,{ { MSI::Dof_Type::VX, MSI::Dof_Type::VY } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ { MSI::Dof_Type::P } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatP );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWGVelocity->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tIWGPressure->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWGVelocity->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWGVelocity->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and populate the set master dof type map
    tIWGVelocity->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) )  = 1;

    // set size and fill the set residual assembly map
    tIWGVelocity->mSet->mResDofAssemblyMap.resize( 2 );
    tIWGVelocity->mSet->mResDofAssemblyMap( 0 ) = { { 0, 15 } };
    tIWGVelocity->mSet->mResDofAssemblyMap( 1 ) = { { 16, 23 } };

    // set size and fill the set jacobian assembly map
    tIWGVelocity->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWGVelocity->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 15 }, { 16, 23 } };
    tIWGVelocity->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 15 }, { 16, 23 } };

    // set size and init the set residual and jacobian
    tIWGVelocity->mSet->mResidual.resize( 1 );
    tIWGVelocity->mSet->mResidual( 0 ).set_size( 24, 1, 0.0 );
    tIWGVelocity->mSet->mJacobian.set_size( 24, 24, 0.0 );

    // build global dof type list
    tIWGVelocity->get_global_dof_type_list();
    tIWGPressure->get_global_dof_type_list();

    // populate the requested master dof type
    tIWGVelocity->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };
    tIWGPressure->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum GEN_DV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWGVelocity->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWGVelocity->set_field_interpolator_manager( &tFIManager );
    tIWGPressure->set_field_interpolator_manager( &tFIManager );

    // check evaluation of the residual/jacobian for IWG Velocity
    //------------------------------------------------------------------------------
    // reset residual and jacobian
    tIWGVelocity->mSet->mResidual( 0 ).fill( 0.0 );
    tIWGVelocity->mSet->mJacobian.fill( 0.0 );

    // evaluate the residual
    tIWGVelocity->compute_residual( 1.0 );

    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tVelocityJacobian;
    Matrix< DDRMat > tVelocityJacobianFD;

    // check jacobian by FD
    bool tCheckVelocityJacobian = tIWGVelocity->check_jacobian( tPerturbation,
                                                                tEpsilon,
                                                                1.0,
                                                                tVelocityJacobian,
                                                                tVelocityJacobianFD );

//    print( tVelocityJacobian({ 0, 15 }, { 0, 15 }),    "tJacobianVV" );
//    print( tVelocityJacobianFD({ 0, 15 }, { 0, 15 }),  "tJacobianFDVV" );
//    print( tVelocityJacobian({ 0, 15 }, { 16, 23 }),   "tJacobianVP" );
//    print( tVelocityJacobianFD({ 0, 15 }, { 16, 23 }), "tJacobianFDVP" );

    // require check is true
    REQUIRE( tCheckVelocityJacobian );

    // check evaluation of the residual/jacobian for IWG Pressure
    //------------------------------------------------------------------------------
    // reset residual and jacobian
    tIWGPressure->mSet->mResidual( 0 ).fill( 0.0 );
    tIWGPressure->mSet->mJacobian.fill( 0.0 );

    // evaluate the residual
    tIWGVelocity->compute_residual( 1.0 );

    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tPressureJacobian;
    Matrix< DDRMat > tPressureJacobianFD;

    // check jacobian by FD
    bool tCheckPressureJacobian = tIWGPressure->check_jacobian( tPerturbation,
                                                                tEpsilon,
                                                                1.0,
                                                                tPressureJacobian,
                                                                tPressureJacobianFD );

//    print( tPressureJacobian(   { 0, 7 }, { 0, 15 }  ), "tJacobianPV" );
//    print( tPressureJacobianFD( { 0, 7 }, { 0, 15 }  ), "tJacobianFDPV" );
//    print( tPressureJacobian(   { 0, 7 }, { 16, 23 } ), "tJacobianPP" );
//    print( tPressureJacobianFD( { 0, 7 }, { 16, 23 } ), "tJacobianFDPP" );

    // require check is true
    REQUIRE( tCheckPressureJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/

TEST_CASE( "IWG_Incompressible_NS_Dirichlet_Nitsche_3D", "[IWG_Incompressible_NS_Dirichlet_Nitsche_3D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    std::shared_ptr< fem::Property > tPropVelocity = std::make_shared< fem::Property >();
    tPropVelocity->set_parameters( { {{ 10.0 }, { 20.0 }, { 30.0 }} } );
    tPropVelocity->set_val_function( tConstValFunction_UTVelocityDirichlet );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterIncFluid = tCMFactory.create_CM( fem::Constitutive_Type::FLUID_INCOMPRESSIBLE );
    tCMMasterIncFluid->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }} );
    tCMMasterIncFluid->set_property( tPropViscosity, "Viscosity" );
    tCMMasterIncFluid->set_property( tPropDensity, "Density" );
    tCMMasterIncFluid->set_space_dim( 3 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitsche = tSPFactory.create_SP( fem::Stabilization_Type::VELOCITY_DIRICHLET_NITSCHE );
    tSPNitsche->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }}, mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropDensity, "Density", mtk::Master_Slave::MASTER );
    tSPNitsche->set_property( tPropViscosity, "Viscosity", mtk::Master_Slave::MASTER );
    tSPNitsche->set_parameters( { {{ 1.0 }}, {{ 1.0 }} } );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWGVelocity = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWGVelocity->set_residual_dof_type( { MSI::Dof_Type::VX } );
    tIWGVelocity->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tIWGVelocity->set_property( tPropVelocity, "Dirichlet" );
    tIWGVelocity->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
    tIWGVelocity->set_stabilization_parameter( tSPNitsche, "DirichletNitsche" );

    std::shared_ptr< fem::IWG > tIWGPressure = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE );
    tIWGPressure->set_residual_dof_type( { MSI::Dof_Type::P } );
    tIWGPressure->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tIWGPressure->set_property( tPropVelocity, "Dirichlet" );
    tIWGPressure->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.75 }, { 0.1 }};

    // set the normal
    Matrix< DDRMat > tNormal = {{ 1.0 }, { 0.0 }, { 0.0 } };
    tIWGVelocity->set_normal( tNormal );
    tIWGPressure->set_normal( tNormal );

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

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 3, tFIRule, &tGI,{ {MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ {MSI::Dof_Type::P } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatP );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWGVelocity->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );
    tIWGPressure->set_set_pointer( reinterpret_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWGVelocity->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWGVelocity->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGVelocity->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and populate the set master dof type map
    tIWGVelocity->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWGVelocity->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and fill the set residual assembly map
    tIWGVelocity->mSet->mResDofAssemblyMap.resize( 2 );
    tIWGVelocity->mSet->mResDofAssemblyMap( 0 ) = { { 0, 47 } };
    tIWGVelocity->mSet->mResDofAssemblyMap( 1 ) = { { 48, 63 } };

    // set size and fill the set jacobian assembly map
    tIWGVelocity->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWGVelocity->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 47 }, { 48, 63 } };
    tIWGVelocity->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 47 }, { 48, 63 } };

    // set size and init the set residual and jacobian
    tIWGVelocity->mSet->mResidual.resize( 1 );
    tIWGVelocity->mSet->mResidual( 0 ).set_size( 64, 1, 0.0 );
    tIWGVelocity->mSet->mJacobian.set_size( 64, 64, 0.0 );

    // build global dof type list
    tIWGVelocity->get_global_dof_type_list();
    tIWGPressure->get_global_dof_type_list();

    // populate the requested master dof type
    tIWGVelocity->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };
    tIWGPressure->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum GEN_DV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWGVelocity->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWGVelocity->set_field_interpolator_manager( &tFIManager );
    tIWGPressure->set_field_interpolator_manager( &tFIManager );

    // check evaluation of the residual/jacobian for IWG velocity
    //------------------------------------------------------------------------------
    // reset residual and jacobian
    tIWGVelocity->mSet->mResidual( 0 ).fill( 0.0 );
    tIWGVelocity->mSet->mJacobian.fill( 0.0 );

    // evaluate the residual
    tIWGVelocity->compute_residual( 1.0 );

    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tVelocityJacobian;
    Matrix< DDRMat > tVelocityJacobianFD;

    // check jacobian by FD
    bool tCheckVelocityJacobian = tIWGVelocity->check_jacobian( tPerturbation,
                                                                tEpsilon,
                                                                1.0,
                                                                tVelocityJacobian,
                                                                tVelocityJacobianFD );

    // require check is true
    REQUIRE( tCheckVelocityJacobian );

    // check evaluation of the residual/jacobian for IWG pressure
    //------------------------------------------------------------------------------
    // reset residual and jacobian
    tIWGPressure->mSet->mResidual( 0 ).fill( 0.0 );
    tIWGPressure->mSet->mJacobian.fill( 0.0 );

    // evaluate the residual
    tIWGPressure->compute_residual( 1.0 );

    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tPressureJacobian;
    Matrix< DDRMat > tPressureJacobianFD;

    // check jacobian by FD
    bool tCheckPressureJacobian = tIWGPressure->check_jacobian( tPerturbation,
                                                                tEpsilon,
                                                                1.0,
                                                                tPressureJacobian,
                                                                tPressureJacobianFD );

    // require check is true
    REQUIRE( tCheckPressureJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/
