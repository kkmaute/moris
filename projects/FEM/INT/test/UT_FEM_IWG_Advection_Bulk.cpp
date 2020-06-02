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

void tConstValFunction_Advection
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

void tFIValFunction_Advection
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val();
}

void tFIDerFunction_Advection
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Advection_Bulk_2D", "[IWG_Advection_Bulk_2D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropHeatCapacity->set_val_function( tFIValFunction_Advection );
    tPropHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropDensity->set_val_function( tFIValFunction_Advection );
    tPropDensity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropConductivity->set_val_function( tFIValFunction_Advection );
    tPropConductivity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity, "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "Heat_Capacity" );
    tCMDiffusion->set_space_dim( 2 );


    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG
    = tSPFactory.create_SP( fem::Stabilization_Type::SUPG_ADVECTION );
    tSPSUPG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }} );
    tSPSUPG->set_property( tPropConductivity, "Conductivity", mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );

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

    arma::Mat< double > tMatrixT;
    tMatrixT.randu( 8, 1 );
    Matrix< DDRMat > tDOFHatT;
    tDOFHatT.matrix_data() = 10.0 * tMatrixT;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI,{ {MSI::Dof_Type::VX, MSI::Dof_Type::VY } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ {MSI::Dof_Type::TEMP } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatT );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( 2 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 15 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 16, 23 } };

    // set size and fill the set jacobian assembly map
    tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 15 }, { 16, 23 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 15 }, { 16, 23 } };

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 24, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 24, 24, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
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
    Matrix< DDRMat > tJacobian;
    Matrix< DDRMat > tJacobianFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                1.0,
                                                tJacobian,
                                                tJacobianFD );

//    // print for debug
//    print( tJacobian({ 0, 7 }, { 0, 15 }), "tJacobian_TV" );
//    print( tJacobianFD({ 0, 7 }, { 0, 15 }), "tJacobianFD_TV" );
//    print( tJacobian({ 0, 7 }, { 16, 23 }), "tJacobian_TT" );
//    print( tJacobianFD({ 0, 7 }, { 16, 23 }), "tJacobianFD_TT" );

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/

TEST_CASE( "IWG_Advection_Bulk_3D", "[IWG_Advection_Bulk_3D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropHeatCapacity->set_val_function( tFIValFunction_Advection );
    tPropHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 10.0 }} } );
    tPropDensity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropDensity->set_val_function( tFIValFunction_Advection );
    tPropDensity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropConductivity = std::make_shared< fem::Property >();
    tPropConductivity->set_parameters( { {{ 5.0 }} } );
    tPropConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropConductivity->set_val_function( tFIValFunction_Advection );
    tPropConductivity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMDiffusion = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMDiffusion->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMDiffusion->set_property( tPropConductivity, "Conductivity" );
    tCMDiffusion->set_property( tPropDensity, "Density" );
    tCMDiffusion->set_property( tPropHeatCapacity, "Heat_Capacity" );
    tCMDiffusion->set_space_dim( 3 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPSUPG
    = tSPFactory.create_SP( fem::Stabilization_Type::SUPG_ADVECTION );
    tSPSUPG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }} );
    tSPSUPG->set_property( tPropConductivity, "Conductivity", mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMDiffusion, "Diffusion" );
    tIWG->set_stabilization_parameter( tSPSUPG, "SUPG" );

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

    arma::Mat< double > tMatrixT;
    tMatrixT.randu( 16, 1 );
    Matrix< DDRMat > tDOFHatT;
    tDOFHatT.matrix_data() = 10.0 * tMatrixT;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 3, tFIRule, &tGI,{ {MSI::Dof_Type::VX, MSI::Dof_Type::VY, MSI::Dof_Type::VZ } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ {MSI::Dof_Type::TEMP } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatT );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );
    tFIs( 1 )->set_space_time( tParamPoint );

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( 2 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 47 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 48, 63 } };

    // set size and fill the set jacobian assembly map
    tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 47 }, { 48, 63 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 47 }, { 48, 63 } };

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 64, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 64, 64, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum PDV_Type > > tDummyDv;
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
    Matrix< DDRMat > tJacobian;
    Matrix< DDRMat > tJacobianFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                1.0,
                                                tJacobian,
                                                tJacobianFD );

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/
