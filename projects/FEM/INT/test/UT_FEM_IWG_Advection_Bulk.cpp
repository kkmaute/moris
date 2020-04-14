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
    tPropHeatCapacity->set_val_function( tConstValFunction_Advection );
    tPropHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunction_Advection );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_property( tPropDensity, "Density" );
    tIWG->set_property( tPropHeatCapacity, "HeatCapacity" );

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
    tSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set master dof type map
    tSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and fill the set residual assembly map
    tSet->mResDofAssemblyMap.resize( 2 );
    tSet->mResDofAssemblyMap( 0 ) = { { 0, 15 } };
    tSet->mResDofAssemblyMap( 1 ) = { { 16, 23 } };

    // set size and fill the set jacobian assembly map
    tSet->mJacDofAssemblyMap.resize( 2 );
    tSet->mJacDofAssemblyMap( 0 ) = { { 0, 15 }, { 16, 23 } };
    tSet->mJacDofAssemblyMap( 1 ) = { { 0, 15 }, { 16, 23 } };

    // set size and init the set residual and jacobian
    tSet->mResidual.resize( 1 );
    tSet->mResidual( 0 ).set_size( 24, 1, 0.0 );
    tSet->mJacobian.set_size( 24, 24, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum GEN_DV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    reinterpret_cast< fem::Set* >( tSet )->mMasterFIManager = &tFIManager;

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

TEST_CASE( "IWG_Advection_Bulk_3D", "[IWG_Advection_Bulk_3D]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropHeatCapacity = std::make_shared< fem::Property >();
    tPropHeatCapacity->set_parameters( { {{ 2.0 }} } );
    tPropHeatCapacity->set_val_function( tConstValFunction_Advection );
    tPropHeatCapacity->set_dof_derivative_functions( { tFIDerFunction_Advection } );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunction_Advection );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::ADVECTION_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_property( tPropDensity, "Density" );
    tIWG->set_property( tPropHeatCapacity, "HeatCapacity" );

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
    tSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    reinterpret_cast< fem::Set* >( tSet )->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and populate the set master dof type map
    tSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 1;

    // set size and fill the set residual assembly map
    tSet->mResDofAssemblyMap.resize( 2 );
    tSet->mResDofAssemblyMap( 0 ) = { { 0, 47 } };
    tSet->mResDofAssemblyMap( 1 ) = { { 48, 63 } };

    // set size and fill the set jacobian assembly map
    tSet->mJacDofAssemblyMap.resize( 2 );
    tSet->mJacDofAssemblyMap( 0 ) = { { 0, 47 }, { 48, 63 } };
    tSet->mJacDofAssemblyMap( 1 ) = { { 0, 47 }, { 48, 63 } };

    // set size and init the set residual and jacobian
    tSet->mResidual.resize( 1 );
    tSet->mResidual( 0 ).set_size( 64, 1, 0.0 );
    tSet->mJacobian.set_size( 64, 64, 0.0 );

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::TEMP } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum GEN_DV > > tDummyDv;
    Field_Interpolator_Manager tFIManager( tDummyDof, tDummyDv, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    reinterpret_cast< fem::Set* >( tSet )->mMasterFIManager = &tFIManager;

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
