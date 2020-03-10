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

void tConstValFunction_NSVBULK
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Incompressible_NS_Velocity_Bulk", "[IWG_Incompressible_NS_Velocity_Bulk]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropViscosity = std::make_shared< fem::Property >();
    tPropViscosity->set_parameters( { {{ 1.0 }} } );
    tPropViscosity->set_val_function( tConstValFunction_NSVBULK );

    std::shared_ptr< fem::Property > tPropDensity = std::make_shared< fem::Property >();
    tPropDensity->set_parameters( { {{ 1.0 }} } );
    tPropDensity->set_val_function( tConstValFunction_NSVBULK );

    std::shared_ptr< fem::Property > tPropGravity = std::make_shared< fem::Property >();
    tPropGravity->set_parameters( { {{ 10.0 }, { 10.0}} } );
    tPropGravity->set_val_function( tConstValFunction_NSVBULK );

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

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::INCOMPRESSIBLE_NS_VELOCITY_BULK );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::VX } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMMasterIncFluid, "IncompressibleFluid" );
    tIWG->set_property( tPropDensity, "Density" );
    tIWG->set_property( tPropGravity, "Gravity" );
    tIWG->set_stabilization_parameter( tSPIncFlow, "IncompressibleFlow" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0 }, { 0.25 }, { 0.1 }};

    // set the normal
    Matrix< DDRMat > tNormal = {{1.0},{0.0}};
    tIWG->set_normal( tNormal );

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
//    Matrix< DDRMat > tDOFHat( 8, 2, 0.0 );
//    Matrix< DDRMat > tDOFHat1( 8, 1, 1.0 );
//    Matrix< DDRMat > tDOFHat2( 8, 1, 2.0 );
//    tDOFHat.get_column( 0 ) = tDOFHat1.matrix_data();
//    tDOFHat.get_column( 1 ) = tDOFHat2.matrix_data();

    arma::Mat< double > tMatrixP;
    tMatrixP.randu( 8, 1 );
    Matrix< DDRMat > tDOFHatP;
    tDOFHatP.matrix_data() = 10.0 * tMatrixP;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI,{ {MSI::Dof_Type::VX, MSI::Dof_Type::VY } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ {MSI::Dof_Type::P } } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );
    tFIs( 1 )->set_coeff( tDOFHatP );

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
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::VX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

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
    tIWG->mSet->mJacobian.set_size( 16, 24, 0.0 );

    // set requested residual dof type flag to true
    tIWG->mResidualDofTypeRequested = true;

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = { { MSI::Dof_Type::VX }, { MSI::Dof_Type::P } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummyDof;
    moris::Cell< moris::Cell< enum GEN_DV > > tDummyDv;
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

    print(tJacobians({ 0, 15 }, { 0, 15 }),   "tJacobiansVV" );
    print(tJacobiansFD({ 0, 15 }, { 0, 15 }), "tJacobiansFDVV" );
    print(tJacobians({ 0, 15 }, { 16, 23 }),   "tJacobiansVP" );
    print(tJacobiansFD({ 0, 15 }, { 16, 23 }), "tJacobiansFDVP" );

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tFIs.clear();

}/*END_TEST_CASE*/
