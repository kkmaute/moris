#include <string>
#include <catch.hpp>

#include "assert.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                               //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Pressure_Dirichlet.hpp"      //FEM//INT//src

#include "op_equal_equal.hpp"

moris::Matrix< moris::DDRMat > tConstValFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                                 moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tGeoValFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                               moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                               moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aGeometryInterpolator->valx()( 0 );
}

moris::Matrix< moris::DDRMat > tFIValFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) + aParameters( 1 )( 0, 0 ) * ( aParameters( 2 ) - aDofFI( 0 )->val() );
}

moris::Matrix< moris::DDRMat > tFIDerFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return -1.0 * aParameters( 1 )( 0, 0 ) * aDofFI( 0 )->N();
}
moris::Matrix< moris::DDRMat > tMValFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return {{ aParameters( 0 )( 0 ), 0.0 },
            { 0.0, aParameters( 0 )( 1 ) }};
}

moris::Matrix< moris::DDRMat > tMValFunction_STRUCDIRICHLET_3D( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return {{ aParameters( 0 )( 0 ), 0.0, 0.0 },
            { 0.0, aParameters( 0 )( 1 ), 0.0 },
            { 0.0, 0.0, aParameters( 0 )( 2 ) }};
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Struc_Dirichlet_Mixed_Pressure", "[IWG_Struc_Dirichlet_Mixed]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create a linear elasticity Dirichlet IWG
    //------------------------------------------------------------------------------

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property >();
    tPropMasterEMod->set_parameters( {{{ 10.0 }}} );
    tPropMasterEMod->set_val_function( tConstValFunction_STRUCDIRICHLET );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( {{{ 0.3 }}} );
    tPropMasterNu->set_val_function( tConstValFunction_STRUCDIRICHLET );

    std::shared_ptr< fem::Property > tPropMasterDirichlet = std::make_shared< fem::Property >();
    tPropMasterDirichlet->set_parameters( {{{ 0.0 }, { 0.0 }}} );
    tPropMasterDirichlet->set_val_function( tConstValFunction_STRUCDIRICHLET );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }, {MSI::Dof_Type::P}} );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_space_dim( 2 );
    tCMMasterStrucLinIso->set_model_type(fem::Model_Type::PLANE_STRESS);
    tCMMasterStrucLinIso->set_model_type(fem::Model_Type::DEVIATORIC);

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_PRESSURE_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::P } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }, {MSI::Dof_Type::P}} );
    tIWG->set_stabilization_parameter( tSPDirichletNitsche, "DirichletNitsche" );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso" );
    tIWG->set_property( tPropMasterDirichlet, "Dirichlet" );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 1.0}, {0.25}, { 0.0}};

    // set the normal
    Matrix< DDRMat > tNormal = {{1.0},{0.0}};
    tIWG->set_normal( tNormal );

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::QUAD,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::CONSTANT,
                                mtk::Interpolation_Order::CONSTANT );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

    // create space coeff xHat
    Matrix< DDRMat > tXHat = {{ 0.0, 0.0 },
                              { 1.0, 0.0 },
                              { 1.0, 1.0 },
                              { 0.0, 1.0}};

    // create time coeff tHat
    Matrix< DDRMat > tTHat = {{ 0.0 }};

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
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 4, 2 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    arma::Mat< double > tMatrixP;
    tMatrixP.randu( 4, 1 );
    Matrix< DDRMat > tDOFHatP;
    tDOFHatP.matrix_data() = 10.0 * tMatrixP;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 2 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI,{ {MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
    tFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI,{ {MSI::Dof_Type::P} } );

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
    tIWG->mSet->mEqnObjDofTypeList.resize( 200, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
    tIWG->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and populate the set master dof type map
    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::P ) ) = 1;

    // set size and fill the set residual assembly map
    tIWG->mSet->mResDofAssemblyMap.resize( 2 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 0 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 0, 3 } };

    // set size and fill the set jacobian assembly map
    tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 0 }, { 0, 0 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 7 }, {8, 11} };

    // set size and init the set residual and jacobian
    tIWG->mSet->mResidual.set_size( 4, 1, 0.0 );
    tIWG->mSet->mJacobian.set_size( 4, 12, 0.0 );

    // set requested residual dof type flag to true
    tIWG->mResidualDofTypeRequested = true;

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // populate the requested master dof type
    tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::UX }, {MSI::Dof_Type::P}};

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager tFIManager( tDummy, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager( &tFIManager );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    tIWG->compute_residual( 1.0 );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Cell< Cell< Matrix< DDRMat > > > tJacobians;
    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

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

}/* END_TEST_CASE */
