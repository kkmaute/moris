#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                      //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                 //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                 //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                //FEM//INT//src

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
    return aParameters( 0 )( 0 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 )( 0 ) * aDofFI( 0 )->N();
}
moris::Matrix< moris::DDRMat > tMValFunction_STRUCDIRICHLET( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return {{ aParameters( 0 )( 0 ), 0.0 },
            { 0.0, aParameters( 0 )( 1 ) }};
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Struc_Dirichlet_Const_Prop", "[moris],[fem],[IWG_Struc_Dirichlet_Const_Prop]" )
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
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMMasterStrucLinIso->set_properties( { tPropMasterEMod, tPropMasterNu } );
    tCMMasterStrucLinIso->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_properties( { tPropMasterEMod }, mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::MASTER );
    tIWG->set_stabilization_parameters( { tSPDirichletNitsche } );
    tIWG->set_constitutive_models( { tCMMasterStrucLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

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

    // create a cell of field interpolators for IWG
    uint tNumFIs = tIWG->get_global_dof_type_list().size();
    Cell< Field_Interpolator* > tFIs( tNumFIs );

    for( uint iDOF = 0; iDOF < tNumFIs; iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_global_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               &tGI,
                                               tIWG->get_global_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG field interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Cell< Cell< Matrix< DDRMat > > > tJacobians;
    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                tJacobians,
                                                tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Geo_Prop", "[moris],[fem],[IWG_Struc_Dirichlet_Geo_Prop]" )
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
    tPropMasterEMod->set_val_function( tGeoValFunction_STRUCDIRICHLET );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property >();
    tPropMasterNu->set_parameters( {{{ 0.3 }}} );
    tPropMasterNu->set_val_function( tConstValFunction_STRUCDIRICHLET );

    std::shared_ptr< fem::Property > tPropMasterDirichlet = std::make_shared< fem::Property >();
    tPropMasterDirichlet->set_parameters( {{{ 0.0 }, { 0.0 }}} );
    tPropMasterDirichlet->set_val_function( tGeoValFunction_STRUCDIRICHLET );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMMasterStrucLinIso->set_properties( { tPropMasterEMod, tPropMasterNu } );
    tCMMasterStrucLinIso->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_properties( { tPropMasterEMod }, mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::MASTER );
    tIWG->set_stabilization_parameters( { tSPDirichletNitsche } );
    tIWG->set_constitutive_models( { tCMMasterStrucLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

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

    // create a cell of field interpolators for IWG
    uint tNumFIs = tIWG->get_global_dof_type_list().size();
    Cell< Field_Interpolator* > tFIs( tNumFIs );

    for( uint iDOF = 0; iDOF < tNumFIs; iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_global_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               &tGI,
                                               tIWG->get_global_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG field interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Cell< Cell< Matrix< DDRMat > > > tJacobians;
    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                tJacobians,
                                                tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Dof_Prop", "[moris],[fem],[IWG_Struc_Dirichlet_Dof_Prop]" )
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
    tPropMasterDirichlet->set_parameters( { {{ 1.0 }} } );
    tPropMasterDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tPropMasterDirichlet->set_val_function( tFIValFunction_STRUCDIRICHLET );
    tPropMasterDirichlet->set_dof_derivative_functions( { tFIDerFunction_STRUCDIRICHLET } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMMasterStrucLinIso->set_properties( { tPropMasterEMod, tPropMasterNu } );
    tCMMasterStrucLinIso->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_properties( { tPropMasterEMod }, mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::MASTER );
    tIWG->set_stabilization_parameters( { tSPDirichletNitsche } );
    tIWG->set_constitutive_models( { tCMMasterStrucLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

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

    // create a cell of field interpolators for IWG
    uint tNumFIs = tIWG->get_global_dof_type_list().size();
    Cell< Field_Interpolator* > tFIs( tNumFIs );

    for( uint iDOF = 0; iDOF < tNumFIs; iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_global_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               &tGI,
                                               tIWG->get_global_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG field interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Cell< Cell< Matrix< DDRMat > > > tJacobians;
    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                tJacobians,
                                                tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

}/* END_TEST_CASE */

TEST_CASE( "IWG_Struc_Dirichlet_Select", "[moris],[fem],[IWG_Struc_Dirichlet_Select]" )
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
    tPropMasterDirichlet->set_parameters( { {{ 1.0 }} } );
    tPropMasterDirichlet->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tPropMasterDirichlet->set_val_function( tFIValFunction_STRUCDIRICHLET );
    tPropMasterDirichlet->set_dof_derivative_functions( { tFIDerFunction_STRUCDIRICHLET } );

    std::shared_ptr< fem::Property > tPropMasterDirichlet2 = std::make_shared< fem::Property >();
    tPropMasterDirichlet2->set_parameters( { {{ 1.0, 0.0 }} } );
    tPropMasterDirichlet2->set_val_function( tMValFunction_STRUCDIRICHLET );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMMasterStrucLinIso->set_properties( { tPropMasterEMod, tPropMasterNu } );
    tCMMasterStrucLinIso->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPDirichletNitsche = tSPFactory.create_SP( fem::Stabilization_Type::DIRICHLET_NITSCHE );
    tSPDirichletNitsche->set_parameters( { {{ 1.0 }} } );
    tSPDirichletNitsche->set_properties( { tPropMasterEMod }, mtk::Master_Slave::MASTER );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::MASTER );
    tIWG->set_stabilization_parameters( { tSPDirichletNitsche } );
    tIWG->set_constitutive_models( { tCMMasterStrucLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet, tPropMasterDirichlet2 }, mtk::Master_Slave::MASTER );

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

    // create a cell of field interpolators for IWG
    uint tNumFIs = tIWG->get_global_dof_type_list().size();
    Cell< Field_Interpolator* > tFIs( tNumFIs );

    for( uint iDOF = 0; iDOF < tNumFIs; iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_global_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                               tFIRule,
                                               &tGI,
                                               tIWG->get_global_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tFIs( iDOF )->set_space_time( tParamPoint );
    }

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG field interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Cell< Cell< Matrix< DDRMat > > > tJacobians;
    Cell< Cell< Matrix< DDRMat > > > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian( tPerturbation,
                                                tEpsilon,
                                                tJacobians,
                                                tJacobiansFD );
    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    for( Field_Interpolator* tFI : tFIs )
    {
        delete tFI;
    }
    tFIs.clear();

}/* END_TEST_CASE */

