#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                   //FEM//INT//src

#include "op_equal_equal.hpp"

moris::Matrix< moris::DDRMat > tConstValFunction_UTIWGDIFFDIR( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                               moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                               moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                               moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tGeoValFunction_UTIWGDIFFDIR( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                             moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aGeometryInterpolator->valx()( 0 );
}

moris::Matrix< moris::DDRMat > tFIValFunction_UTIWGDIFFDIR( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                            moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction_UTIWGDIFFDIR( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                            moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                            moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Diff_Dirichlet_Const_Prop", "[moris],[fem],[IWG_Diff_Dirichlet_Const_Prop]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
    tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
    tPropMasterConductivity->set_val_function( tConstValFunction_UTIWGDIFFDIR );

    std::shared_ptr< fem::Property > tPropMasterDirichlet = std::make_shared< fem::Property > ();
    tPropMasterDirichlet->set_parameters( { {{ 1.0 }} } );
    tPropMasterDirichlet->set_val_function( tConstValFunction_UTIWGDIFFDIR );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_properties( { tPropMasterConductivity } );
    tCMMasterDiffLinIso->set_space_dim( 3 );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_models( { tCMMasterDiffLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

    // set the normal
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG->set_normal( tNormal );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

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
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 1 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG field interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian  by FD
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

}/*END_TEST_CASE*/

TEST_CASE( "IWG_Diff_Dirichlet_Geo_Prop", "[moris],[fem],[IWG_Diff_Dirichlet_Geo_Prop]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
    tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
    tPropMasterConductivity->set_val_function( tGeoValFunction_UTIWGDIFFDIR );

    std::shared_ptr< fem::Property > tPropMasterDirichlet = std::make_shared< fem::Property > ();
    tPropMasterDirichlet->set_parameters( { {{ 1.0 }} } );
    tPropMasterDirichlet->set_val_function( tGeoValFunction_UTIWGDIFFDIR );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_properties( { tPropMasterConductivity } );
    tCMMasterDiffLinIso->set_space_dim( 3 );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_models( { tCMMasterDiffLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

    // set the normal
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG->set_normal( tNormal );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

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
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 1 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG geometry interpolators
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian  by FD
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

TEST_CASE( "IWG_Diff_Dirichlet_Dof_Prop", "[moris],[fem],[IWG_Diff_Dirichlet_Dof_Prop]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
    tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
    tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropMasterConductivity->set_val_function( tFIValFunction_UTIWGDIFFDIR );
    tPropMasterConductivity->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFDIR } );

    std::shared_ptr< fem::Property > tPropMasterDirichlet = std::make_shared< fem::Property > ();
    tPropMasterDirichlet->set_parameters( { {{ 1.0 }} } );
    tPropMasterDirichlet->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropMasterDirichlet->set_val_function( tFIValFunction_UTIWGDIFFDIR );
    tPropMasterDirichlet->set_dof_derivative_functions( { tFIDerFunction_UTIWGDIFFDIR } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_properties( { tPropMasterConductivity } );
    tCMMasterDiffLinIso->set_space_dim( 3 );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_DIRICHLET );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_models( { tCMMasterDiffLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_properties( { tPropMasterDirichlet }, mtk::Master_Slave::MASTER );

    // set the normal
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG->set_normal( tNormal );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 0.35}, {-0.25}, { 0.75}, { 0.0 }};

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // create a space geometry interpolation rule
    Interpolation_Rule tGIRule( mtk::Geometry_Type::HEX,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR,
                                Interpolation_Type::LAGRANGE,
                                mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

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
                                 Interpolation_Type::CONSTANT,
                                 mtk::Interpolation_Order::CONSTANT );

    // create random coefficients
    arma::Mat< double > tMatrix;
    tMatrix.randu( 8, 1 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    // build global dof type list
    tIWG->build_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tFIs );

    // set IWG geometry interpolator
    tIWG->set_geometry_interpolator( &tGI );

    // check evaluation of the residual for IWG
    //------------------------------------------------------------------------------
    // evaluate the residual
    Cell< Matrix< DDRMat > > tResidual;
    tIWG->compute_residual( tResidual );

    // check evaluation of the jacobian  by FD
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

