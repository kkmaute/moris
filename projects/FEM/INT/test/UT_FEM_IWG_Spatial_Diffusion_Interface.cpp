#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src


#include "op_equal_equal.hpp"

#define protected public
#define private   public
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                                //FEM//INT//src


moris::Matrix< moris::DDRMat > tConstValFunction_UTInterface( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                              moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tFIValFunction_UTInterface( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                           moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                           moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                           moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction_UTInterface( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                           moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                           moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                           moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "IWG_Diff_Interface", "[moris],[fem],[IWG_Diff_Interface]" )
{

    // create a spatial diffusion bulk IWG
    //------------------------------------------------------------------------------

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterConductivity = std::make_shared< fem::Property > ();
    tPropMasterConductivity->set_parameters( { {{ 1.0 }} } );
    tPropMasterConductivity->set_val_function( tConstValFunction_UTInterface );

    std::shared_ptr< fem::Property > tPropSlaveConductivity = std::make_shared< fem::Property > ();
    tPropSlaveConductivity->set_parameters( { {{ 1.0 }} } );
    tPropSlaveConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropSlaveConductivity->set_val_function( tFIValFunction_UTInterface );
    tPropSlaveConductivity->set_dof_derivative_functions( { tFIDerFunction_UTInterface } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_properties( { tPropMasterConductivity } );
    tCMMasterDiffLinIso->set_space_dim( 3 );

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMSlaveDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMSlaveDiffLinIso->set_properties( { tPropSlaveConductivity } );
    tCMSlaveDiffLinIso->set_space_dim( 3 );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );
    tIWG->set_constitutive_models( { tCMMasterDiffLinIso }, mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_models( { tCMSlaveDiffLinIso }, mtk::Master_Slave::SLAVE );

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
    Cell< Field_Interpolator* > tMasterFIs( 1 );

    // create the field interpolator
    tMasterFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tMasterFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tMasterFIs( 0 )->set_space_time( tParamPoint );

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tSlaveFIs( 1 );

    // create the field interpolator
    tSlaveFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { MSI::Dof_Type::TEMP } );

    // set the coefficients uHat
    tSlaveFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tSlaveFIs( 0 )->set_space_time( tParamPoint );

    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define aperturbation relative size
    real tPerturbation = 1E-6;

    MSI::Equation_Set * tSet = new fem::Set();

    tIWG->set_set_pointer(static_cast<fem::Set*>(tSet));

    tIWG->mSet->mEqnObjDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    tIWG->mSet->mDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 0;
    tIWG->mSet->mDofTypeMap( static_cast< int >(MSI::Dof_Type::VX) ) = 1;
    tIWG->mSet->mDofTypeMap( static_cast< int >(MSI::Dof_Type::LS1) ) = 2;
    tIWG->mSet->mDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 3;

    // build global dof type list
    tIWG->get_global_dof_type_list();

    // set IWG field interpolators
    tIWG->set_dof_field_interpolators( tMasterFIs );
    tIWG->set_dof_field_interpolators( tSlaveFIs, mtk::Master_Slave::SLAVE );

    // set IWG geometry interpolator
    tIWG->set_geometry_interpolator( &tGI );
    tIWG->set_geometry_interpolator( &tGI, mtk::Master_Slave::SLAVE );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
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
    bool tCheckJacobian = tIWG->check_jacobian_double( tPerturbation,
                                                       tEpsilon,
                                                       tJacobians,
                                                       tJacobiansFD );

//    // print for debug
//    print( tJacobians( 0 )( 0 ),"tJacobians00");
//    print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD00");
//
//    print( tJacobians( 0 )( 1 ),"tJacobians01");
//    print( tJacobiansFD( 0 )( 1 ),"tJacobiansFD01");
//
//    print( tJacobians( 1 )( 0 ),"tJacobians10");
//    print( tJacobiansFD( 1 )( 0 ),"tJacobiansFD10");
//
//    print( tJacobians( 1 )( 1 ),"tJacobians11");
//    print( tJacobiansFD( 1 )( 1 ),"tJacobiansFD11");

    // require check is true
    REQUIRE( tCheckJacobian );

    for( Field_Interpolator* tFI : tMasterFIs )
    {
        delete tFI;
    }
    tMasterFIs.clear();

    for( Field_Interpolator* tFI : tSlaveFIs )
    {
        delete tFI;
    }
    tSlaveFIs.clear();

}/* END_TEST_CASE */
