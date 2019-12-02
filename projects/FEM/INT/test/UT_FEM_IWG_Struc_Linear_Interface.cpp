#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src
#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                                //FEM//INT//src

#include "op_equal_equal.hpp"


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

TEST_CASE( "IWG_Struc_Linear_Interface", "[moris],[fem],[IWG_Struc_Linear_Interface]" )
{
    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property > ();
    tPropMasterEMod->set_parameters( {{{ 1000000.0 }}} );
    tPropMasterEMod->set_val_function( tConstValFunction_UTInterface );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property > ();
    tPropMasterNu->set_parameters( {{{ 0.0 }}} );
    tPropMasterNu->set_val_function( tConstValFunction_UTInterface );

    std::shared_ptr< fem::Property > tPropSlaveEMod = std::make_shared< fem::Property > ();
    tPropSlaveEMod->set_parameters( {{{ 1000000.0 }}} );
    tPropSlaveEMod->set_val_function( tConstValFunction_UTInterface );

    std::shared_ptr< fem::Property > tPropSlaveNu = std::make_shared< fem::Property > ();
    tPropSlaveNu->set_parameters( {{{ 0.0 }}} );
    tPropSlaveNu->set_val_function( tConstValFunction_UTInterface );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMMasterStrucLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterStrucLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterStrucLinIso->set_space_dim( 2 );

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMSlaveStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_space_dim( 2 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
    tSPNitscheInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
    tSPMasterWeightInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPMasterWeightInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
    tSPSlaveWeightInterface->set_property( tPropMasterEMod, "Material", mtk::Master_Slave::MASTER );
    tSPSlaveWeightInterface->set_property( tPropSlaveEMod, "Material", mtk::Master_Slave::SLAVE );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::STRUC_LINEAR_INTERFACE );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::UX, MSI::Dof_Type::UY } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }}, mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
    tIWG->set_stabilization_parameter( tSPMasterWeightInterface, "MasterWeightInterface" );
    tIWG->set_stabilization_parameter( tSPSlaveWeightInterface, "SlaveWeightInterface" );
    tIWG->set_constitutive_model( tCMMasterStrucLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveStrucLinIso, "ElastLinIso", mtk::Master_Slave::SLAVE );

    // set the normal
    Matrix< DDRMat > tNormal = {{1.0},{0.0},{0.0}};
    tIWG->set_normal( tNormal );

    // create evaluation point xi, tau
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tParamPoint = {{ 0.35 }, { 0.35 }, { 0.0 }};

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
                              { 0.0, 1.0 }};

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
    Cell< Field_Interpolator* > tMasterFIs( tIWG->get_dof_type_list().size() );

    for( uint iDOF = 0; iDOF < tIWG->get_dof_type_list().size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_dof_type_list()( iDOF ).size();

        // create the field interpolator
        tMasterFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                                     tFIRule,
                                                     &tGI,
                                                     tIWG->get_dof_type_list()( iDOF ) );

        // set the coefficients uHat
        tMasterFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tMasterFIs( iDOF )->set_space_time( tParamPoint );
    }

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tSlaveFIs( tIWG->get_dof_type_list( mtk::Master_Slave::SLAVE ).size() );

    for( uint iDOF = 0; iDOF < tIWG->get_dof_type_list( mtk::Master_Slave::SLAVE ).size(); iDOF++ )
    {
        // get the number of DOF
        uint tNumOfFields = tIWG->get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ).size();

        // create the field interpolator
        tSlaveFIs( iDOF ) = new Field_Interpolator( tNumOfFields,
                                                    tFIRule,
                                                    &tGI,
                                                    tIWG->get_dof_type_list( mtk::Master_Slave::SLAVE )( iDOF ) );

        // set the coefficients uHat
        tSlaveFIs( iDOF )->set_coeff( tDOFHat );

        //set the evaluation point xi, tau
        tSlaveFIs( iDOF )->set_space_time( tParamPoint );
    }

    // define an epsilon environment
    real tEpsilon = 1E-4;

    // define a perturbation relative size
    real tPerturbation = 1E-4;

    SECTION( "IWG_Spatial_Diffusion : check residual and jacobian with constant property" )
    {
        // build global dof type list
        tIWG->build_global_dof_type_list();

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

//        // print for debug
//        print( tJacobians( 0 )( 0 ),"tJacobians00");
//        print( tJacobiansFD( 0 )( 0 ),"tJacobiansFD00");
//
//        print( tJacobians( 0 )( 1 ),"tJacobians01");
//        print( tJacobiansFD( 0 )( 1 ),"tJacobiansFD01");
//
//        print( tJacobians( 1 )( 0 ),"tJacobians10");
//        print( tJacobiansFD( 1 )( 0 ),"tJacobiansFD10");
//
//        print( tJacobians( 1 )( 1 ),"tJacobians11");
//        print( tJacobiansFD( 1 )( 1 ),"tJacobiansFD11");

        // require check is true
        REQUIRE( tCheckJacobian );

    }/* END_SECTION */

    // clean up
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
