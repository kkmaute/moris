#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src
                               //FEM//INT//src

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                               //FEM//INT//src
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"      //FEM//INT//src

#include "op_equal_equal.hpp"

void tConstValFunction_UTInterface
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
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
    tCMMasterStrucLinIso->set_model_type(fem::Model_Type::PLANE_STRESS);

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveStrucLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMSlaveStrucLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY }} );
    tCMSlaveStrucLinIso->set_property( tPropSlaveEMod, "YoungsModulus" );
    tCMSlaveStrucLinIso->set_property( tPropSlaveNu, "PoissonRatio" );
    tCMSlaveStrucLinIso->set_space_dim( 2 );
    tCMSlaveStrucLinIso->set_model_type(fem::Model_Type::PLANE_STRESS);

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

    MSI::Equation_Set * tSet = new fem::Set();

    tIWG->set_set_pointer(static_cast<fem::Set*>(tSet));

    tIWG->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;

    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >(MSI::Dof_Type::UX) ) = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >(MSI::Dof_Type::UX) ) = 0;

    tIWG->mSet->mResDofAssemblyMap.resize( 2 );
    tIWG->mSet->mJacDofAssemblyMap.resize( 2 );
    tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, 7 } };
    tIWG->mSet->mResDofAssemblyMap( 1 ) = { { 8, 15 } };
    tIWG->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 7 },{ 8, 15 } };
    tIWG->mSet->mJacDofAssemblyMap( 1 ) = { { 0, 7 },{ 8, 15 } };

    tIWG->mSet->mResidual.resize( 1 );
    tIWG->mSet->mResidual( 0 ).set_size( 16, 1 , 0.0 );
    tIWG->mSet->mJacobian.set_size( 16, 16, 0.0 );

    tIWG->mResidualDofTypeRequested = true;

    // build global dof type list
    tIWG->get_global_dof_type_list();

    tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::UX }};
    tIWG->mRequestedSlaveGlobalDofTypes  = {{ MSI::Dof_Type::UX }};

    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager tMasterFIManager( tDummy, tSet, mtk::Master_Slave::MASTER );
    Field_Interpolator_Manager tSlaveFIManager( tDummy, tSet, mtk::Master_Slave::SLAVE );

    // populate the field interpolator manager
    tMasterFIManager.mFI = tMasterFIs;
    tMasterFIManager.mIPGeometryInterpolator = &tGI;
    tMasterFIManager.mIGGeometryInterpolator = &tGI;
    tSlaveFIManager.mFI  = tSlaveFIs;
    tSlaveFIManager.mIPGeometryInterpolator = &tGI;
    tSlaveFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIWG->mSet->mMasterFIManager = &tMasterFIManager;
    tIWG->mSet->mSlaveFIManager  = &tSlaveFIManager;

    // set IWG field interpolator manager
    tIWG->set_field_interpolator_manager( &tMasterFIManager );
    tIWG->set_field_interpolator_manager( &tSlaveFIManager, mtk::Master_Slave::SLAVE );

    // check evaluation of the residual for IWG Helmholtz Bulk ?
    //------------------------------------------------------------------------------
    // evaluate the residual
    tIWG->compute_residual( 1.0 );

    // check evaluation of the jacobian  by FD
    //------------------------------------------------------------------------------
    // init the jacobian for IWG and FD evaluation
    Matrix< DDRMat > tJacobians;
    Matrix< DDRMat > tJacobiansFD;

    // check jacobian by FD
    bool tCheckJacobian = tIWG->check_jacobian_double( tPerturbation,
                                                       tEpsilon,
                                                       1.0,
                                                       tJacobians,
                                                       tJacobiansFD );

//        // print for debug
//        print( tJacobians,"tJacobians");
//        print( tJacobiansFD,"tJacobiansFD");

    // require check is true
    REQUIRE( tCheckJacobian );

    // clean up
    tMasterFIs.clear();
    tSlaveFIs.clear();

}/* END_TEST_CASE */




