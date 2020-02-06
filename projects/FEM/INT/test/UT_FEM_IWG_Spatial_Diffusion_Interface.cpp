#include <string>
#include <catch.hpp>

#include "assert.hpp"

#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                     //FEM//INT/src


#include "op_equal_equal.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                        //FEM//INT//src
#include "cl_FEM_Property.hpp"                                  //FEM//INT//src
#include "cl_FEM_IWG_Factory.hpp"                                //FEM//INT//src
#include "cl_FEM_SP_Factory.hpp"                                //FEM//INT//src
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
    tPropMasterConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropMasterConductivity->set_val_function( tFIValFunction_UTInterface );
    tPropMasterConductivity->set_dof_derivative_functions( { tFIDerFunction_UTInterface } );

    std::shared_ptr< fem::Property > tPropSlaveConductivity = std::make_shared< fem::Property > ();
    tPropSlaveConductivity->set_parameters( { {{ 5.0 }} } );
    tPropSlaveConductivity->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tPropSlaveConductivity->set_val_function( tFIValFunction_UTInterface );
    tPropSlaveConductivity->set_dof_derivative_functions( { tFIDerFunction_UTInterface } );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMMasterDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMMasterDiffLinIso->set_property( tPropMasterConductivity, "Conductivity" );
    tCMMasterDiffLinIso->set_space_dim( 3 );

    std::shared_ptr< fem::Constitutive_Model > tCMSlaveDiffLinIso = tCMFactory.create_CM( fem::Constitutive_Type::DIFF_LIN_ISO );
    tCMSlaveDiffLinIso->set_dof_type_list( {{ MSI::Dof_Type::TEMP }} );
    tCMSlaveDiffLinIso->set_property( tPropSlaveConductivity, "Conductivity" );
    tCMSlaveDiffLinIso->set_space_dim( 3 );

    // define stabilization parameters
    fem::SP_Factory tSPFactory;

    std::shared_ptr< fem::Stabilization_Parameter > tSPNitscheInterface = tSPFactory.create_SP( fem::Stabilization_Type::NITSCHE_INTERFACE );
    tSPNitscheInterface->set_parameters( { {{ 1.0 }} } );
    tSPNitscheInterface->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
    tSPNitscheInterface->set_property( tPropSlaveConductivity, "Material", mtk::Master_Slave::SLAVE );

    std::shared_ptr< fem::Stabilization_Parameter > tSPMasterWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::MASTER_WEIGHT_INTERFACE );
    tSPMasterWeightInterface->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
    tSPMasterWeightInterface->set_property( tPropSlaveConductivity, "Material", mtk::Master_Slave::SLAVE );

    std::shared_ptr< fem::Stabilization_Parameter > tSPSlaveWeightInterface = tSPFactory.create_SP( fem::Stabilization_Type::SLAVE_WEIGHT_INTERFACE );
    tSPSlaveWeightInterface->set_property( tPropMasterConductivity, "Material", mtk::Master_Slave::MASTER );
    tSPSlaveWeightInterface->set_property( tPropSlaveConductivity, "Material", mtk::Master_Slave::SLAVE );

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG = tIWGFactory.create_IWG( fem::IWG_Type::SPATIALDIFF_INTERFACE );
    tIWG->set_residual_dof_type( { MSI::Dof_Type::TEMP } );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::MASTER );
    tIWG->set_dof_type_list( {{ MSI::Dof_Type::TEMP }}, mtk::Master_Slave::SLAVE );
    tIWG->set_stabilization_parameter( tSPNitscheInterface, "NitscheInterface" );
    tIWG->set_stabilization_parameter( tSPMasterWeightInterface, "MasterWeightInterface" );
    tIWG->set_stabilization_parameter( tSPSlaveWeightInterface, "SlaveWeightInterface" );
    tIWG->set_constitutive_model( tCMMasterDiffLinIso, "DiffLinIso", mtk::Master_Slave::MASTER );
    tIWG->set_constitutive_model( tCMSlaveDiffLinIso, "DiffLinIso", mtk::Master_Slave::SLAVE );

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

    tIWG->mSet->mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    tIWG->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mSlaveDofTypeMap .set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tIWG->mSet->mMasterDofTypeMap( static_cast< int >(MSI::Dof_Type::TEMP) ) = 0;
    tIWG->mSet->mSlaveDofTypeMap ( static_cast< int >(MSI::Dof_Type::TEMP) ) = 0;

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

    tIWG->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::TEMP }};
    tIWG->mRequestedSlaveGlobalDofTypes  = {{ MSI::Dof_Type::TEMP }};

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

//    // print for debug
//    print( tJacobians,  "tJacobians");
//    print( tJacobiansFD,"tJacobiansFD");

    // require check is true
    REQUIRE( tCheckJacobian );

    tMasterFIs.clear();
    tSlaveFIs.clear();

}/* END_TEST_CASE */
