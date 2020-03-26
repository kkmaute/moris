/*
 * UT_FEM_IQI_Strain_Energy.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include <string>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private   public
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#include "cl_FEM_IQI.hpp"         //FEM/INT/src
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#undef protected
#undef private

#include "cl_FEM_Field_Interpolator.hpp"                   //FEM//INT//src
#include "cl_FEM_Property.hpp"                   //FEM//INT//src
#include "cl_FEM_CM_Factory.hpp"                   //FEM//INT//src
#include "cl_FEM_IQI_Factory.hpp"                   //FEM//INT//src
#include "cl_MTK_Enums.hpp" //MTK/src
#include "cl_FEM_Enums.hpp"                                //FEM//INT/src
#include "cl_MSI_Design_Variable_Interface.hpp"                   //FEM//INT//src
#include "FEM/MSI/test/MSI_Test_Proxy/cl_MSI_Design_Variable_Interface_Proxy.hpp"

#include "op_equal_equal.hpp"

void tConstValFunction_UTIQISTRAINENERGY
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 );
}


void tFIValDvFunction_UTIQISTRAINENERGY
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( GEN_DV::DENSITY0 )->val();
}

void tFIDerDvFunction_UTIQISTRAINENERGY
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
  moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
  moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix = aParameters( 0 ) * aFIManager->get_field_interpolators_for_type( GEN_DV::DENSITY0 )->N();
}

using namespace moris;
using namespace fem;


TEST_CASE( "IQI_Strain_Energy", "[moris],[fem],[IQI_Strain_Energy]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create the properties
    std::shared_ptr< fem::Property > tPropMasterEMod = std::make_shared< fem::Property > ();
    tPropMasterEMod->set_parameters( { {{ 1.0 }} } );
    tPropMasterEMod->set_dv_type_list( {{ GEN_DV::DENSITY0 }} );
    tPropMasterEMod->set_val_function( tFIValDvFunction_UTIQISTRAINENERGY );
    tPropMasterEMod->set_dv_derivative_functions( { tFIDerDvFunction_UTIQISTRAINENERGY } );

    std::shared_ptr< fem::Property > tPropMasterNu = std::make_shared< fem::Property > ();
    tPropMasterNu->set_parameters( { {{ 0.3 }} } );
    tPropMasterNu->set_val_function( tConstValFunction_UTIQISTRAINENERGY );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMMasterElastLinIso = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMMasterElastLinIso->set_dof_type_list( {{ MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }} );
    tCMMasterElastLinIso->set_property( tPropMasterEMod, "YoungsModulus" );
    tCMMasterElastLinIso->set_property( tPropMasterNu, "PoissonRatio" );
    tCMMasterElastLinIso->set_space_dim( 3 );

    // define the IWGs
    fem::IQI_Factory tIQIFactory;

    std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::STRAIN_ENERGY );
    tIQI->set_constitutive_model( tCMMasterElastLinIso, "Elast", mtk::Master_Slave::MASTER );
    tIQI->set_IQI_phase_type( Phase_Type::PHASE0 );

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
    arma::Mat< double > tXMatrix;
    tXMatrix.randu( 8, 3 );
    Matrix< DDRMat > tXHat;
    tXHat.matrix_data() = 10.0 * tXMatrix;

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
    tMatrix.randu( 8, 3 );
    Matrix< DDRMat > tDOFHat;
    tDOFHat.matrix_data() = 10.0 * tMatrix;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tFIs( 1 );

    // create the field interpolator
    tFIs( 0 ) = new Field_Interpolator( 3, tFIRule, &tGI, { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } );

    // set the coefficients uHat
    tFIs( 0 )->set_coeff( tDOFHat );

    //set the evaluation point xi, tau
    tFIs( 0 )->set_space_time( tParamPoint );

    // create random dv coefficients
    arma::Mat< double > tMatrixDv;
    tMatrixDv.randu( 8, 1 );
    Matrix< DDRMat > tDvHat;
    tDvHat.matrix_data() = 10.0 * tMatrixDv;

    // create a cell of field interpolators for IWG
    Cell< Field_Interpolator* > tDvFIs( 1 );

    // create the field interpolator
    tDvFIs( 0 ) = new Field_Interpolator( 1, tFIRule, &tGI, { GEN_DV::DENSITY0 } );

    // set the coefficients uHat
    tDvFIs( 0 )->set_coeff( tDvHat );

    //set the evaluation point xi, tau
    tDvFIs( 0 )->set_space_time( tParamPoint );

    // create a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();

    // create a GEN/MSI interface
    MSI::Design_Variable_Interface * tGENMSIInterface = new MSI::Design_Variable_Interface_Proxy();
    tSet->set_dv_interface( tGENMSIInterface );

    // set fem set pointer for IQI
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set mUniqueDofTypeList
    tIQI->mSet->mUniqueDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    // set size for the set mUniqueDvTypeList
    tIQI->mSet->mUniqueDvTypeList.resize( 4, GEN_DV::END_ENUM );

    // set size and populate the set dof type map
    tIQI->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIQI->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set dof type map
    tIQI->mSet->mUniqueDvTypeMap.set_size( static_cast< int >( GEN_DV::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mUniqueDvTypeMap( static_cast< int >( GEN_DV::DENSITY0 ) ) = 0;

    // set size and populate the set master dof type map
    tIQI->mSet->mMasterDvTypeMap.set_size( static_cast< int >( GEN_DV::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mMasterDvTypeMap( static_cast< int >( GEN_DV::DENSITY0 ) ) = 0;

    // set size and populate residual assembly map
    tIQI->mSet->mResDofAssemblyMap.resize( 1 );
    tIQI->mSet->mResDofAssemblyMap( 0 ) = { { 0, 23 } };

    // set size and populate jacobian assembly map
    tIQI->mSet->mJacDofAssemblyMap.resize( 1 );
    tIQI->mSet->mJacDofAssemblyMap( 0 ) = { { 0, 23 } };

    // set size and fill the set dv assembly map
    tIQI->mSet->mDvAssemblyMap.resize( 1 );
    tIQI->mSet->mDvAssemblyMap( 0 ) = { { 0, 7 } };

    // set size and init the set residual and jacobian
    tIQI->mSet->mResidual.resize( 1 );
    tIQI->mSet->mResidual( 0 ).set_size( 24, 1, 0.0 );

    // set size and init the set residual and jacobian
    tIQI->mSet->mQI.resize( 1 );
    tIQI->mSet->mQI( 0 ).set_size( 1, 1, 0.0 );

    // populate the requested master dof type
    tIQI->mRequestedMasterGlobalDofTypes = {{ MSI::Dof_Type::UX }};

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager tFIManager( tDummy, tSet );

    // populate the field interpolator manager
    tFIManager.mFI   = tFIs;
    tFIManager.mDvFI = tDvFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set
    tIQI->mSet->mMasterFIManager = &tFIManager;

    // set IWG field interpolator manager
    tIQI->set_field_interpolator_manager( &tFIManager );

    // check evaluation of the quantity of interest
    //------------------------------------------------------------------------------
    // evaluate the quantity of interest
    Matrix< DDRMat > tQI;
    tIQI->compute_QI( tQI );
//    print( tQI, "tQI" );
    tIQI->compute_QI( 1.0 );
//    print( tIQI->mSet->get_QI()( 0 ), "QI" );

    // check evaluation of the derivative of the quantity of interest wrt to dof
    //------------------------------------------------------------------------------
    // evaluate the quantity of interest derivatives wrt to dof
    Matrix< DDRMat > tdQIdu;
    Matrix< DDRMat > tdQIduFD;
    bool tCheckdQIdu = tIQI->check_dQIdu_FD( 1.0,
                                             tPerturbation,
                                             tEpsilon,
                                             tdQIdu,
                                             tdQIduFD );
//    // print for debug
//    print( tdQIdu,   "tdQIdu" );
//    print( tdQIduFD, "tdQIduFD" );

    // require check is true
    REQUIRE( tCheckdQIdu );

    // check evaluation of the derivative of the quantity of interest wrt to dv
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tdQIdpMatFD;
    tIQI->compute_dQIdp_FD_material( 1.0, tPerturbation, tdQIdpMatFD );
//    print( tdQIdpMatFD, "tdQIdpMatFD" );

    // assume that every x, y, z coords are active pdv
    moris::Cell< Matrix< DDSMat > > tIsActive( 8 );
    tIsActive( 0 ) = { {1}, {1}, {1} };
    tIsActive( 1 ) = { {1}, {1}, {1} };
    tIsActive( 2 ) = { {1}, {1}, {1} };
    tIsActive( 3 ) = { {1}, {1}, {1} };
    tIsActive( 4 ) = { {1}, {1}, {1} };
    tIsActive( 5 ) = { {1}, {1}, {1} };
    tIsActive( 6 ) = { {1}, {1}, {1} };
    tIsActive( 7 ) = { {1}, {1}, {1} };

    Matrix< DDRMat > tdQIdpGeoFD;
    tIQI->compute_dQIdp_FD_geometry( 1.0, tPerturbation, tIsActive, tdQIdpGeoFD );
//    print( tdQIdpGeoFD, "tdQIdpGeoFD" );

}/*END_TEST_CASE*/



