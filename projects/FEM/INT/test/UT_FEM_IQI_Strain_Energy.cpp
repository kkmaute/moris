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

#include "op_equal_equal.hpp"

moris::Matrix< moris::DDRMat > tConstValFunction_UTIQISTRAINENERGY( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                    moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                                    moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                                    moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tFIValFunction_UTIQISTRAINENERGY( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                                 moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->val();
}

moris::Matrix< moris::DDRMat > tFIDerFunction_UTIQISTRAINENERGY( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDofFI,
                                                                 moris::Cell< moris::fem::Field_Interpolator* > & aDvFI,
                                                                 moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aDofFI( 0 )->N();
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
    tPropMasterEMod->set_val_function( tConstValFunction_UTIQISTRAINENERGY );

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
    tIQI->set_constitutive_model( tCMMasterElastLinIso, "ElastLinIso", mtk::Master_Slave::MASTER );

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

//    // create space coeff xHat
//    Matrix< DDRMat > tXHat = {{ 0.0, 0.0, 0.0 },
//                              { 1.0, 0.0, 0.0 },
//                              { 1.0, 1.0, 0.0 },
//                              { 0.0, 1.0, 0.0 },
//                              { 0.0, 0.0, 1.0 },
//                              { 1.0, 0.0, 1.0 },
//                              { 1.0, 1.0, 1.0 },
//                              { 0.0, 1.0, 1.0 }};

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

    // set a fem set pointer
    MSI::Equation_Set * tSet = new fem::Set();
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set IG GI in the set
    static_cast< fem::Set* >( tSet )->mMasterIGGeometryInterpolator = &tGI;

    // set size for the set EqnObjDofTypeList
    tIQI->mSet->mEqnObjDofTypeList.resize( 4, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIQI->mSet->mDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set master dof type map
    tIQI->mSet->mMasterDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIQI->mSet->mMasterDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate residual assembly map
    tIQI->mSet->mResDofAssemblyMap.resize( 1 );
    tIQI->mSet->mResDofAssemblyMap( 0 ) = { { 0, 23 } };

    // create a field interpolator manager
    moris::Cell< moris::Cell< enum MSI::Dof_Type > > tDummy;
    Field_Interpolator_Manager tFIManager( tDummy, tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFIs;

    // set IWG field interpolator manager
    tIQI->set_field_interpolator_manager( &tFIManager );

    // set IWG field interpolators
    tIQI->set_geometry_interpolator( &tGI );

    // check evaluation of the quantity of interest
    //------------------------------------------------------------------------------
    // evaluate the quantity of interest
    Matrix< DDRMat > tQI;
    tIQI->compute_QI( tQI );
    //print( tQI, "tQI" );

    // check evaluation of the derivative of the quantity of interest wrt to dof
    //------------------------------------------------------------------------------
    // evaluate the quantity of interest derivatives wrt to dof
    Matrix< DDRMat > tdQIdDof;
    Matrix< DDRMat > tdQIdDofFD;
    bool tCheckdQIdDof = tIQI->check_dQIdDof_FD( tPerturbation,
                                                 tEpsilon,
                                                 tdQIdDof,
                                                 tdQIdDofFD );

    // require check is true
    REQUIRE( tCheckdQIdDof );
    //print( tdQIdDof,   "tdQIdDof" );
    //print( tdQIdDofFD, "tdQIdDofFD" );

    // check evaluation of the derivative of the quantity of interest wrt to dv
    //------------------------------------------------------------------------------
    Matrix< DDRMat > tdQIdpMatFD;
    Matrix< DDRMat > tdQIdpGeoFD;
    tIQI->compute_dQIdDv_FD( tdQIdpMatFD,
                             tdQIdpGeoFD,
                             tPerturbation );
    //print( tdQIdpMatFD, "tdQIdpMatFD" );
    //print( tdQIdpGeoFD, "tdQIdpGeoFD" );


}/*END_TEST_CASE*/



