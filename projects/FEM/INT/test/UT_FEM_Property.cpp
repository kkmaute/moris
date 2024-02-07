/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Property.cpp
 *
 */

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

#define protected public
#define private   public
#include "cl_FEM_Set.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"                   //FEM//INT//src
#undef protected
#undef private

void tValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  { { 1.0 } };
}

void tDerFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    moris::Matrix< moris::DDRMat > tPropertyDer( 1, 1, 2.0);
    aPropMatrix =  tPropertyDer;
}

void tValFunction2
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 0 )
                         + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val()
                         + aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::UX )->val();
}

void tDerFunction2
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}

void tDerFunction3
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::UX )->N();
}

void tConstValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 0 );
}

void tGeoValFunction
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 0 ) * aFIManager->get_IP_geometry_interpolator()->valx();
}

void tValFunction3
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 0 )
                         + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->val()
                         + aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::UX )->val()
                         + aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::LS1 )->val()
                         + aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::LS2 )->val();
}
void tDerFunction3_TEMP
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::TEMP )->N();
}
void tDerFunction3_UX
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::MSI::Dof_Type::UX )->N();
}
void tDerFunction3_LS1
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 1 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::LS1 )->N();
}
void tDerFunction3_LS2
( moris::Matrix< moris::DDRMat >                 & aPropMatrix,
        moris::Vector< moris::Matrix< moris::DDRMat > >  & aParameters,
        moris::fem::Field_Interpolator_Manager         * aFIManager )
{
    aPropMatrix =  aParameters( 2 ) * aFIManager->get_field_interpolators_for_type( moris::gen::PDV_Type::LS2 )->N();
}

using namespace moris;
using namespace fem;

TEST_CASE( "Property", "[moris],[fem],[Property]" )
{
    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    // set coeffs
    Vector< Matrix< DDRMat > > tCoeff;

    // create a property object
    fem::Property tProperty;

    // set dof types
    Vector< Vector< MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }};
    tProperty.set_dof_type_list( tDofTypes );

    //set dv types
    Vector< Vector< gen::PDV_Type > > tDvTypes;
    tProperty.set_dv_type_list( tDvTypes );

    // set parameter
    tProperty.set_parameters( tCoeff );

    //set value function
    tProperty.set_val_function( tValFunction );

    // set dof derivative functions
    Vector< PropertyFunc > tDofDerFunc = { tDerFunction };
    tProperty.set_dof_derivative_functions( tDofDerFunc );

    // set dv derivative function
    Vector< PropertyFunc > tDvDerFunc;
    tProperty.set_dv_derivative_functions( tDvDerFunc );

    //check dof dependencies
    CHECK( equal_to( static_cast< uint >( tProperty.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

    // set field interpolators
    Vector< Field_Interpolator* > tFieldInterpolator( 1 );
    tFieldInterpolator( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP });

    // create a field interpolator manager
    fem::Set tSet; // dummy set
    tSet.mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;

    Field_Interpolator_Manager tFIManager( moris::Vector< moris::Vector< enum MSI::Dof_Type > >( 0 ), &tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFieldInterpolator;
    tFIManager.mIPGeometryInterpolator = &tGeomInterpolator;
    tProperty.set_field_interpolator_manager( &tFIManager );

    // check the property value
    Matrix< DDRMat > tPropertyValue = tProperty.val();
    CHECK( equal_to( tPropertyValue( 0, 0 ), 1.0 ) );

    // check that property depends on TEMP
    REQUIRE( tProperty.check_dof_dependency( { MSI::Dof_Type::TEMP } ) );

    // check that property does not depend on UX
    REQUIRE( !tProperty.check_dof_dependency( { MSI::Dof_Type::UX } ) );

    // check the property derivative wrt to TEMP (in dependencies)
    Matrix< DDRMat > tPropertyDerivative = tProperty.dPropdDOF( { MSI::Dof_Type::TEMP } );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 2.0 ) );

    // clean up
    tFieldInterpolator.clear();

}/* TEST CASE */

TEST_CASE( "Property_with_dependency", "[moris],[fem],[Property_with_dependency]" )
{
    //create a quad4 space element
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;
    //create a line time element
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat );

    // set the property coefficients
    Vector< Matrix< DDRMat > > tCoeff( 3 );
    tCoeff( 0 ) = {{ 1.0 }};
    tCoeff( 1 ) = {{ 2.0 }};
    tCoeff( 2 ) = {{ 3.0 }};

    // create a property object
    fem::Property tProperty;

    // set dof types
    Vector< Vector< MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }};
    tProperty.set_dof_type_list( tDofTypes );

    // set dv types
    Vector< Vector< gen::PDV_Type > > tDvTypes;
    tProperty.set_dv_type_list( tDvTypes );

    // set parameters
    tProperty.set_parameters( tCoeff );

    // set the value function
    tProperty.set_val_function( tValFunction2 );

    // set dof derivative functions
    Vector< PropertyFunc > tDofDerFunc = { tDerFunction2, tDerFunction3 };
    tProperty.set_dof_derivative_functions( tDofDerFunc );

    // set dv derivative functions
    Vector< PropertyFunc > tDvDerFunc;
    tProperty.set_dv_derivative_functions( tDvDerFunc );

    // create an interpolation rule
    mtk::Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );
    // create a TEMP field interpolator
    uint tNumberOfFields = 1;
    Vector< Field_Interpolator* > tFieldInterpolator( 2, nullptr );
    tFieldInterpolator( 0 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { MSI::Dof_Type::TEMP } );
    tFieldInterpolator( 1 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { MSI::Dof_Type::UX } );

    // set coefficients for field interpolators
    Matrix< DDRMat > tUHat0( 8, 1, 2.0 );
    Matrix< DDRMat > tUHat1( 8, 1, 4.0 );
    tFieldInterpolator( 0 )->set_coeff( tUHat0 );
    tFieldInterpolator( 1 )->set_coeff( tUHat1 );

    // set evaluation point for field interpolators
    Matrix< DDRMat > tParamPoint = {{ 0.0 }, { 0.0 }, { 0.0 }};
    tFieldInterpolator( 0 )->set_space_time( tParamPoint );
    tFieldInterpolator( 1 )->set_space_time( tParamPoint );

    // create a field interpolator manager
    fem::Set tSet; // dummy set
    tSet.mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 1;

    Field_Interpolator_Manager tFIManager( moris::Vector< moris::Vector< enum MSI::Dof_Type > >( 0 ), &tSet );

    // populate the field interpolator manager
    tFIManager.mFI = tFieldInterpolator;
    tFIManager.mIPGeometryInterpolator = &tGeomInterpolator;
    tProperty.set_field_interpolator_manager( &tFIManager );

    //evaluate the property
    Matrix< DDRMat > tPropertyValue = tProperty.val();
    CHECK( equal_to( tPropertyValue( 0, 0 ), 17.0 ) );

    // check that property depends on TEMP
    REQUIRE( tProperty.check_dof_dependency( { MSI::Dof_Type::TEMP } ) );

    // check that property depends on UX
    REQUIRE( tProperty.check_dof_dependency( { MSI::Dof_Type::UX } ) );

    // check that property does not depend on LS1
    REQUIRE( !tProperty.check_dof_dependency( { MSI::Dof_Type::LS1 } ) );

    // evaluate the property derivative wrt to TEMP (in dependencies)
    Matrix< DDRMat > tPropertyDerivative = tProperty.dPropdDOF( {MSI::Dof_Type::TEMP} );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.25 ) );

    // evaluate the property derivative wrt to UX (in dependencies)
    tPropertyDerivative = tProperty.dPropdDOF( {MSI::Dof_Type::UX} );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.375 ) );

    // clean up
    tFieldInterpolator.clear();

}/* TEST_CASE */

TEST_CASE( "Property_with_dof_dv_dependency", "[moris],[fem],[Property_with_dof_dv_dependency]" )
{
    //create a quad4 space element
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 1.25;
    tXHat( 2, 0 ) = 4.5; tXHat( 2, 1 ) = 4.0;
    tXHat( 3, 0 ) = 1.0; tXHat( 3, 1 ) = 3.25;
    //create a line time element
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat );

    // set the property coefficients
    Vector< Matrix< DDRMat > > tCoeff( 3 );
    tCoeff( 0 ) = {{ 1.0 }};
    tCoeff( 1 ) = {{ 2.0 }};
    tCoeff( 2 ) = {{ 3.0 }};

    // create a property object
    fem::Property tProperty;

    // set dof types
    Vector< Vector< MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }};
    tProperty.set_dof_type_list( tDofTypes );

    // set dv types
    Vector< Vector< gen::PDV_Type > > tDvTypes = {{ gen::PDV_Type::LS1 }, { gen::PDV_Type::LS2 }};
    tProperty.set_dv_type_list( tDvTypes );

    // set parameters
    tProperty.set_parameters( tCoeff );

    // set value function
    tProperty.set_val_function( tValFunction3 );

    // set dof derivative functions
    Vector< PropertyFunc > tDofDerFunc = { tDerFunction3_TEMP, tDerFunction3_UX };
    tProperty.set_dof_derivative_functions( tDofDerFunc );

    // set dv derivative functions
    Vector< PropertyFunc > tDvDerFunc = { tDerFunction3_LS1, tDerFunction3_LS2 };
    tProperty.set_dv_derivative_functions( tDvDerFunc );

    // create an interpolation rule
    mtk::Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );
    // create a dof field interpolators
    uint tNumberOfFields = 1;
    Vector< Field_Interpolator* > tDofFI( 2, nullptr );
    tDofFI( 0 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { MSI::Dof_Type::TEMP } );
    tDofFI( 1 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { MSI::Dof_Type::UX } );

    // create a dv field interpolators
    Vector< Field_Interpolator* > tDvFI( 2, nullptr );
    tDvFI( 0 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { gen::PDV_Type::LS1 } );
    tDvFI( 1 ) = new Field_Interpolator ( tNumberOfFields,
            tInterpolationRule,
            &tGeomInterpolator,
            { gen::PDV_Type::LS2 } );

    // set coefficients for field interpolators
    Matrix< DDRMat > tUHat0( 8, 1, 2.0 );
    Matrix< DDRMat > tUHat1( 8, 1, 4.0 );
    tDofFI( 0 )->set_coeff( tUHat0 );
    tDofFI( 1 )->set_coeff( tUHat1 );
    Matrix< DDRMat > tUHat2( 8, 1, 2.0 );
    Matrix< DDRMat > tUHat3( 8, 1, 4.0 );
    tDvFI( 0 )->set_coeff( tUHat2 );
    tDvFI( 1 )->set_coeff( tUHat3 );

    // set evaluation point for field interpolators
    Matrix< DDRMat > tParamPoint = {{ 0.0 }, { 0.0 }, { 0.0 }};
    tDofFI( 0 )->set_space_time( tParamPoint );
    tDofFI( 1 )->set_space_time( tParamPoint );
    tDvFI( 0 )->set_space_time( tParamPoint );
    tDvFI( 1 )->set_space_time( tParamPoint );

    // create a field interpolator manager
    fem::Set tSet; // dummy set
    tSet.mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )   = 1;
    tSet.mLeaderDvTypeMap.set_size( static_cast< int >(gen::PDV_Type::UNDEFINED) + 1, 1, -1 );
    tSet.mLeaderDvTypeMap( static_cast< int >( gen::PDV_Type::LS1 ) ) = 0;
    tSet.mLeaderDvTypeMap( static_cast< int >( gen::PDV_Type::LS2 ) )   = 1;
    Field_Interpolator_Manager tFIManager( moris::Vector< moris::Vector< enum MSI::Dof_Type > >( 0 ), &tSet );

    // populate the field interpolator manager
    tFIManager.mFI   = tDofFI;
    tFIManager.mDvFI = tDvFI;
    tFIManager.mIPGeometryInterpolator = &tGeomInterpolator;
    tProperty.set_field_interpolator_manager( &tFIManager );

    //evaluate the property
    Matrix< DDRMat > tPropertyValue = tProperty.val();
    CHECK( equal_to( tPropertyValue( 0, 0 ), 33.0 ) );

    // check that property depends on TEMP
    REQUIRE( tProperty.check_dof_dependency( { MSI::Dof_Type::TEMP } ) );

    // check that property depends on UX
    REQUIRE( tProperty.check_dof_dependency( { MSI::Dof_Type::UX } ) );

    // check that property does not depend on LS1
    REQUIRE( !tProperty.check_dof_dependency( { MSI::Dof_Type::LS1 } ) );

    // evaluate the property derivative wrt to TEMP (in dependencies)
    Matrix< DDRMat > tPropertyDerivative = tProperty.dPropdDOF( {MSI::Dof_Type::TEMP} );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.25 ) );

    // evaluate the property derivative wrt to UX (in dependencies)
    tPropertyDerivative = tProperty.dPropdDOF( {MSI::Dof_Type::UX} );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.375 ) );

    // check that property depends on dv LS1
    REQUIRE( tProperty.check_dv_dependency( { gen::PDV_Type::LS1 } ) );

    // check that property depends on dv LS2
    REQUIRE( tProperty.check_dv_dependency( { gen::PDV_Type::LS2 } ) );

    // check that property does not depend on dv UNDEFINED
    REQUIRE( !tProperty.check_dv_dependency( { gen::PDV_Type::UNDEFINED } ) );

    // evaluate the property derivative wrt to LS1 (in dependencies)
    tPropertyDerivative = tProperty.dPropdDV( {gen::PDV_Type::LS1} );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.25 ) );

    // evaluate the property derivative wrt to LS2 (in dependencies)
    tPropertyDerivative = tProperty.dPropdDV( {gen::PDV_Type::LS2 } );
    CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.375 ) );

    // clean up
    // delete field interpolators
    tDofFI.clear();
    tDvFI.clear();

}/* TEST_CASE */

TEST_CASE( "Property_geometry", "[moris],[fem],[Property_geometry]" )
{
    //create a quad4 space element
    Matrix< DDRMat > tXHat( 4, 2 );
    tXHat( 0, 0 ) = 0.0; tXHat( 0, 1 ) = 0.0;
    tXHat( 1, 0 ) = 3.0; tXHat( 1, 1 ) = 0.0;
    tXHat( 2, 0 ) = 3.0; tXHat( 2, 1 ) = 3.0;
    tXHat( 3, 0 ) = 0.0; tXHat( 3, 1 ) = 3.0;

    //create a line time element
    Matrix< DDRMat > tTHat( 2, 1 );
    tTHat( 0 ) = 0.0;
    tTHat( 1 ) = 5.0;

    //create a space geometry interpolation rule
    mtk::Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    //create a space and a time geometry interpolator
    Geometry_Interpolator tGeomInterpolator( tGeomInterpRule );

    //set the coefficients xHat, tHat
    tGeomInterpolator.set_coeff( tXHat, tTHat );

    //set the evaluation point tParamPoint
    Matrix< DDRMat > tParamPoint = {{ 0.0 }, { 0.0 }, { 0.0 }};
    tGeomInterpolator.set_space_time( tParamPoint );

    // create property coeffs
    Vector< Matrix< DDRMat > > tCoeff( 1 );
    tCoeff( 0 ) = {{ 1.0 }};

    // create a property object
    fem::Property tProperty;

    // set dof types
    Vector< Vector< MSI::Dof_Type > > tDofTypes = {{ MSI::Dof_Type::TEMP }};
    tProperty.set_dof_type_list( tDofTypes );

    // set dv types
    Vector< Vector< gen::PDV_Type > > tDvTypes;
    tProperty.set_dv_type_list( tDvTypes );

    // set parameters
    tProperty.set_parameters( tCoeff );

    // set value function
    tProperty.set_val_function( tGeoValFunction );

    // set dof derivative functions
    Vector< PropertyFunc > tDofDerFunc;
    tProperty.set_dof_derivative_functions( tDofDerFunc );

    // set dv derivative functions
    Vector< PropertyFunc > tDvDerFunc;
    tProperty.set_dv_derivative_functions( tDvDerFunc );

    //check dof dependencies
    CHECK( equal_to( static_cast< uint >( tProperty.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

    // set field interpolators
    Vector< Field_Interpolator* > tFieldInterpolator( 1 );
    tFieldInterpolator( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP });

    // create a field interpolator manager
    fem::Set tSet; // dummy set
    tSet.mLeaderDofTypeMap.set_size( static_cast< int >(MSI::Dof_Type::END_ENUM) + 1, 1, -1 );
    tSet.mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::TEMP ) ) = 0;
    Field_Interpolator_Manager tFIManager( moris::Vector< moris::Vector< enum MSI::Dof_Type > >( 0 ), &tSet );

    // populate the field interpolator manager
    tFIManager.mFI  = tFieldInterpolator;
    tFIManager.mIPGeometryInterpolator = &tGeomInterpolator;
    tProperty.set_field_interpolator_manager( &tFIManager );

    // evaluate the property
    Matrix< DDRMat > tPropertyValue = tProperty.val();

    //check property value
    CHECK( equal_to( tPropertyValue( 0, 0 ), 1.5 ) );

    // clean up
    // delete field interpolators
    tFieldInterpolator.clear();

}/* TEST CASE */
