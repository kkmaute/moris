
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

moris::Matrix< moris::DDRMat > tValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                             moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    moris::Matrix< moris::DDRMat > tPropertyVal( 1, 1, 1.0 );
    return tPropertyVal;
}

moris::Matrix< moris::DDRMat > tDerFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                             moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                             moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    moris::Matrix< moris::DDRMat > tPropertyDer( 1, 1, 2.0);
    return tPropertyDer;
}

moris::Matrix< moris::DDRMat > tValFunction2( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                              moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) + aParameters( 1 ) * aFieldInterpolator( 0 )->val() + aParameters( 2 ) * aFieldInterpolator( 1 )->val();
}

moris::Matrix< moris::DDRMat > tDerFunction2( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                              moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 1 ) * aFieldInterpolator( 0 )->N();
}

moris::Matrix< moris::DDRMat > tDerFunction3( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                              moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                              moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 2 ) * aFieldInterpolator( 1 )->N();
}

moris::Matrix< moris::DDRMat > tConstValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                  moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                  moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

moris::Matrix< moris::DDRMat > tGeoValFunction( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 ) * aGeometryInterpolator->valx();
}

namespace moris
{
    namespace fem
    {

        TEST_CASE( "Property", "[moris],[fem],[Property]" )
        {
            //create a space geometry interpolation rule
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            //create a space and a time geometry interpolator
            Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

            // set coeffs
            Cell< Matrix< DDRMat > > tCoeff;

            // create a property object
            fem::Property tProperty( fem::Property_Type::TEMP_DIRICHLET,
                                     {{ MSI::Dof_Type::TEMP }},
                                     tCoeff,
                                     tValFunction,
                                     { tDerFunction },
                                     tGeomInterpolator );

            // check property type
            CHECK( equal_to( static_cast< uint >( tProperty.get_property_type() ), 1 ) );

            //check dof dependencies
            CHECK( equal_to( static_cast< uint >( tProperty.get_dof_type_list()( 0 )( 0 ) ), 3 ) );


            // set field interpolators
            Cell< Field_Interpolator* > tFieldInterpolator( 1 );
            tFieldInterpolator( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP });
            tProperty.set_field_interpolators( tFieldInterpolator );

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
            // delete field interpolators
            for( Field_Interpolator* tFI : tFieldInterpolator )
            {
                delete tFI;
            }
            tFieldInterpolator.clear();

            // delete geometry interpolator
            delete tGeomInterpolator;
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
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            //create a space and a time geometry interpolator
            Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

            //set the coefficients xHat, tHat
            tGeomInterpolator->set_coeff( tXHat, tTHat );

            // set the property coefficients
            Cell< Matrix< DDRMat > > tCoeff( 3 );
            tCoeff( 0 ) = {{ 1.0 }};
            tCoeff( 1 ) = {{ 2.0 }};
            tCoeff( 2 ) = {{ 3.0 }};

            // create a property object
            fem::Property tProperty( fem::Property_Type::CONDUCTIVITY,
                                     {{ MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX }},
                                     tCoeff,
                                     tValFunction2,
                                     { tDerFunction2, tDerFunction3 },
                                     tGeomInterpolator );

            // create an interpolation rule
            Interpolation_Rule tInterpolationRule ( mtk::Geometry_Type::QUAD,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR,
                                                    Interpolation_Type::LAGRANGE,
                                                    mtk::Interpolation_Order::LINEAR );
            // create a TEMP field interpolator
            uint tNumberOfFields = 1;
            Cell< Field_Interpolator* > tFieldInterpolator( 2, nullptr );
            tFieldInterpolator( 0 ) = new Field_Interpolator ( tNumberOfFields,
                                                               tInterpolationRule,
                                                               tGeomInterpolator,
                                                               { MSI::Dof_Type::TEMP } );
            tFieldInterpolator( 1 ) = new Field_Interpolator ( tNumberOfFields,
                                                               tInterpolationRule,
                                                               tGeomInterpolator,
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

            // set coeffs and field interpolators
            tProperty.set_field_interpolators( tFieldInterpolator );

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
            // delete geometry interpolator
            delete tGeomInterpolator;

            // delete field interpolators
            for( Field_Interpolator* tFI : tFieldInterpolator )
            {
                delete tFI;
            }
            tFieldInterpolator.clear();

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
            Interpolation_Rule tGeomInterpRule( mtk::Geometry_Type::QUAD,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR,
                                                Interpolation_Type::LAGRANGE,
                                                mtk::Interpolation_Order::LINEAR );

            //create a space and a time geometry interpolator
            Geometry_Interpolator* tGeomInterpolator = new Geometry_Interpolator( tGeomInterpRule );

            //set the coefficients xHat, tHat
            tGeomInterpolator->set_coeff( tXHat, tTHat );

            //set the evaluation point tParamPoint
            Matrix< DDRMat > tParamPoint = {{ 0.0 }, { 0.0 }, { 0.0 }};
            tGeomInterpolator->set_space_time( tParamPoint );

            // create property coeffs
            Cell< Matrix< DDRMat > > tCoeff( 1 );
            tCoeff( 0 ) = {{ 1.0 }};

            // create a property object
            fem::Property tProperty( fem::Property_Type::TEMP_DIRICHLET,
                                     {{ MSI::Dof_Type::TEMP }},
                                     tCoeff,
                                     tGeoValFunction,
                                     Cell< PropertyFunc >( 0 ),
                                     tGeomInterpolator );

            // check property type
            CHECK( equal_to( static_cast< uint >( tProperty.get_property_type() ), 1 ) );

            //check dof dependencies
            CHECK( equal_to( static_cast< uint >( tProperty.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // set field interpolators
            Cell< Field_Interpolator* > tFieldInterpolator( 1 );
            tFieldInterpolator( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP });
            tProperty.set_field_interpolators( tFieldInterpolator );

            // evaluate the property
            Matrix< DDRMat > tPropertyValue = tProperty.val();

            //check property value
            CHECK( equal_to( tPropertyValue( 0, 0 ), 1.5 ) );

            // clean up
            delete tFieldInterpolator( 0 );
            delete tGeomInterpolator;
        }/* TEST CASE */

    }/* namespace fem */
}/* namespace moris */
