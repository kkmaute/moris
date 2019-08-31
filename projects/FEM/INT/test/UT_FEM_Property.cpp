
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property.hpp"              //FEM/INT/src
#include "cl_FEM_Geometry_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_Field_Interpolator.hpp"         //FEM/INT/src
#include "cl_FEM_IWG.hpp"         //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"         //FEM/INT/src

namespace moris
{
    namespace fem
    {

        Matrix< DDRMat > tValFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                       moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyVal( 1, 1, 1.0);
            return tPropertyVal;
        }

        Matrix< DDRMat > tDerFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                       moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer( 1, 1, 2.0);
            return tPropertyDer;
        }

        Matrix< DDRMat > tValFunction2( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyVal( 1, 1, 0.0);
            tPropertyVal = aCoeff( 0 ) + aCoeff( 1 ) * aFieldInterpolator( 0 )->val() + aCoeff( 2 ) * aFieldInterpolator( 1 )->val();
            return tPropertyVal;
        }

        Matrix< DDRMat > tDerFunction2( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer;
            tPropertyDer = aCoeff( 1 ) * aFieldInterpolator( 0 )->N();
            return tPropertyDer;
        }

        Matrix< DDRMat > tDerFunction3( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                        moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            Matrix< DDRMat > tPropertyDer;
            tPropertyDer = aCoeff( 2 ) * aFieldInterpolator( 1 )->N();
            return tPropertyDer;
        }

        Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                            moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            return aCoeff( 0 );
        }

        TEST_CASE( "Property", "[moris],[fem],[Property]" )
        {
            // create a list of dof dependencies for the property
            Cell< Cell< MSI::Dof_Type > >  tActiveDofTypes = {{ MSI::Dof_Type::TEMP }};

            // create the function pointers for the value
            fem::PropertyFunc tValFunction0 = tValFunction;

            // create the cell of function pointers for the derivatives
            fem::PropertyFunc tDerFunction0 = tDerFunction;
            Cell< fem::PropertyFunc > tDerFunctions( 1 );
            tDerFunctions( 0 ) = tDerFunction0;

            // create a property object
            fem::Property tProperty( fem::Property_Type::TEMP_DIRICHLET,
                                     tActiveDofTypes,
                                     tValFunction0,
                                     tDerFunctions );

            // check property type
            CHECK( equal_to( static_cast< uint >( tProperty.get_property_type() ), 1 ) );

            //check dof dependencies
            CHECK( equal_to( static_cast< uint >( tProperty.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // set coeffs and field interpolators
            Cell< Matrix< DDRMat > > tCoeff;
            tProperty.set_coefficients( tCoeff );
            Cell< Field_Interpolator* > tFieldInterpolator( 1 );
            tFieldInterpolator( 0 ) = new Field_Interpolator( 1, { MSI::Dof_Type::TEMP });
            tProperty.set_field_interpolators( tFieldInterpolator );

            // evaluate the property
            Matrix< DDRMat > tPropertyValue = tProperty.val();

            //check property value
            CHECK( equal_to( tPropertyValue( 0, 0 ), 1.0 ) );

            // evaluate the property derivative wrt to TEMP (in dependencies)
            Matrix< DDRMat > tPropertyDerivative = tProperty.dPropdDOF( { MSI::Dof_Type::TEMP } );

            //check property value
            CHECK( equal_to( tPropertyDerivative( 0, 0 ), 2.0 ) );

            // evaluate the property derivative wrt to UX (not in dependencies)
            tPropertyDerivative = tProperty.dPropdDOF( { MSI::Dof_Type::UX } );

            //check property value
            CHECK( equal_to( tPropertyDerivative( 0, 0 ), 0.0 ) );

            // clean up
            delete tFieldInterpolator( 0 );
        }

        TEST_CASE( "Property_with_dependency", "[moris],[fem],[Property_with_dependency]" )
        {
            // create a list of dof dependencies for the property
            Cell< Cell< MSI::Dof_Type > >  tActiveDofTypes = {{ MSI::Dof_Type::TEMP },
                                                              { MSI::Dof_Type::UX }};

            // create the function pointers for the value
            fem::PropertyFunc tValFunction0 = tValFunction2;

            // create the cell of function pointers for the derivatives
            fem::PropertyFunc tDerFunction0 = tDerFunction2;
            fem::PropertyFunc tDerFunction1 = tDerFunction3;
            Cell< fem::PropertyFunc > tDerFunctions( 2, nullptr );
            tDerFunctions( 0 ) = tDerFunction0;
            tDerFunctions( 1 ) = tDerFunction1;

            // create a property object
            fem::Property tProperty( fem::Property_Type::CONDUCTIVITY,
                                     tActiveDofTypes,
                                     tValFunction0,
                                     tDerFunctions );

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

            // set the property coefficients
            Cell< Matrix< DDRMat > > tCoeff( 3 );
            tCoeff( 0 ) = {{ 1.0 }};
            tCoeff( 1 ) = {{ 2.0 }};
            tCoeff( 2 ) = {{ 3.0 }};

            // set coeffs and field interpolators
            tProperty.set_field_interpolators( tFieldInterpolator );
            tProperty.set_coefficients( tCoeff );

            //evaluate the property
            Matrix< DDRMat > tPropertyValue = tProperty.val();
            //print( tPropertyValue, "tPropertyValue" );

            // evaluate the property derivative wrt to TEMP (in dependencies)
            Matrix< DDRMat > tPropertyDerivative = tProperty.dPropdDOF( MSI::Dof_Type::TEMP );
            //print( tPropertyDerivative, "tPropertyDerivative" );

            // evaluate the property derivative wrt to UX (in dependencies)
            tPropertyDerivative = tProperty.dPropdDOF( MSI::Dof_Type::UX );
            //print( tPropertyDerivative, "tPropertyDerivative" );

            // evaluate the property derivative wrt to LS1 (not in dependencies)
            tPropertyDerivative = tProperty.dPropdDOF( MSI::Dof_Type::LS1 );
            //print( tPropertyDerivative, "tPropertyDerivative" );

            // clean up
            delete tGeomInterpolator;
            delete tFieldInterpolator( 0 );
            delete tFieldInterpolator( 1 );

        }/* TEST_CASE */

        TEST_CASE( "Property_creation_model", "[moris],[fem],[Property_creation_model]" )
        {
            // list of property type
            Cell< fem::Property_Type > tPropertyTypeList = {{fem::Property_Type::CONDUCTIVITY},
                                                            {fem::Property_Type::TEMP_DIRICHLET}};

            // list of property dependencies
            Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 2 );
            tPropertyDofList( 0 ) = {{ MSI::Dof_Type::TEMP},
                                     { MSI::Dof_Type::UX  }};
            tPropertyDofList( 1 ) = {{ MSI::Dof_Type::TEMP},
                                     { MSI::Dof_Type::UX  }};

            // list of the property coefficients
            Cell< Cell< Matrix< DDRMat > > > tCoeffList( 2 );
            tCoeffList( 0 ).resize( 3 );
            tCoeffList( 0 )( 0 )= {{ 1.0 }};
            tCoeffList( 0 )( 1 )= {{ 2.0 }};
            tCoeffList( 0 )( 2 )= {{ 3.0 }};
            tCoeffList( 1 ).resize( 3 );
            tCoeffList( 1 )( 0 )= {{ 1.0 }};
            tCoeffList( 1 )( 1 )= {{ 2.0 }};
            tCoeffList( 1 )( 2 )= {{ 3.0 }};

            // cast free function into std::function
            std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                              moris::Cell< Field_Interpolator* > & aFieldInterpolator) > tValFunction0 = tValFunction2;

            // create the list with function pointers for the value
            Cell< std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                                    moris::Cell< Field_Interpolator* > & aFieldInterpolator) > > tValFuncList( 2, nullptr );
            tValFuncList( 0 ) = tValFunction0;
            tValFuncList( 1 ) = tValFunction0;

            // cast free function into std::function
            std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                              moris::Cell< Field_Interpolator* > & aFieldInterpolator) > tDerFunction0 = tDerFunction2;
            std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                              moris::Cell< Field_Interpolator* > & aFieldInterpolator) > tDerFunction1 = tDerFunction3;

            // create the list with cell of function pointers for the derivatives
            Cell< Cell< std::function< Matrix< DDRMat > ( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                                          moris::Cell< Field_Interpolator* > & aFieldInterpolator) > > > tDerFuncList( 2 );
            tDerFuncList( 0 ).resize( 2, nullptr );
            tDerFuncList( 1 ).resize( 2, nullptr );
            tDerFuncList( 0 )( 0 ) = tDerFunction0;
            tDerFuncList( 0 )( 1 ) = tDerFunction1;
            tDerFuncList( 1 )( 0 ) = tDerFunction0;
            tDerFuncList( 1 )( 1 ) = tDerFunction1;

            // model property map
            sint tMaxEnum = 0;
            for( uint iProp = 0; iProp < tPropertyTypeList.size(); iProp++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( tPropertyTypeList( iProp ) ) );
            }
            tMaxEnum++;

            // creation of properties
            Matrix< DDSMat > tPropertyTypeMap( tMaxEnum, 1, -1 );
            Cell< Property* > tModelProperties( tPropertyTypeList.size() );
            for( uint iProp = 0; iProp < tPropertyTypeList.size(); iProp++ )
            {
                tModelProperties( iProp ) = new Property( tPropertyTypeList( iProp ),
                                                          tPropertyDofList( iProp ),
                                                          tValFuncList( iProp ),
                                                          tDerFuncList( iProp ) );

                tPropertyTypeMap( static_cast< int >( tPropertyTypeList( iProp ) ), 0 ) = iProp;
            }
            //print(tPropertyTypeMap,"tPropertyTypeMap");

            // clean up
            for( Property* tProperty : tModelProperties )
            {
                delete tProperty;
            }

        }/* TEST_CASE */

        TEST_CASE( "Property_model_diffusion", "[moris],[fem],[Property_model_diffusion]" )
        {
            // list of property type
            Cell< fem::Property_Type > tPropertyTypeList = {{ fem::Property_Type::CONDUCTIVITY   },
                                                            { fem::Property_Type::TEMP_DIRICHLET },
                                                            { fem::Property_Type::TEMP_NEUMANN   }};

            // list of property dependencies
            Cell< Cell< Cell< MSI::Dof_Type > > > tPropertyDofList( 3 );

            // list of the property coefficients
            Cell< Cell< Matrix< DDRMat > > > tCoeffList( 3 );
            tCoeffList( 0 ).resize( 1 );
            tCoeffList( 0 )( 0 )= {{ 1.0 }};
            tCoeffList( 1 ).resize( 1 );
            tCoeffList( 1 )( 0 )= {{ 5.0 }};
            tCoeffList( 2 ).resize( 1 );
            tCoeffList( 2 )( 0 )= {{ 20.0 }};

            // cast free function into std::function
            fem::PropertyFunc tValFunction0 = tConstValFunction;

            // create the list with function pointers for the value
            Cell< fem::PropertyFunc > tValFuncList( 3, tValFunction0 );

            // create the list with cell of function pointers for the derivatives
            Cell< Cell< fem::PropertyFunc > > tDerFuncList( 3 );

            // model property map
            sint tMaxEnum = 0;
            for( uint iProp = 0; iProp < tPropertyTypeList.size(); iProp++ )
            {
                tMaxEnum = std::max( tMaxEnum, static_cast< int >( tPropertyTypeList( iProp ) ) );
            }
            tMaxEnum++;

            // creation of properties
            Matrix< DDSMat > tPropertyTypeMap( tMaxEnum, 1, -1 );
            Cell< Property* > tModelProperties( tPropertyTypeList.size() );
            for( uint iProp = 0; iProp < tPropertyTypeList.size(); iProp++ )
            {
                tModelProperties( iProp ) = new Property( tPropertyTypeList( iProp ),
                                                          tPropertyDofList( iProp ),
                                                          tValFuncList( iProp ),
                                                          tDerFuncList( iProp ) );

                tPropertyTypeMap( static_cast< int >( tPropertyTypeList( iProp ) ), 0 ) = iProp;
            }
            //print(tPropertyTypeMap,"tPropertyTypeMap");

            // create a list of IWG type
            Cell< Cell< fem::IWG_Type > >tIWGTypeList( 3 );
            tIWGTypeList( 0 ).resize( 1, fem::IWG_Type::SPATIALDIFF_BULK );
            tIWGTypeList( 1 ).resize( 1, fem::IWG_Type::SPATIALDIFF_DIRICHLET );
            tIWGTypeList( 2 ).resize( 1, fem::IWG_Type::SPATIALDIFF_NEUMANN );

            // number of IWGs to be created
            uint tNumOfIWGs = tIWGTypeList.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs for the problem considered
            Cell< Cell< IWG* > > tIWGs( tNumOfIWGs );

            // loop over the IWG types
            for( uint i = 0; i < tNumOfIWGs; i++)
            {
                tIWGs( i ).resize( tIWGTypeList( i ).size(), nullptr );

                for( uint Ki = 0; Ki < tIWGTypeList( i ).size(); Ki++)
                {
                    // create an IWG with the factory for the ith IWG type
                    tIWGs( i )( Ki ) = tIWGFactory.create_IWGs( tIWGTypeList( i )( Ki ) );

                    // get the IWG properties
                    Cell< fem::Property_Type > tIWGPropertyType;
                    tIWGPropertyType = tIWGs( i )( Ki )->get_property_type_list();

                    // loop over property type
                    Cell< Property* > tIWGProperties( tIWGPropertyType.size() );
                    for( uint iProp = 0; iProp < tIWGPropertyType.size(); iProp++ )
                    {
                        // collect the properties that are active for the IWG
                        uint propIndex = tPropertyTypeMap( static_cast< int >( tIWGPropertyType( iProp ) ) );
                        tIWGProperties( iProp ) = tModelProperties( propIndex );
                    }

                    // set the IWG properties
                    tIWGs( i )( Ki )->set_properties( tIWGProperties );
                }
            }

            // clean up
            for( Property* tProperty : tModelProperties )
            {
                delete tProperty;
            }

            for( uint iIWGSet = 0; iIWGSet < tNumOfIWGs; iIWGSet++ )
            {
                for( uint iIWG = 0; iIWG < tIWGs( iIWGSet ).size(); iIWG++ )
                {
                    delete tIWGs( iIWGSet )( iIWG );
                }
            }

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
