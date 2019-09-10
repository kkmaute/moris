
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >    & aCoeff,
                                            moris::Cell< Field_Interpolator* > & aFieldInterpolator )
        {
            return aCoeff( 0 );
        }

        TEST_CASE( "Property_User_Defined_Info", "[moris],[fem],[Property_User_Defined_Info]" )
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
            PropertyFunc tValFunction0 = tConstValFunction;

            // create the list with function pointers for the value
            Cell< PropertyFunc > tValFuncList( 3, tValFunction0 );

            // create the list with cell of function pointers for the derivatives
            Cell< Cell< PropertyFunc > > tDerFuncList( 3 );

            // collect properties info
            fem::Property_User_Defined_Info tPropertyUserDefinedInfo( tPropertyTypeList,
                                                                      tPropertyDofList,
                                                                      tCoeffList,
                                                                      tValFuncList,
                                                                      tDerFuncList );

            //check property map content
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 3, 0 ), 0 ) );
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 1, 0 ), 1 ) );
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 2, 0 ), 2 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
