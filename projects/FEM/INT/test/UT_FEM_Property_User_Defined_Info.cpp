
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Property_User_Defined_Info.hpp"              //FEM/INT/src

moris::Matrix< moris::DDRMat > tConstValFunction_UTPropContainer( moris::Cell< moris::Matrix< moris::DDRMat > >  & aParameters,
                                                                  moris::Cell< moris::fem::Field_Interpolator* > & aFieldInterpolator,
                                                                  moris::fem::Geometry_Interpolator              * aGeometryInterpolator )
{
    return aParameters( 0 );
}

namespace moris
{
    namespace fem
    {

        TEST_CASE( "Property_User_Defined_Info", "[moris],[fem],[Property_User_Defined_Info]" )
        {

            // create property user defined info
            fem::Property_User_Defined_Info tPropertyUserDefinedInfo( { fem::Property_Type::CONDUCTIVITY },
                                                                      {{ MSI::Dof_Type::TEMP }},
                                                                      {{{ 1.0 }}},
                                                                      tConstValFunction_UTPropContainer,
                                                                      { tConstValFunction_UTPropContainer } );

            //check property type
            CHECK( equal_to( static_cast< uint >( tPropertyUserDefinedInfo.get_property_type() ), 3 ) );

            //check number of dof type and dof type
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPropertyUserDefinedInfo.get_property_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // check number of parameters and parameter values
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_param_list().size(), 1.0 ) );
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_param_list()( 0 )( 0, 0 ), 1.0 ) );

            // check number of valFunc
            fem::PropertyFunc tValFunc = tPropertyUserDefinedInfo.get_property_valFunc();

            // check number of derFunc
            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_derFunc_list().size(), 1 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
