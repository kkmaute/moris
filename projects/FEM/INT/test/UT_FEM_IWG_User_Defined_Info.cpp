
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        TEST_CASE( "IWG_User_Defined_Info", "[moris],[fem],[IWG_User_Defined_Info]" )
        {

            // build IWG user defined info
            fem::IWG_User_Defined_Info tIWGUserInfo( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                     { MSI::Dof_Type::TEMP },
                                                     {{ MSI::Dof_Type::TEMP }},
                                                     { fem::Property_Type::TEMP_DIRICHLET },
                                                     { fem::Constitutive_Type::DIFF_LIN_ISO } );

            //check IWG type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_IWG_type() ), 8 ) );

            //check residual dof type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_residual_dof_type()( 0 ) ), 3 ) );

            //check number of dof type and dof type
            CHECK( equal_to( tIWGUserInfo.get_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            //check number of property type and property type
            CHECK( equal_to( tIWGUserInfo.get_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_property_type_list()( 0 ) ), 1 ) );

            //check number of constitutive type and constitutive type
            CHECK( equal_to( tIWGUserInfo.get_constitutive_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_constitutive_type_list()( 0 ) ), 1 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
