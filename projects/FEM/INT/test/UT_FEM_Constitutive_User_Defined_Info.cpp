
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Constitutive_User_Defined_Info.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        TEST_CASE( "Constitutive_User_Defined_Info", "[moris],[fem],[Constitutive_User_Defined_Info]" )
        {

            // build constitutive user defined info
            fem::Constitutive_User_Defined_Info tConstitutiveUserInfo( fem::Constitutive_Type::DIFF_LIN_ISO,
                                                                       {{ MSI::Dof_Type::TEMP }},
                                                                       { fem::Property_Type::CONDUCTIVITY } );

            //check constitutive type
            CHECK( equal_to( static_cast< uint >( tConstitutiveUserInfo.get_constitutive_type() ), 1 ) );

            //check number of dof type and dof type
            CHECK( equal_to( tConstitutiveUserInfo.get_constitutive_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tConstitutiveUserInfo.get_constitutive_dof_type_list()( 0 )( 0 ) ), 3 ) );

            //check number of property type and property type
            CHECK( equal_to( tConstitutiveUserInfo.get_constitutive_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tConstitutiveUserInfo.get_constitutive_property_type_list()( 0 ) ), 3 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
