
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_Penalty_User_Defined_Info.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        // This test case tests all the member functions in Penalty_User_Defined_Info.
        TEST_CASE( "Penalty_User_Defined_Info_Constructor", "[moris],[fem],[Penalty_User_Defined_Info_Constructor]" )
        {

            // build penalty user defined info
            fem::Penalty_User_Defined_Info tPenaltyUserInfo( fem::Penalty_Type::DIRICHLET_NITSCHE,
                                                             {{ MSI::Dof_Type::TEMP }},
                                                             {{ MSI::Dv_Type::LS1 }},
                                                             { fem::Property_Type::TEMP_DIRICHLET },
                                                             { fem::Constitutive_Type::DIFF_LIN_ISO },
                                                             {{ MSI::Dof_Type::TEMP }},
                                                             {{ MSI::Dv_Type::LS1 }},
                                                             { fem::Property_Type::TEMP_DIRICHLET },
                                                             { fem::Constitutive_Type::DIFF_LIN_ISO } );

            // check penalty type
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_penalty_type() ), 1 ) );

            // check master number of dof type and dof type
            CHECK( equal_to( tPenaltyUserInfo.get_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // check slave number of dof type and dof type
            CHECK( equal_to( tPenaltyUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE )( 0 )( 0 ) ), 3 ) );

            // check master number of property type and property type
            CHECK( equal_to( tPenaltyUserInfo.get_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_property_type_list()( 0 ) ), 1 ) );

            // check slave number of property type and property type
            CHECK( equal_to( tPenaltyUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );

            //check master number of constitutive type and constitutive type
            CHECK( equal_to( tPenaltyUserInfo.get_constitutive_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_constitutive_type_list()( 0 ) ), 1 ) );

            //check slave number of constitutive type and constitutive type
            CHECK( equal_to( tPenaltyUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );
        }/* TEST_CASE */

        TEST_CASE( "Penalty_User_Defined_Info_Set", "[moris],[fem],[Penalty_User_Defined_Info_Set]" )
        {

            // build penalty user defined info
            fem::Penalty_Type              tPenaltyType     = fem::Penalty_Type::DIRICHLET_NITSCHE;
            Cell< MSI::Dof_Type >          tResidualDofType = { MSI::Dof_Type::TEMP };
            Cell< Cell< MSI::Dof_Type > >  tDofTypes        = {{ MSI::Dof_Type::TEMP }};
            Cell< Cell< MSI::Dv_Type > >   tDvTypes         = {{ MSI::Dv_Type::LS1 }};
            Cell< fem::Property_Type>      tPropTypes       = { fem::Property_Type::TEMP_DIRICHLET };
            Cell< fem::Constitutive_Type > tConstTypes      = { fem::Constitutive_Type::DIFF_LIN_ISO };

            fem::Penalty_User_Defined_Info tPenaltyUserInfo;
            tPenaltyUserInfo.set_penalty_type( tPenaltyType );
            tPenaltyUserInfo.set_dof_type_list( tDofTypes );
            tPenaltyUserInfo.set_dv_type_list( tDvTypes );
            tPenaltyUserInfo.set_property_type_list( tPropTypes );
            tPenaltyUserInfo.set_constitutive_type_list( tConstTypes );
            tPenaltyUserInfo.set_dof_type_list( tDofTypes, mtk::Master_Slave::SLAVE );
            tPenaltyUserInfo.set_dv_type_list( tDvTypes, mtk::Master_Slave::SLAVE );
            tPenaltyUserInfo.set_property_type_list( tPropTypes, mtk::Master_Slave::SLAVE );
            tPenaltyUserInfo.set_constitutive_type_list( tConstTypes, mtk::Master_Slave::SLAVE );

            // check penalty type
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_penalty_type() ), 1 ) );

            // check master number of dof type and dof type
            CHECK( equal_to( tPenaltyUserInfo.get_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // check slave number of dof type and dof type
            CHECK( equal_to( tPenaltyUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE )( 0 )( 0 ) ), 3 ) );

            // check master number of property type and property type
            CHECK( equal_to( tPenaltyUserInfo.get_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_property_type_list()( 0 ) ), 1 ) );

            // check slave number of property type and property type
            CHECK( equal_to( tPenaltyUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );

            //check master number of constitutive type and constitutive type
            CHECK( equal_to( tPenaltyUserInfo.get_constitutive_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_constitutive_type_list()( 0 ) ), 1 ) );

            //check slave number of constitutive type and constitutive type
            CHECK( equal_to( tPenaltyUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tPenaltyUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );
        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
