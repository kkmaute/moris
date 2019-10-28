
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        // This test case tests all the member functions in IWG_User_Defined_Info.
        TEST_CASE( "IWG_User_Defined_Info_Constructor", "[moris],[fem],[IWG_User_Defined_Info_Constructor]" )
        {

            // build IWG user defined info
            fem::IWG_User_Defined_Info tIWGUserInfo( fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                     { MSI::Dof_Type::TEMP },
                                                     {{ MSI::Dof_Type::TEMP }},
                                                     {{ MSI::Dv_Type::LS1 }},
                                                     { fem::Property_Type::TEMP_DIRICHLET },
                                                     { fem::Constitutive_Type::DIFF_LIN_ISO },
                                                     { fem::Penalty_Type::DIRICHLET_NITSCHE },
                                                     {{ MSI::Dof_Type::TEMP }},
                                                     {{ MSI::Dv_Type::LS1 }},
                                                     { fem::Property_Type::TEMP_DIRICHLET },
                                                     { fem::Constitutive_Type::DIFF_LIN_ISO } );

            // check IWG type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_IWG_type() ), 8 ) );

            // check residual dof type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_residual_dof_type()( 0 ) ), 3 ) );

            // check master number of dof type and dof type
            CHECK( equal_to( tIWGUserInfo.get_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // check slave number of dof type and dof type
            CHECK( equal_to( tIWGUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE )( 0 )( 0 ) ), 3 ) );

            // check master number of property type and property type
            CHECK( equal_to( tIWGUserInfo.get_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_property_type_list()( 0 ) ), 1 ) );

            // check penalty type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_penalty_type_list().size() ), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_penalty_type_list()( 0 ) ), 1 ) );

            // check slave number of property type and property type
            CHECK( equal_to( tIWGUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );

            //check master number of constitutive type and constitutive type
            CHECK( equal_to( tIWGUserInfo.get_constitutive_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_constitutive_type_list()( 0 ) ), 1 ) );

            //check slave number of constitutive type and constitutive type
            CHECK( equal_to( tIWGUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );
        }/* TEST_CASE */

        TEST_CASE( "IWG_User_Defined_Info_Set", "[moris],[fem],[IWG_User_Defined_Info_Set]" )
        {

            // build IWG user defined info
            fem::IWG_Type                  tIWGType         = fem::IWG_Type::SPATIALDIFF_DIRICHLET;
            Cell< MSI::Dof_Type >          tResidualDofType = { MSI::Dof_Type::TEMP };
            Cell< Cell< MSI::Dof_Type > >  tDofTypes        = {{ MSI::Dof_Type::TEMP }};
            Cell< Cell< MSI::Dv_Type > >   tDvTypes         = {{ MSI::Dv_Type::LS1 }};
            Cell< fem::Property_Type>      tPropTypes       = { fem::Property_Type::TEMP_DIRICHLET };
            Cell< fem::Constitutive_Type > tConstTypes      = { fem::Constitutive_Type::DIFF_LIN_ISO };
            Cell< fem::Penalty_Type >      tPenaltyTypes    = { fem::Penalty_Type::DIRICHLET_NITSCHE };

            fem::IWG_User_Defined_Info tIWGUserInfo;
            tIWGUserInfo.set_IWG_type( tIWGType );
            tIWGUserInfo.set_residual_dof_type( tResidualDofType );
            tIWGUserInfo.set_dof_type_list( tDofTypes );
            tIWGUserInfo.set_dv_type_list( tDvTypes );
            tIWGUserInfo.set_property_type_list( tPropTypes );
            tIWGUserInfo.set_constitutive_type_list( tConstTypes );
            tIWGUserInfo.set_penalty_type_list( tPenaltyTypes );
            tIWGUserInfo.set_dof_type_list( tDofTypes, mtk::Master_Slave::SLAVE );
            tIWGUserInfo.set_dv_type_list( tDvTypes, mtk::Master_Slave::SLAVE );
            tIWGUserInfo.set_property_type_list( tPropTypes, mtk::Master_Slave::SLAVE );
            tIWGUserInfo.set_constitutive_type_list( tConstTypes, mtk::Master_Slave::SLAVE );

            // check IWG type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_IWG_type() ), 8 ) );

            // check residual dof type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_residual_dof_type()( 0 ) ), 3 ) );

            // check master number of dof type and dof type
            CHECK( equal_to( tIWGUserInfo.get_dof_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_dof_type_list()( 0 )( 0 ) ), 3 ) );

            // check slave number of dof type and dof type
            CHECK( equal_to( tIWGUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_dof_type_list( mtk::Master_Slave::SLAVE )( 0 )( 0 ) ), 3 ) );

            // check master number of property type and property type
            CHECK( equal_to( tIWGUserInfo.get_property_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_property_type_list()( 0 ) ), 1 ) );

            // check penalty type
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_penalty_type_list().size() ), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_penalty_type_list()( 0 ) ), 1 ) );

            // check slave number of property type and property type
            CHECK( equal_to( tIWGUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_property_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );

            //check master number of constitutive type and constitutive type
            CHECK( equal_to( tIWGUserInfo.get_constitutive_type_list().size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_constitutive_type_list()( 0 ) ), 1 ) );

            //check slave number of constitutive type and constitutive type
            CHECK( equal_to( tIWGUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE ).size(), 1 ) );
            CHECK( equal_to( static_cast< uint >( tIWGUserInfo.get_constitutive_type_list( mtk::Master_Slave::SLAVE )( 0 ) ), 1 ) );
        }/* TEST_CASE */


    }/* namespace fem */
}/* namespace moris */
