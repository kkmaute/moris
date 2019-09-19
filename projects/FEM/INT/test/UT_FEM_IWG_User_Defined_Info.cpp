
#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "cl_FEM_IWG_User_Defined_Info.hpp"              //FEM/INT/src
#include "cl_FEM_IWG.hpp"              //FEM/INT/src
#include "cl_FEM_IWG_Factory.hpp"              //FEM/INT/src

namespace moris
{
    namespace fem
    {

        TEST_CASE( "IWG_User_Defined_Info", "[moris],[fem],[IWG_User_Defined_Info]" )
        {

            // build IWG user defined info
            moris::Cell< moris::Cell< fem::IWG_User_Defined_Info > > tIWGUserInfo( 1 );
            tIWGUserInfo( 0 ).resize( 4 );

            tIWGUserInfo( 0 )( 0 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_BULK, 3, { MSI::Dof_Type::TEMP },
                                                                 {{ MSI::Dof_Type::TEMP }},
                                                                 { fem::Property_Type::CONDUCTIVITY },
                                                                 moris::Cell< fem::Constitutive_Type >( 0 ) );
            tIWGUserInfo( 0 )( 1 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::SPATIALDIFF_DIRICHLET, 3, { MSI::Dof_Type::TEMP },
                                                                 {{ MSI::Dof_Type::TEMP }},
                                                                 { fem::Property_Type::CONDUCTIVITY, fem::Property_Type::TEMP_DIRICHLET },
                                                                 moris::Cell< fem::Constitutive_Type >( 0 ) );
            tIWGUserInfo( 0 )( 2 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::HELMHOLTZ, 3, { MSI::Dof_Type::VX },
                                                                 {{ MSI::Dof_Type::VX }},
                                                                 moris::Cell< fem::Property_Type >( 0 ),
                                                                 moris::Cell< fem::Constitutive_Type >( 0 ) );
            tIWGUserInfo( 0 )( 3 ) = fem::IWG_User_Defined_Info( fem::IWG_Type::LSNORMAL, 3, { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ },
                                                                 {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ }, { MSI::Dof_Type::LS1}},
                                                                 moris::Cell< fem::Property_Type >( 0 ),
                                                                 moris::Cell< fem::Constitutive_Type >( 0 ) );

            // number of sets
            uint tNumSets = tIWGUserInfo.size();

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs
            moris::Cell< moris::Cell< fem::IWG* > > tIWGs( tNumSets );

            // loop over the sets
            for( uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                // number of IWG
                uint tNumIWG = tIWGUserInfo( iSet ).size();

                // set size for the cell of IWGs for the set
                tIWGs( iSet ).resize( tNumIWG, nullptr );

                // loop over the IWG types for the set
                for( uint iIWG = 0; iIWG < tNumIWG; iIWG++ )
                {
                    // create an IWG with the factory for the IWG type
                    tIWGs( iSet )( iIWG ) = tIWGFactory.create_IWGs( tIWGUserInfo( iSet )( iIWG ).get_IWG_type() );

                    // set residual dof type
                    tIWGs( iSet )( iIWG )->set_space_dim( tIWGUserInfo( iSet )( iIWG ).get_IWG_space_dim() );

                    // set residual dof type
                    tIWGs( iSet )( iIWG )->set_residual_dof_type( tIWGUserInfo( iSet )( iIWG ).get_residual_dof_type() );

                    // set active dof type
                    tIWGs( iSet )( iIWG )->set_dof_type_list( tIWGUserInfo( iSet )( iIWG ).get_dof_type_list() );

                    // set active property type
                    tIWGs( iSet )( iIWG )->set_property_type_list( tIWGUserInfo( iSet )( iIWG ).get_property_type_list() );
                }
            }

            // clean up
            for( uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                for( IWG* tIWG : tIWGs( iSet ) )
                {
                    delete tIWG;
                }
            }

            //check property map content
//            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 3, 0 ), 0 ) );
//            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 1, 0 ), 1 ) );
//            CHECK( equal_to( tPropertyUserDefinedInfo.get_property_map()( 2, 0 ), 2 ) );

        }/* TEST_CASE */

    }/* namespace fem */
}/* namespace moris */
