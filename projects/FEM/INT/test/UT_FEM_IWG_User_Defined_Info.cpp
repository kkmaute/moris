
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
            // list of IWG types
            moris::Cell< moris::Cell< fem::IWG_Type > >  tIWGTypeList = {{ fem::IWG_Type::SPATIALDIFF_BULK ,
                                                                           fem::IWG_Type::SPATIALDIFF_DIRICHLET,
                                                                           fem::IWG_Type::HELMHOLTZ,
                                                                           fem::IWG_Type::LSNORMAL }};
            uint tNumSets = tIWGTypeList.size();

            // list of residual dof type
            moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > tResidualDofType( tNumSets );
            tResidualDofType( 0 ).resize( tIWGTypeList( 0 ).size() );
            tResidualDofType( 0 )( 0 ) = { MSI::Dof_Type::TEMP };
            tResidualDofType( 0 )( 1 ) = { MSI::Dof_Type::TEMP };
            tResidualDofType( 0 )( 2 ) = { MSI::Dof_Type::VX };
            tResidualDofType( 0 )( 3 ) = { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ };

            // list of IWG master dof dependencies
            moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > tMasterDofTypes( tNumSets );
            tMasterDofTypes( 0 ).resize( tIWGTypeList( 0 ).size() );
            tMasterDofTypes( 0 )( 0 ) = {{ MSI::Dof_Type::TEMP }};
            tMasterDofTypes( 0 )( 1 ) = {{ MSI::Dof_Type::TEMP }};
            tMasterDofTypes( 0 )( 2 ) = {{ MSI::Dof_Type::VX }};
            tMasterDofTypes( 0 )( 3 ) = {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ },
                                         { MSI::Dof_Type::LS1}};

            // list of IWG slave dof dependencies
            moris::Cell< moris::Cell< moris::Cell< moris::Cell< MSI::Dof_Type > > > > tSlaveDofTypes;

            // list of IWG master property dependencies
            moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > tMasterPropTypes( tNumSets );
            tMasterPropTypes( 0 ).resize( tIWGTypeList( 0 ).size() );
            tMasterPropTypes( 0 )( 0 ) = { fem::Property_Type::CONDUCTIVITY };
            tMasterPropTypes( 0 )( 1 ) = { fem::Property_Type::CONDUCTIVITY,
                                           fem::Property_Type::TEMP_DIRICHLET };

            // list of IWG master property dependencies
            moris::Cell< moris::Cell< moris::Cell< fem::Property_Type > > > tSlavePropTypes;

            // build an IWG user defined info
            IWG_User_Defined_Info tIWGUserDefinedInfo( tIWGTypeList,
                                                       tResidualDofType,
                                                       tMasterDofTypes, tMasterPropTypes,
                                                       tSlaveDofTypes,  tSlavePropTypes );

            // a factory to create the IWGs
            fem::IWG_Factory tIWGFactory;

            // create a cell of IWGs
            moris::Cell< moris::Cell< fem::IWG* > > tIWGs( tNumSets );

            // loop over the sets
            for( uint iSet = 0; iSet < tNumSets; iSet++ )
            {
                // set size for the cell of IWGs for the set
                tIWGs( iSet ).resize( tIWGTypeList( iSet ).size(), nullptr );

                // loop over the IWG types for the set
                for( uint iIWG = 0; iIWG < tIWGTypeList( iSet ).size(); iIWG++ )
                {
                    // create an IWG with the factory for the IWG type
                    tIWGs( iSet )( iIWG ) = tIWGFactory.create_IWGs( tIWGTypeList( iSet )( iIWG ) );

                    // set residual dof type
                    tIWGs( iSet )( iIWG )->set_residual_dof_type( tIWGUserDefinedInfo.get_residual_dof_type()( iSet )( iIWG ) );

                    // set active dof type
                    tIWGs( iSet )( iIWG )->set_dof_type_list( tIWGUserDefinedInfo.get_dof_type_list()( iSet )( iIWG ) );

                    // set active property type
                    tIWGs( iSet )( iIWG )->set_property_type_list( tIWGUserDefinedInfo.get_property_type_list()( iSet )( iIWG ) );
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
