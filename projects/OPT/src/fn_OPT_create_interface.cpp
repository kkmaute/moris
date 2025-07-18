/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_OPT_create_interface.cpp
 *
 */

#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_OPT_Interface_Manager.hpp"
#include "cl_Submodule_Parameter_Lists.hpp"

namespace moris::opt
{
    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Criteria_Interface >
    create_interface(
            Submodule_Parameter_Lists                        aParameterLists,
            Vector< std::shared_ptr< Criteria_Interface > > aInterfaces )
    {
        // Get number of interfaces
        uint tNumInterfaces = aParameterLists.size() + aInterfaces.size() - 1;

        // Single interface without manager
        if ( tNumInterfaces == 0 )
        {
            if ( aParameterLists.size() > 0 )
            {
                return create_interface( aParameterLists( 0 ) );
            }
            else
            {
                return aInterfaces( 0 );
            }
        }

        // Multiple interfaces, create interface manager
        else
        {
            uint tNumCreatedInterfaces = aInterfaces.size();
            aInterfaces.resize( tNumInterfaces );
            for ( uint tInterfaceIndex = tNumCreatedInterfaces; tInterfaceIndex < tNumInterfaces; tInterfaceIndex++ )
            {
                aInterfaces( tInterfaceIndex ) = create_interface( aParameterLists( tInterfaceIndex + 1 ) );
            }
            return std::make_shared< Interface_Manager >( aParameterLists( 0 ), aInterfaces );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Criteria_Interface >
    create_interface( const Parameter_List& aParameterList )
    {
        std::string tInterfaceType = aParameterList.get< std::string >( "field_type" );
        if ( !tInterfaceType.compare( "user_defined" ) )
        {
            return std::make_shared< Interface_User_Defined >( aParameterList );
        }
        else
        {
            MORIS_ERROR( false,
                    "%s is not recognized as a valid Criteria_Interface type in fn_OPT_create_interface.",
                    tInterfaceType.c_str() );
            return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::opt
