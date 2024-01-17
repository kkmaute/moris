/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Design_Factory.cpp
 *
 */

#include "cl_GEN_Design_Factory.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Library_IO.hpp"
#include "fn_GEN_create_field.hpp"

namespace moris::ge
{
    //--------------------------------------------------------------------------------------------------------------

    Design_Factory::Design_Factory(
            Cell< ParameterList >         aParameterLists,
            Matrix< DDRMat >&             aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMesh,
            Node_Manager&                 aNodeManager )
    {
        // Count maximum number of possible designs
        uint tGeometryIndex = 0;
        uint tPropertyIndex = 0;
        uint tFieldIndex    = 0;
        for ( const ParameterList& iParameterList : aParameterLists )
        {
            // Count up the type of design
            if ( iParameterList.exists( "design_type" ) )
            {
                if ( iParameterList.get< std::string >( "design_type" ) == "geometry" )
                {
                    tGeometryIndex++;
                }
                else if ( iParameterList.get< std::string >( "design_type" ) == "property" )
                {
                    tPropertyIndex++;
                }
                else
                {
                    MORIS_ERROR( false, "GEN design parameter list with unknown design type detected." );
                }
            }

            // Count if this design requires a field or not
            if ( iParameterList.exists( "field_type" ) )
            {
                tFieldIndex++;
            }
        }

        // Perform resizes
        mFields.resize( tFieldIndex );
        mGeometries.resize( tGeometryIndex );
        mProperties.resize( tPropertyIndex );
        uint tNumberOfDesignsLeft = tGeometryIndex + tPropertyIndex;

        // Re-initialize counters
        tGeometryIndex = 0;
        tPropertyIndex = 0;
        tFieldIndex    = 0;

        // Loop until all designs are built
        while ( tNumberOfDesignsLeft > 0 )
        {
            // Track if at least one field was found in this loop that can be built, based on its dependencies
            bool tSomethingHasBeenBuilt = false;

            // Loop over all designs
            for ( ParameterList& iParameterList : aParameterLists )
            {
                // If we can build this field or not
                bool tCanBuild = true;

                // Check if a field is required
                if ( iParameterList.exists( "field_type" ) )
                {
                    // Get field dependency names
                    Cell< std::string > tDependencyNames =
                            string_to_cell< std::string >( iParameterList.get< std::string >( "dependencies" ) );

                    // Cell of field dependencies
                    Cell< std::shared_ptr< Field > > tDependencyFields( tDependencyNames.size() );

                    // Loop over dependencies
                    for ( uint iDependencyIndex = 0; iDependencyIndex < tDependencyNames.size(); iDependencyIndex++ )
                    {
                        // Dependency detected, must be found first before building
                        tCanBuild = false;

                        // Loop over fields
                        for ( const auto& iField : mFields )
                        {
                            // Check if field has already been built and name match is found
                            if ( iField and iField->get_name() == tDependencyNames( iDependencyIndex ) )
                            {
                                // Set dependency
                                tDependencyFields( iDependencyIndex ) = iField;
                                tCanBuild                             = true;
                                break;
                            }
                        }

                        // Dependency not found, can't build yet, exit
                        if ( not tCanBuild )
                        {
                            break;
                        }
                    }

                    // Build field if we can
                    if ( tCanBuild )
                    {
                        // Build field
                        mFields( tFieldIndex++ ) = create_field( iParameterList, aADVs, tDependencyFields, aLibrary, aMesh );

                        // Remove this parameter to signal we don't need to build again
                        iParameterList.erase( "field_type" );

                        // Set that a field has been built, so it is fine if another outer loop is needed
                        tSomethingHasBeenBuilt = true;
                    }
                }

                // Check if design needs building
                if ( iParameterList.exists( "design_type" ) )
                {
                    // Get design type
                    std::string tDesignType = iParameterList.get< std::string >( "design_type" );

                    // Geometry
                    if ( tDesignType == "geometry" )
                    {
                        // Base pointer
                        std::shared_ptr< Geometry > tGeometry;

                        // Get geometry type
                        std::string tGeometryType = iParameterList.get< std::string >( "geometry_type" );

                        // Level-set field
                        if ( tGeometryType == "level_set" and tCanBuild )
                        {
                            tGeometry = std::make_shared< Level_Set_Geometry >( mFields( tFieldIndex - 1 ), Level_Set_Parameters( iParameterList ), aNodeManager );
                        }
                        else if ( tGeometryType == "surface_mesh" )
                        {
                            tGeometry              = std::make_shared< Surface_Mesh_Geometry >( Surface_Mesh_Parameters( iParameterList ) );
                            tSomethingHasBeenBuilt = true;
                        }
                        else
                        {
                            MORIS_ERROR( false, "GEN does not recognize geometry type: %s", tGeometryType.c_str() );
                        }

                        // Assign new geometry
                        mGeometries( tGeometryIndex++ ) = tGeometry;
                    }

                    // Property
                    else if ( tDesignType == "property" )
                    {
                        // Create new property
                        auto tProperty = std::make_shared< Property >( mFields( tFieldIndex - 1 ), Property_Parameters( iParameterList ), aNodeManager );

                        // Assign new property
                        mProperties( tPropertyIndex++ ) = tProperty;
                    }
                    else
                    {
                        MORIS_ERROR( false, "GEN does not recognize design type: %s", tDesignType.c_str() );
                    }

                    // Indicate that this design has now been built
                    iParameterList.erase( "design_type" );

                    // Decrement number of designs left
                    tNumberOfDesignsLeft--;
                }
            }

            // Increment and check loop counter
            MORIS_ERROR( tSomethingHasBeenBuilt, "While creating GEN fields, a field dependency was not found or a circular dependency was detected. Exiting." );
        }

        // Resize final containers
        mGeometries.resize( tGeometryIndex );
        mGeometries.shrink_to_fit();
        mProperties.resize( tPropertyIndex );
        mProperties.shrink_to_fit();
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< Geometry > > Design_Factory::get_geometries()
    {
        return mGeometries;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< Property > > Design_Factory::get_properties()
    {
        return mProperties;
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::ge
