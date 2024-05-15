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
#include "cl_Parameter_List.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_Library_IO.hpp"
#include "fn_GEN_create_field.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_GEN_Voxel_Geometry.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Design_Factory::Design_Factory(
            Vector< Parameter_List >      aParameterLists,
            Matrix< DDRMat >&             aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMesh,
            Node_Manager&                 aNodeManager )
    {
        // Count maximum number of possible designs
        uint tGeometryIndex = 0;
        uint tPropertyIndex = 0;
        uint tFieldIndex    = 0;
        for ( const Parameter_List& iParameterList : aParameterLists )
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
            for ( Parameter_List& iParameterList : aParameterLists )
            {
                // If we can build this field or not
                bool tCanBuild = true;

                // Check if a field is required
                if ( iParameterList.exists( "field_type" ) )
                {
                    // Get field dependency names
                    Vector< std::string > tDependencyNames =
                            string_to_cell< std::string >( iParameterList.get< std::string >( "dependencies" ) );

                    // Cell of field dependencies
                    Vector< std::shared_ptr< Field > > tDependencyFields( tDependencyNames.size() );

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
                        // Get geometry type
                        std::string tGeometryType = iParameterList.get< std::string >( "geometry_type" );

                        // Level-set field
                        if ( tGeometryType == "level_set" and tCanBuild )
                        {
                            mGeometries( tGeometryIndex++ ) = std::make_shared< Level_Set_Geometry >( mFields( tFieldIndex - 1 ), Level_Set_Parameters( iParameterList ), aNodeManager );
                        }
                        else if ( tGeometryType == "voxel" )
                        {
                            // Get voxel-specific info
                            std::string      tVoxelFieldName    = iParameterList.get< std::string >( "voxel_field_file" );
                            Matrix< DDRMat > tDomainDimensions  = string_to_mat< DDRMat >( iParameterList.get< std::string >( "domain_dimensions" ) );
                            Matrix< DDRMat > tDomainOffset      = string_to_mat< DDRMat >( iParameterList.get< std::string >( "domain_offset" ) );
                            Matrix< DDRMat > tGrainIdToValueMap = string_to_mat< DDRMat >( iParameterList.get< std::string >( "grain_id_value_map" ) );

                            // Create voxel input
                            auto tVoxelInput = std::make_shared< Voxel_Input >(
                                    tVoxelFieldName,
                                    tDomainDimensions,
                                    tDomainOffset,
                                    tGrainIdToValueMap,
                                    aNodeManager );

                            // Get number of voxel IDs
                            uint tNumberOfVoxelIDs = tVoxelInput->get_num_voxel_IDs();

                            // Resize geometries
                            mGeometries.resize( mGeometries.size() + tNumberOfVoxelIDs );

                            // Create each single voxel geometry
                            for ( uint iVoxelIndex = 0; iVoxelIndex < tNumberOfVoxelIDs; iVoxelIndex++ )
                            {
                                mGeometries( tGeometryIndex++ ) = std::make_shared< Voxel_Geometry >( tVoxelInput, iVoxelIndex );
                            }

                            // Let the factory know it is still working
                            tSomethingHasBeenBuilt = true;
                        }
                        else if ( tGeometryType == "surface_mesh" )
                        {
                            mGeometries( tGeometryIndex++ ) = std::make_shared< Surface_Mesh_Geometry >(
                                    aMesh,
                                    aADVs,
                                    Surface_Mesh_Parameters( iParameterList ),
                                    aNodeManager,
                                    aLibrary );
                            tSomethingHasBeenBuilt = true;
                        }
                        else
                        {
                            MORIS_ERROR( false, "GEN does not recognize geometry type: %s", tGeometryType.c_str() );
                        }
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

    Vector< std::shared_ptr< Geometry > > Design_Factory::get_geometries()
    {
        return mGeometries;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< std::shared_ptr< Property > > Design_Factory::get_properties()
    {
        return mProperties;
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::gen
