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
            const Cell< ParameterList >&  aGeometryParameterLists,
            const Cell< ParameterList >&  aPropertyParameterLists,
            Matrix< DDRMat >&             aADVs,
            std::shared_ptr< Library_IO > aLibrary,
            mtk::Mesh*                    aMesh )
            : mDesigns( aGeometryParameterLists.size() + aPropertyParameterLists.size(), nullptr )
            , mGeometries( aGeometryParameterLists.size(), nullptr )
            , mProperties( aPropertyParameterLists.size(), nullptr )
    {
        // Combine parameter lists
        Cell< ParameterList > tDesignParameterLists;
        tDesignParameterLists.append( aGeometryParameterLists );
        tDesignParameterLists.append( aPropertyParameterLists );

        // Initialize counters
        uint tGeometryIndex = 0;
        uint tPropertyIndex = 0;
        uint tDesignIndex = 0;
        uint tLoopCount = 0;
        uint tNumberOfDesignsLeft = tDesignParameterLists.size();

        // Loop until all designs are built
        while ( tNumberOfDesignsLeft > 0 )
        {
            // Loop over all designs
            for ( ParameterList& iDesignParameterList : tDesignParameterLists )
            {
                // Check if design needs building
                if ( iDesignParameterList.exists( "design_type" ) )
                {
                    // Check if a field is required
                    if ( iDesignParameterList.exists( "field_type" ) )
                    {
                        // Get dependency names
                        Cell< std::string > tDependencyNames =
                                string_to_cell< std::string >( iDesignParameterList.get< std::string >( "dependencies" ) );

                        // Cell of field dependencies
                        Cell< std::shared_ptr< Field > > tDependencyFields( tDependencyNames.size() );

                        // If we can build this field or not
                        bool tCanBuild = true;

                        // Loop over dependencies
                        for ( uint iDependencyIndex = 0; iDependencyIndex < tDependencyNames.size(); iDependencyIndex++ )
                        {
                            // Dependency detected, must be found first before building
                            tCanBuild = false;

                            // Loop over designs
                            for ( const std::shared_ptr< Design_Field >& iDesign : mDesigns )
                            {
                                // Check if design has already been built and name match is found
                                if ( iDesign and iDesign->get_name() == tDependencyNames( iDependencyIndex ) )
                                {
                                    // Set dependency
                                    tDependencyFields( iDependencyIndex ) = iDesign->get_field();
                                    tCanBuild = true;
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
                            // Field object
                            std::shared_ptr< Field > tField = create_field( iDesignParameterList, aADVs, tDependencyFields, aLibrary, aMesh );

                            // Get design type
                            std::string tDesignType = iDesignParameterList.get< std::string >( "design_type" );

                            // Geometry
                            if ( tDesignType == "geometry" )
                            {
                                // Base pointer
                                std::shared_ptr< Level_Set_Geometry > tGeometry;

                                // Get geometry type
                                std::string tGeometryType = iDesignParameterList.get< std::string >( "geometry_type" );

                                // Level-set field
                                if ( tGeometryType == "level_set" )
                                {
                                    tGeometry = std::make_shared< Level_Set_Geometry >( tField, Level_Set_Parameters( iDesignParameterList ) );
                                }
                                else
                                {
                                    MORIS_ERROR( false, "GEN does not recognize geometry type: %s", tGeometryType.c_str() );
                                }

                                // Assign new geometry
                                mGeometries( tGeometryIndex++ ) = tGeometry;
                                mDesigns( tDesignIndex++ ) = tGeometry;
                            }

                            // Property
                            else if ( tDesignType == "property" )
                            {
                                // Create new property
                                auto tProperty = std::make_shared< Property >( tField, Property_Parameters( iDesignParameterList ) );

                                // Assign new property
                                mProperties( tPropertyIndex++ ) = tProperty;
                                mDesigns( tDesignIndex++ ) = tProperty;
                            }
                            else
                            {
                                MORIS_ERROR( false, "GEN does not recognize design type: %s", tDesignType.c_str() );
                            }

                            // Indicate that this design has now been built
                            iDesignParameterList.erase( "design_type" );

                            // Decrement number of designs left
                            tNumberOfDesignsLeft--;
                        }
                    }
                }
            }

            // Increment and check loop counter
            MORIS_ERROR( tLoopCount++ < tDesignParameterLists.size(), "While creating GEN fields, a field dependency was not found or a circular dependency was detected. Exiting." );
        }

        // Resize final containers
        mGeometries.resize( tGeometryIndex );
        mGeometries.shrink_to_fit();
        mProperties.resize( tPropertyIndex );
        mProperties.shrink_to_fit();
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< Level_Set_Geometry > > Design_Factory::get_geometries()
    {
        return mGeometries;
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< std::shared_ptr< Property > > Design_Factory::get_properties()
    {
        return mProperties;
    }

    //--------------------------------------------------------------------------------------------------------------
}
