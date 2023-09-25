/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_properties.cpp
 *
 */

#include "fn_GEN_create_properties.hpp"
#include "fn_Parsing_Tools.hpp"

#include "fn_GEN_create_field.hpp"

namespace moris::ge
{
    template< typename Vector_Type >
    Cell< std::shared_ptr< Property > >
    create_properties(
            Cell< ParameterList >               aPropertyParameterLists,
            Vector_Type&                        aADVs,
            Cell< std::shared_ptr< Level_Set_Geometry > > aGeometries,
            std::shared_ptr< Library_IO >       aLibrary )
    {
        // Initialize
        uint                                     tNumProperties = aPropertyParameterLists.size();
        Cell< std::shared_ptr< Property > >      tProperties( tNumProperties );
        Cell< std::string >                      tPropertyNames( tNumProperties );
        Cell< Cell< std::string > >              tNeededFieldNames( tNumProperties );
        Cell< Cell< std::shared_ptr< Field > > > tNeededFields( tNumProperties );

        // Fill names, dependencies
        for ( uint tPropertyIndex = 0; tPropertyIndex < tNumProperties; tPropertyIndex++ )
        {
            tPropertyNames( tPropertyIndex ) = aPropertyParameterLists( tPropertyIndex ).get< std::string >( "name" );
            tNeededFieldNames( tPropertyIndex ) =
                    string_to_cell< std::string >( aPropertyParameterLists( tPropertyIndex ).get< std::string >( "dependencies" ) );
            tNeededFields( tPropertyIndex ).resize( tNeededFieldNames( tPropertyIndex ).size() );
        }

        // Build based on dependencies (this is not optimally efficient, but doesn't need to be)
        bool tBuild;
        uint tNumPropertiesLeft = tNumProperties;
        uint tLoopCount         = 0;
        while ( tNumPropertiesLeft > 0 )
        {
            for ( uint tBuildPropertyIndex = 0; tBuildPropertyIndex < tNumProperties; tBuildPropertyIndex++ )
            {
                // Check if property needs to be built
                if ( tProperties( tBuildPropertyIndex ) == nullptr )
                {
                    tBuild = true;

                    // Check if property dependencies are built
                    if ( tNeededFieldNames( tBuildPropertyIndex ).size() > 0 )
                    {
                        // Loop over dependencies
                        for ( uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames( tBuildPropertyIndex ).size(); tDependencyIndex++ )
                        {
                            // Checking each property by name
                            for ( uint tCheckPropertyIndex = 0; tCheckPropertyIndex < tNumProperties; tCheckPropertyIndex++ )
                            {
                                // Name match found
                                if ( tNeededFieldNames( tBuildPropertyIndex )( tDependencyIndex ) == tPropertyNames( tCheckPropertyIndex ) )
                                {
                                    // Dependency is built already
                                    if ( tProperties( tCheckPropertyIndex ) )
                                    {
                                        tNeededFields( tBuildPropertyIndex )( tDependencyIndex ) = tProperties( tCheckPropertyIndex )->get_field();
                                    }

                                    // Dependency is not built, cannot build current property
                                    else
                                    {
                                        tBuild = false;
                                    }
                                }
                            }
                        }
                    }

                    // Build
                    if ( tBuild )
                    {
                        // Loop over dependencies, this time to add generic field dependencies
                        for ( uint tDependencyIndex = 0; tDependencyIndex < tNeededFieldNames( tBuildPropertyIndex ).size(); tDependencyIndex++ )
                        {
                            // Checking each field by name
                            for ( uint tCheckFieldIndex = 0; tCheckFieldIndex < aGeometries.size(); tCheckFieldIndex++ )
                            {
                                // Name match found
                                if ( tNeededFieldNames( tBuildPropertyIndex )( tDependencyIndex ) == aGeometries( tCheckFieldIndex )->get_name() )
                                {
                                    tNeededFields( tBuildPropertyIndex )( tDependencyIndex ) = aGeometries( tCheckFieldIndex )->get_field();
                                }
                            }
                        }

                        // Build property and decrement remaining properties to build
                        tProperties( tBuildPropertyIndex ) = create_property(
                                aPropertyParameterLists( tBuildPropertyIndex ),
                                aADVs,
                                tNeededFields( tBuildPropertyIndex ),
                                aLibrary );
                        tNumPropertiesLeft--;
                    }
                }
            }
            tLoopCount++;
            MORIS_ERROR( tLoopCount <= tNumProperties, "While creating GEN properties, a circular property dependency was detected. Exiting." );
        }

        return tProperties;
    }

    //--------------------------------------------------------------------------------------------------------------

    template< typename Vector_Type >
    std::shared_ptr< Property >
    create_property(
            ParameterList                    aPropertyParameterList,
            Vector_Type&                     aADVs,
            Cell< std::shared_ptr< Field > > aFieldDependencies,
            std::shared_ptr< Library_IO >    aLibrary )
    {
        // Property type/name
        std::string tPropertyType = aPropertyParameterList.get< std::string >( "type" );

        // Create field
        std::shared_ptr< Field > tField = create_field( aPropertyParameterList, aADVs, aFieldDependencies, aLibrary );

        // Property parameters
        Property_Parameters tParameters;
        tParameters.mName                     = aPropertyParameterList.get< std::string >( "name" );
        tParameters.mNumRefinements           = aPropertyParameterList.get< std::string >( "number_of_refinements" );
        tParameters.mRefinementMeshIndices    = aPropertyParameterList.get< std::string >( "refinement_mesh_index" );
        tParameters.mRefinementFunctionIndex  = aPropertyParameterList.get< sint >( "refinement_function_index" );
        tParameters.mDiscretizationIndex      = aPropertyParameterList.get< sint >( "discretization_mesh_index" );
        tParameters.mDiscretizationLowerBound = aPropertyParameterList.get< real >( "discretization_lower_bound" );
        tParameters.mDiscretizationUpperBound = aPropertyParameterList.get< real >( "discretization_upper_bound" );

        map< std::string, PDV_Type > tPDVTypeMap = get_pdv_type_map();
        tParameters.mDependencyNames             = string_to_cell< std::string >( aPropertyParameterList.get< std::string >( "dependencies" ) );
        tParameters.mPDVType                     = tPDVTypeMap[ aPropertyParameterList.get< std::string >( "pdv_type" ) ];
        tParameters.mPDVMeshSetIndices           = string_to_mat< DDUMat >( aPropertyParameterList.get< std::string >( "pdv_mesh_set_indices" ) );
        tParameters.mPDVMeshSetNames             = string_to_cell< std::string >( aPropertyParameterList.get< std::string >( "pdv_mesh_set_names" ) );

        // Build Property
        return std::make_shared< Property >( tField, tParameters );
    }

    //--------------------------------------------------------------------------------------------------------------
    // Explicit template instantiation
    //--------------------------------------------------------------------------------------------------------------

    template Cell< std::shared_ptr< Property > > create_properties(
            Cell< ParameterList >               aPropertyParameterLists,
            Matrix< DDRMat >&                   aADVs,
            Cell< std::shared_ptr< Level_Set_Geometry > > aGeometries,
            std::shared_ptr< Library_IO >       aLibrary );

    template Cell< std::shared_ptr< Property > > create_properties(
            Cell< ParameterList >               aPropertyParameterLists,
            sol::Dist_Vector*&                  aADVs,
            Cell< std::shared_ptr< Level_Set_Geometry > > aGeometries,
            std::shared_ptr< Library_IO >       aLibrary );

    //--------------------------------------------------------------------------------------------------------------

}
