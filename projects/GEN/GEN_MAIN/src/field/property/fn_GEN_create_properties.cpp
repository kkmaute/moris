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
#include "st_GEN_Property_Parameters.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Constant_Property.hpp"
#include "cl_GEN_Scaled_Field.hpp"
#include "cl_GEN_User_Defined_Property.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------
        // Get the length of different vector types
        //--------------------------------------------------------------------------------------------------------------

        uint
        get_length( Matrix< DDRMat >& aVector )
        {
            return aVector.length();
        }

        uint
        get_length( sol::Dist_Vector* aVector )
        {
            return aVector->vec_local_length();
        }

        //--------------------------------------------------------------------------------------------------------------
        // Definitions
        //--------------------------------------------------------------------------------------------------------------

        template< typename Vector_Type >
        Vector< std::shared_ptr< Property > >
        create_properties(
                Vector< ParameterList >               aPropertyParameterLists,
                Vector_Type&                        aADVs,
                Vector< std::shared_ptr< Geometry > > aGeometries,
                std::shared_ptr< Library_IO >       aLibrary )
        {
            // Initialize
            uint                                     tNumProperties = aPropertyParameterLists.size();
            Vector< std::shared_ptr< Property > >      tProperties( tNumProperties );
            Vector< std::string >                      tPropertyNames( tNumProperties );
            Vector< Vector< std::string > >              tNeededFieldNames( tNumProperties );
            Vector< Vector< std::shared_ptr< Field > > > tNeededFields( tNumProperties );

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
                                            tNeededFields( tBuildPropertyIndex )( tDependencyIndex ) = tProperties( tCheckPropertyIndex );
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
                                        tNeededFields( tBuildPropertyIndex )( tDependencyIndex ) = aGeometries( tCheckFieldIndex );
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
                MORIS_ERROR( tLoopCount <= tNumProperties, "In fn_GEN_create_properties, a circular property dependency was detected. Exiting." );
            }

            return tProperties;
        }

        //--------------------------------------------------------------------------------------------------------------

        template< typename Vector_Type >
        std::shared_ptr< Property >
        create_property(
                ParameterList                    aPropertyParameterList,
                Vector_Type&                     aADVs,
                Vector< std::shared_ptr< Field > > aFieldDependencies,
                std::shared_ptr< Library_IO >    aLibrary )
        {
            // Property type/name
            std::string tPropertyType = aPropertyParameterList.get< std::string >( "type" );

            // Property inputs
            Matrix< DDUMat > tPropertyVariableIndices( 0, 0 );
            Matrix< DDUMat > tADVIndices( 0, 0 );
            bool             tFillVariables = false;
            bool             tFillADVs      = false;

            // Determine if variable or ADV indices need to be filled (specified by "all")
            if ( aPropertyParameterList.get< std::string >( "field_variable_indices" ) == "all" )
            {
                tFillVariables = true;
            }
            else
            {
                string_to_mat( aPropertyParameterList.get< std::string >( "field_variable_indices" ), tPropertyVariableIndices );
            }
            if ( aPropertyParameterList.get< std::string >( "adv_indices" ) == "all" )
            {
                tFillADVs = true;
            }
            else
            {
                string_to_mat( aPropertyParameterList.get< std::string >( "adv_indices" ), tADVIndices );
            }

            // Perform fill
            if ( tFillVariables and tFillADVs )
            {
                uint tNumADVs = get_length( aADVs );
                tPropertyVariableIndices.resize( tNumADVs, 1 );
                tADVIndices.resize( tNumADVs, 1 );
                for ( uint tIndex = 0; tIndex < tNumADVs; tIndex++ )
                {
                    tPropertyVariableIndices( tIndex ) = tIndex;
                    tADVIndices( tIndex )              = tIndex;
                }
            }
            else if ( tFillVariables )
            {
                tPropertyVariableIndices.resize( tADVIndices.length(), 1 );
                for ( uint tIndex = 0; tIndex < tADVIndices.length(); tIndex++ )
                {
                    tPropertyVariableIndices( tIndex ) = tIndex;
                }
            }
            else if ( tFillADVs )
            {
                tADVIndices.resize( tPropertyVariableIndices.length(), 1 );
                for ( uint tIndex = 0; tIndex < tPropertyVariableIndices.length(); tIndex++ )
                {
                    tADVIndices( tIndex ) = tIndex;
                }
            }

            // Get constant parameters
            Matrix< DDRMat > tConstants = string_to_mat< DDRMat >( aPropertyParameterList.get< std::string >( "constant_parameters" ) );

            // Property parameters
            Property_Field_Parameters tParameters;
            tParameters.mName                     = aPropertyParameterList.get< std::string >( "name" );
            tParameters.mNumRefinements           = aPropertyParameterList.get< std::string >( "number_of_refinements" );
            tParameters.mRefinementMeshIndices    = aPropertyParameterList.get< std::string >( "refinement_mesh_index" );
            tParameters.mRefinementFunctionIndex  = aPropertyParameterList.get< sint >( "refinement_function_index" );
            tParameters.mDiscretizationMeshIndex  = aPropertyParameterList.get< sint >( "discretization_mesh_index" );
            tParameters.mDiscretizationLowerBound = aPropertyParameterList.get< real >( "discretization_lower_bound" );
            tParameters.mDiscretizationUpperBound = aPropertyParameterList.get< real >( "discretization_upper_bound" );

            map< std::string, PDV_Type > tPDVTypeMap = get_pdv_type_map();
            tParameters.mDependencyNames             = string_to_cell< std::string >( aPropertyParameterList.get< std::string >( "dependencies" ) );
            tParameters.mPDVType                     = tPDVTypeMap[ aPropertyParameterList.get< std::string >( "pdv_type" ) ];
            tParameters.mPDVMeshSetIndices           = string_to_mat< DDUMat >( aPropertyParameterList.get< std::string >( "pdv_mesh_set_indices" ) );
            tParameters.mPDVMeshSetNames             = string_to_cell< std::string >( aPropertyParameterList.get< std::string >( "pdv_mesh_set_names" ) );

            // Build Property
            if ( tPropertyType == "constant" )
            {
                return std::make_shared< Constant_Property >(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstants,
                        tParameters );
            }
            else if ( tPropertyType == "scaled_field" )
            {
                return std::make_shared< Scaled_Field >(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstants,
                        aFieldDependencies( 0 ),
                        tParameters );
            }
            else if ( tPropertyType == "user_defined" )
            {
                // Check if library is given
                MORIS_ERROR( aLibrary != nullptr, "Library must be given in order to create a user-defined property." );

                // Get sensitivity function if needed
                std::string          tSensitivityFunctionName = aPropertyParameterList.get< std::string >( "sensitivity_function_name" );
                Sensitivity_Function tSensitivityFunction =
                        ( tSensitivityFunctionName == "" ? nullptr : aLibrary->load_function< Sensitivity_Function >( tSensitivityFunctionName ) );

                return std::make_shared< User_Defined_Property >(
                        aADVs,
                        tPropertyVariableIndices,
                        tADVIndices,
                        tConstants,
                        aLibrary->load_function< Field_Function >( aPropertyParameterList.get< std::string >( "field_function_name" ) ),
                        tSensitivityFunction,
                        tParameters );
            }
            else
            {
                MORIS_ERROR( false, "%s is not recognized as a valid Property type in fn_GEN_create_property.", tPropertyType.c_str() );
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        // Explicit template instantiation
        //--------------------------------------------------------------------------------------------------------------

        template Vector< std::shared_ptr< Property > > create_properties(
                Vector< ParameterList >               aPropertyParameterLists,
                Matrix< DDRMat >&                   aADVs,
                Vector< std::shared_ptr< Geometry > > aGeometries,
                std::shared_ptr< Library_IO >       aLibrary );

        template Vector< std::shared_ptr< Property > > create_properties(
                Vector< ParameterList >               aPropertyParameterLists,
                sol::Dist_Vector*&                  aADVs,
                Vector< std::shared_ptr< Geometry > > aGeometries,
                std::shared_ptr< Library_IO >       aLibrary );

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
