/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_geometries.cpp
 *
 */

#include "fn_GEN_create_geometries.hpp"
#include "st_GEN_Geometry_Parameters.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Superellipse.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Superellipsoid.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Geometry.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_GEN_Single_Grain.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Multigeometry.hpp"
#include "cl_GEN_Swiss_Cheese_Slice.hpp"
#include "cl_GEN_Mesh_Field_Geometry.hpp"
#include "cl_GEN_Geometry_SDF.hpp"
#include "cl_GEN_Image_SDF_Geometry.hpp"

namespace moris
{
    namespace ge
    {
        //--------------------------------------------------------------------------------------------------------------
        // Get the length of different vector types
        //--------------------------------------------------------------------------------------------------------------
        inline uint
        get_length( Matrix< DDRMat >& aVector )
        {
            return aVector.length();
        }

        inline uint
        get_length( sol::Dist_Vector* aVector )
        {
            return aVector->vec_local_length();
        }

        //--------------------------------------------------------------------------------------------------------------
        // Definitions
        //--------------------------------------------------------------------------------------------------------------

        template< typename Vector_Type >
        inline Cell< std::shared_ptr< Geometry > >
        create_geometries(
                Cell< ParameterList >         aGeometryParameterLists,
                Vector_Type&                  aADVs,
                std::shared_ptr< Library_IO > aLibrary,
                mtk::Mesh*                    aMTKMesh )
        {
            // Create geometry cell
            Cell< std::shared_ptr< Geometry > >      tGeometries( 0 );
            Cell< std::shared_ptr< Multigeometry > > tMultigeometries( 0 );

            // Create individual geometries
            for ( uint tGeometryIndex = 0; tGeometryIndex < aGeometryParameterLists.size(); tGeometryIndex++ )
            {
                // Create geometry
                std::shared_ptr< Geometry > tGeometry = create_geometry(
                        aGeometryParameterLists( tGeometryIndex ),
                        aADVs,
                        aLibrary,
                        nullptr,
                        aMTKMesh );

                if ( aGeometryParameterLists( tGeometryIndex ).get< bool >( "multilinear_intersections" ) )
                {
                    tGeometry->set_intersection_interpolation( "multilinear" );
                }

                // Determine if to add to multigeometry
                bool        tMultigeometryFound = false;
                std::string tGeometryName       = tGeometry->get_name();
                if ( tGeometryName != "" )
                {
                    // Loop to see if this multigeometry ID exists already
                    for ( uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++ )
                    {
                        if ( tMultigeometries( tMultigeometryIndex )->get_name() == tGeometryName )
                        {
                            tMultigeometryFound = true;
                            tMultigeometries( tMultigeometryIndex )->add_geometry( tGeometry );
                            break;
                        }
                    }

                    // Check for creating new multigeometry
                    if ( not tMultigeometryFound )
                    {
                        for ( uint tCreatedGeometryIndex = 0; tCreatedGeometryIndex < tGeometries.size(); tCreatedGeometryIndex++ )
                        {
                            if ( tGeometries( tCreatedGeometryIndex )->get_name() == tGeometryName )
                            {
                                tMultigeometryFound = true;
                                tMultigeometries.push_back( std::make_shared< Multigeometry >(
                                        Cell< std::shared_ptr< Geometry > >( { tGeometries( tCreatedGeometryIndex ), tGeometry } ) ) );
                                tGeometries.erase( tCreatedGeometryIndex );
                                break;
                            }
                        }
                    }
                }

                // If no multigeometry, add as regular geometry
                if ( not tMultigeometryFound )
                {
                    tGeometries.push_back( tGeometry );
                }

                // TODO generalize this
                if ( aGeometryParameterLists( tGeometryIndex ).get< std::string >( "type" ) == "voxel" )
                {
                    uint tNumVoxelIDs = reinterpret_cast< Voxel_Input* >( tGeometry.get() )->get_num_voxel_Ids();

                    for ( uint tVoxelID = 1; tVoxelID <= tNumVoxelIDs; tVoxelID++ )
                    {
                        std::shared_ptr< Geometry > tGeometrySingleGrain =
                                create_geometry(
                                        aGeometryParameterLists( tGeometryIndex ),
                                        aADVs,
                                        aLibrary,
                                        tGeometry,
                                        aMTKMesh,
                                        tVoxelID );

                        tGeometries.push_back( tGeometrySingleGrain );
                    }
                }
            }

            // Add multigeometries at the end
            for ( uint tMultigeometryIndex = 0; tMultigeometryIndex < tMultigeometries.size(); tMultigeometryIndex++ )
            {
                tGeometries.push_back( tMultigeometries( tMultigeometryIndex ) );
            }

            return tGeometries;
        }

        //--------------------------------------------------------------------------------------------------------------

        template< typename Vector_Type >
        inline std::shared_ptr< Geometry >
        create_geometry(
                ParameterList                 aGeometryParameterList,
                Vector_Type&                  aADVs,
                std::shared_ptr< Library_IO > aLibrary,
                std::shared_ptr< Geometry >   aGeometry,
                mtk::Mesh*                    aMTKMesh,
                uint                          aIndex )
        {
            // Geometry type
            std::string tGeometryType = aGeometryParameterList.get< std::string >( "type" );

            // Geometry inputs
            Matrix< DDUMat > tGeometryVariableIndices( 0, 0 );
            Matrix< DDUMat > tADVIndices( 0, 0 );
            Matrix< DDRMat > tConstants( 0, 0 );

            // If not a swiss cheese, get ADV inputs
            if ( tGeometryType.compare( 0, 12, "swiss_cheese" ) )
            {
                bool tFillVariables = false;
                bool tFillADVs      = false;

                // Determine if variable or ADV indices need to be filled (specified by "all")
                if ( aGeometryParameterList.get< std::string >( "field_variable_indices" ) == "all" )
                {
                    tFillVariables = true;
                }
                else
                {
                    string_to_mat( aGeometryParameterList.get< std::string >( "field_variable_indices" ), tGeometryVariableIndices );
                }
                if ( aGeometryParameterList.get< std::string >( "adv_indices" ) == "all" )
                {
                    tFillADVs = true;
                }
                else
                {
                    string_to_mat( aGeometryParameterList.get< std::string >( "adv_indices" ), tADVIndices );
                }

                // Perform fill
                if ( tFillVariables and tFillADVs )
                {
                    uint tNumADVs = get_length( aADVs );
                    tGeometryVariableIndices.resize( tNumADVs, 1 );
                    tADVIndices.resize( tNumADVs, 1 );
                    for ( uint tIndex = 0; tIndex < tNumADVs; tIndex++ )
                    {
                        tGeometryVariableIndices( tIndex ) = tIndex;
                        tADVIndices( tIndex )              = tIndex;
                    }
                }
                else if ( tFillVariables )
                {
                    tGeometryVariableIndices.resize( tADVIndices.length(), 1 );
                    for ( uint tIndex = 0; tIndex < tADVIndices.length(); tIndex++ )
                    {
                        tGeometryVariableIndices( tIndex ) = tIndex;
                    }
                }
                else if ( tFillADVs )
                {
                    tADVIndices.resize( tGeometryVariableIndices.length(), 1 );
                    for ( uint tIndex = 0; tIndex < tGeometryVariableIndices.length(); tIndex++ )
                    {
                        tADVIndices( tIndex ) = tIndex;
                    }
                }

                // Constant parameters
                tConstants = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "constant_parameters" ) );
            }

            // Geometry parameters
            Geometry_Field_Parameters tParameters;
            tParameters.mName                     = aGeometryParameterList.get< std::string >( "name" );
            tParameters.mNumRefinements           = aGeometryParameterList.get< std::string >( "number_of_refinements" );
            tParameters.mRefinementMeshIndices    = aGeometryParameterList.get< std::string >( "refinement_mesh_index" );
            tParameters.mRefinementFunctionIndex  = aGeometryParameterList.get< sint >( "refinement_function_index" );
            tParameters.mDiscretizationMeshIndex  = aGeometryParameterList.get< sint >( "discretization_mesh_index" );
            tParameters.mDiscretizationLowerBound = aGeometryParameterList.get< real >( "discretization_lower_bound" );
            tParameters.mDiscretizationUpperBound = aGeometryParameterList.get< real >( "discretization_upper_bound" );
            tParameters.mIsocontourThreshold      = aGeometryParameterList.get< real >( "isocontour_threshold" );
            tParameters.mIsocontourTolerance      = aGeometryParameterList.get< real >( "isocontour_tolerance" );
            tParameters.mIntersectionTolerance    = aGeometryParameterList.get< real >( "intersection_tolerance" );

            // Build Geometry
            if ( tGeometryType == "circle" )
            {
                return std::make_shared< Circle >( aADVs, tGeometryVariableIndices, tADVIndices, tConstants, tParameters );
            }
            else if ( tGeometryType == "superellipse" )
            {
                return std::make_shared< Superellipse >( aADVs, tGeometryVariableIndices, tADVIndices, tConstants, tParameters );
            }
            else if ( tGeometryType == "sphere" )
            {
                return std::make_shared< Sphere >( aADVs, tGeometryVariableIndices, tADVIndices, tConstants, tParameters );
            }
            else if ( tGeometryType == "superellipsoid" )
            {
                return std::make_shared< Superellipsoid >( aADVs, tGeometryVariableIndices, tADVIndices, tConstants, tParameters );
            }
            else if ( tGeometryType == "plane" )
            {
                return std::make_shared< Plane >( aADVs, tGeometryVariableIndices, tADVIndices, tConstants, tParameters );
            }
            else if ( tGeometryType == "nodal_field" )
            {
                MORIS_ERROR( aMTKMesh != nullptr, "Mesh is a null ptr for nodal field geometry" );

                return std::make_shared< ge::Mesh_Field_Geometry >( aMTKMesh, tParameters.mName, mtk::EntityRank::NODE );
            }
            else if ( tGeometryType == "nodal_field_from_file" )
            {
                std::string tFileName    = aGeometryParameterList.get< std::string >( "file_name" );
                std::string tFieldName   = aGeometryParameterList.get< std::string >( "field_name" );
                std::string tFieldFormat = aGeometryParameterList.get< std::string >( "file_format" );
                real        tOffset      = aGeometryParameterList.get< real >( "offset" );

                return std::make_shared< ge::Mesh_Field_Geometry >(
                        aMTKMesh,
                        tFileName,
                        tFieldName,
                        tFieldFormat,
                        tOffset,
                        mtk::EntityRank::NODE,
                        tParameters );
            }
            else if ( tGeometryType == "sdf_field" )
            {
                std::string      tObjectPath   = aGeometryParameterList.get< std::string >( "sdf_object_path" );
                Matrix< DDRMat > tObjectOffset = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "sdf_object_offset" ) );
                real             tSDFShift     = aGeometryParameterList.get< real >( "sdf_shift" );

                return std::make_shared< ge::Geometry_SDF >(
                        tObjectPath,
                        tObjectOffset,
                        tSDFShift,
                        tParameters );
            }
            else if ( tGeometryType == "image_sdf" )
            {
                // Get voxel-specific info
                std::string      tImageFileName    = aGeometryParameterList.get< std::string >( "image_file" );
                Matrix< DDRMat > tDomainDimensions = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "image_dimensions" ) );
                Matrix< DDRMat > tDomainOffset     = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "image_offset" ) );

                real tSDFScaling = aGeometryParameterList.get< real >( "image_sdf_scaling" );
                real tSDFShift   = aGeometryParameterList.get< real >( "image_sdf_shift" );
                real tSDFDefault = aGeometryParameterList.get< real >( "image_sdf_default" );
                bool tsDFInterp  = aGeometryParameterList.get< bool >( "image_sdf_interpolate" );

                return std::make_shared< Image_SDF_Geometry >(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstants,
                        tImageFileName,
                        tDomainDimensions,
                        tDomainOffset,
                        tSDFScaling,
                        tSDFShift,
                        tSDFDefault,
                        tsDFInterp,
                        tParameters );
            }
            else if ( tGeometryType == "user_defined" )
            {

                // Get intersection mode
                std::string tIntersectionModeString = aGeometryParameterList.get< std::string >( "intersection_mode" );

                map< std::string, Intersection_Mode > tIntersectionModeMap = get_intersection_mode_map();

                tParameters.mIntersectionMode = tIntersectionModeMap[ tIntersectionModeString ];

                // Check if library is given
                MORIS_ERROR( aLibrary != nullptr, "Library must be given in order to create a user-defined geometry." );

                // Get sensitivity function if needed
                std::string          tSensitivityFunctionName = aGeometryParameterList.get< std::string >( "sensitivity_function_name" );
                Sensitivity_Function tSensitivityFunction =
                        ( tSensitivityFunctionName == "" ? nullptr : aLibrary->load_function< Sensitivity_Function >( tSensitivityFunctionName ) );

                // Create user-defined geometry
                return std::make_shared< User_Defined_Geometry >(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstants,
                        aLibrary->load_function< Field_Function >( aGeometryParameterList.get< std::string >( "field_function_name" ) ),
                        tSensitivityFunction,
                        tParameters );
            }
            else if ( tGeometryType == "voxel" and aGeometry )
            {
                return std::make_shared< Single_Grain >(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstants,
                        aGeometry,
                        aIndex,
                        tParameters );
            }
            else if ( tGeometryType == "voxel" )
            {
                // Get voxel-specific info
                std::string      tVoxelFieldName    = aGeometryParameterList.get< std::string >( "voxel_field_file" );
                Matrix< DDRMat > tDomainDimensions  = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "domain_dimensions" ) );
                Matrix< DDRMat > tDomainOffset      = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "domain_offset" ) );
                Matrix< DDRMat > tGrainIdToValueMap = string_to_mat< DDRMat >( aGeometryParameterList.get< std::string >( "grain_id_value_map" ) );

                return std::make_shared< Voxel_Input >(
                        aADVs,
                        tGeometryVariableIndices,
                        tADVIndices,
                        tConstants,
                        tVoxelFieldName,
                        tDomainDimensions,
                        tDomainOffset,
                        tGrainIdToValueMap,
                        tParameters );
            }
            else if ( tGeometryType == "swiss_cheese_slice" )
            {
                // Check for definition
                uint tNumXHoles      = (uint)aGeometryParameterList.get< sint >( "number_of_x_holes" );
                uint tNumYHoles      = (uint)aGeometryParameterList.get< sint >( "number_of_y_holes" );
                real tTargetXSpacing = aGeometryParameterList.get< real >( "target_x_spacing" );
                real tTargetYSpacing = aGeometryParameterList.get< real >( "target_y_spacing" );

                MORIS_ERROR( ( tNumXHoles > 1 && tNumYHoles > 1 ) || ( tTargetXSpacing && tTargetYSpacing ),
                        "In a swiss cheese parameter list, you must specify either a number of holes > 1 %s",
                        "or a target spacing in each direction.\n" );

                if ( tNumXHoles )
                {
                    return std::make_shared< Swiss_Cheese_Slice >(
                            aGeometryParameterList.get< real >( "left_bound" ),
                            aGeometryParameterList.get< real >( "right_bound" ),
                            aGeometryParameterList.get< real >( "bottom_bound" ),
                            aGeometryParameterList.get< real >( "top_bound" ),
                            tNumXHoles,
                            tNumYHoles,
                            aGeometryParameterList.get< real >( "hole_x_semidiameter" ),
                            aGeometryParameterList.get< real >( "hole_y_semidiameter" ),
                            aGeometryParameterList.get< real >( "superellipse_exponent" ),
                            aGeometryParameterList.get< real >( "superellipse_scaling" ),
                            aGeometryParameterList.get< real >( "superellipse_regularization" ),
                            aGeometryParameterList.get< real >( "superellipse_shift" ),
                            aGeometryParameterList.get< real >( "row_offset" ),
                            tParameters );
                }
                else
                {
                    return std::make_shared< Swiss_Cheese_Slice >(
                            aGeometryParameterList.get< real >( "left_bound" ),
                            aGeometryParameterList.get< real >( "right_bound" ),
                            aGeometryParameterList.get< real >( "bottom_bound" ),
                            aGeometryParameterList.get< real >( "top_bound" ),
                            tTargetXSpacing,
                            tTargetYSpacing,
                            aGeometryParameterList.get< real >( "hole_x_semidiameter" ),
                            aGeometryParameterList.get< real >( "hole_y_semidiameter" ),
                            aGeometryParameterList.get< real >( "superellipse_exponent" ),
                            aGeometryParameterList.get< real >( "superellipse_scaling" ),
                            aGeometryParameterList.get< real >( "superellipse_regularization" ),
                            aGeometryParameterList.get< real >( "superellipse_shift" ),
                            aGeometryParameterList.get< real >( "row_offset" ),
                            aGeometryParameterList.get< bool >( "allow_less_than_target_spacing" ),
                            tParameters );
                }
            }
            else
            {
                MORIS_ERROR( false, "%s is not recognized as a valid Geometry type in fn_GEN_create_geometry.", tGeometryType.c_str() );
                return nullptr;
            }
        }

        //--------------------------------------------------------------------------------------------------------------
        // Explicit template instantiation
        //--------------------------------------------------------------------------------------------------------------

        template Cell< std::shared_ptr< Geometry > > create_geometries(
                Cell< ParameterList >         aGeometryParameterLists,
                Matrix< DDRMat >&             aADVs,
                std::shared_ptr< Library_IO > aLibrary,
                mtk::Mesh*                    aMTKMesh );

        template Cell< std::shared_ptr< Geometry > > create_geometries(
                Cell< ParameterList >         aGeometryParameterLists,
                sol::Dist_Vector*&            aADVs,
                std::shared_ptr< Library_IO > aLibrary,
                mtk::Mesh*                    aMTKMesh );

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace ge
}    // namespace moris
