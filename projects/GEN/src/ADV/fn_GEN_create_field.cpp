/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_field.cpp
 *
 */

#include "fn_GEN_create_field.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Constant_Property.hpp"
#include "cl_GEN_Scaled_Field.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Superellipse.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Superellipsoid.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_GEN_Single_Grain.hpp"
#include "cl_GEN_Combined_Field.hpp"
#include "cl_GEN_Swiss_Cheese_Slice.hpp"
#include "cl_GEN_Mesh_Field_Geometry.hpp"
#include "cl_GEN_Geometry_SDF.hpp"
#include "cl_GEN_Image_SDF_Geometry.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< Field >
    create_field(
            ParameterList                    aFieldParameterList,
            Matrix< DDRMat >&                aADVs,
            Cell< std::shared_ptr< Field > > aFieldDependencies,
            std::shared_ptr< Library_IO >    aLibrary,
            mtk::Mesh*                       aMTKMesh,
            uint                             aIndex )
    {
        // Geometry type
        std::string tFieldType = aFieldParameterList.get< std::string >( "field_type" );

        // Geometry inputs
        Matrix< DDUMat > tVariableIndices( 0, 0 );
        Matrix< DDUMat > tADVIndices( 0, 0 );
        Matrix< DDRMat > tConstants( 0, 0 );

        // If not a swiss cheese, get ADV inputs
        if ( tFieldType.compare( 0, 12, "swiss_cheese" ) )
        {
            bool tFillVariables = false;
            bool tFillADVs      = false;

            // Determine if variable or ADV indices need to be filled (specified by "all")
            if ( aFieldParameterList.get< std::string >( "field_variable_indices" ) == "all" )
            {
                tFillVariables = true;
            }
            else
            {
                string_to_mat( aFieldParameterList.get< std::string >( "field_variable_indices" ), tVariableIndices );
            }
            if ( aFieldParameterList.get< std::string >( "adv_indices" ) == "all" )
            {
                tFillADVs = true;
            }
            else
            {
                string_to_mat( aFieldParameterList.get< std::string >( "adv_indices" ), tADVIndices );
            }

            // Perform fill
            if ( tFillVariables and tFillADVs )
            {
                uint tNumADVs = aADVs.length();
                tVariableIndices.resize( tNumADVs, 1 );
                tADVIndices.resize( tNumADVs, 1 );
                for ( uint tIndex = 0; tIndex < tNumADVs; tIndex++ )
                {
                    tVariableIndices( tIndex ) = tIndex;
                    tADVIndices( tIndex )              = tIndex;
                }
            }
            else if ( tFillVariables )
            {
                tVariableIndices.resize( tADVIndices.length(), 1 );
                for ( uint tIndex = 0; tIndex < tADVIndices.length(); tIndex++ )
                {
                    tVariableIndices( tIndex ) = tIndex;
                }
            }
            else if ( tFillADVs )
            {
                tADVIndices.resize( tVariableIndices.length(), 1 );
                for ( uint tIndex = 0; tIndex < tVariableIndices.length(); tIndex++ )
                {
                    tADVIndices( tIndex ) = tIndex;
                }
            }

            // Constant parameters
            tConstants = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "constant_parameters" ) );
        }

        // Name of the field
        std::string tName = aFieldParameterList.get< std::string >( "name" );

        // Build Geometry
        if ( tFieldType == "constant" )
        {
            return std::make_shared< Constant_Property >(
                    aADVs,
                    tVariableIndices,
                    tADVIndices,
                    tConstants,
                    tName );
        }
        else if ( tFieldType == "scaled_field" )
        {
            return std::make_shared< Scaled_Field >(
                    aFieldDependencies( 0 ),
                    aADVs,
                    tVariableIndices,
                    tADVIndices,
                    tConstants,
                    tName );
        }
        else if ( tFieldType == "circle" )
        {
            return std::make_shared< Circle >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "superellipse" )
        {
            return std::make_shared< Superellipse >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "sphere" )
        {
            return std::make_shared< Sphere >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "superellipsoid" )
        {
            return std::make_shared< Superellipsoid >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "plane" )
        {
            return std::make_shared< Plane >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "nodal_field" )
        {
            MORIS_ERROR( aMTKMesh != nullptr, "Mesh is a null ptr for nodal field geometry" );

            return std::make_shared< ge::Mesh_Field_Geometry >( aMTKMesh, tName, mtk::EntityRank::NODE );
        }
        else if ( tFieldType == "nodal_field_from_file" )
        {
            std::string tFileName    = aFieldParameterList.get< std::string >( "file_name" );
            std::string tFieldName   = aFieldParameterList.get< std::string >( "field_name" );
            std::string tFieldFormat = aFieldParameterList.get< std::string >( "file_format" );
            real        tOffset      = aFieldParameterList.get< real >( "offset" );

            return std::make_shared< ge::Mesh_Field_Geometry >(
                    aMTKMesh,
                    tFileName,
                    tFieldName,
                    tFieldFormat,
                    tOffset,
                    mtk::EntityRank::NODE );
        }
        else if ( tFieldType == "sdf_field" )
        {
            std::string      tObjectPath   = aFieldParameterList.get< std::string >( "sdf_object_path" );
            Matrix< DDRMat > tObjectOffset = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "sdf_object_offset" ) );
            real             tSDFShift     = aFieldParameterList.get< real >( "sdf_shift" );

            return std::make_shared< ge::Geometry_SDF >(
                    tObjectPath,
                    tObjectOffset,
                    tSDFShift );
        }
        else if ( tFieldType == "image_sdf" )
        {
            // Get voxel-specific info
            std::string      tImageFileName    = aFieldParameterList.get< std::string >( "image_file" );
            Matrix< DDRMat > tDomainDimensions = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "image_dimensions" ) );
            Matrix< DDRMat > tDomainOffset     = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "image_offset" ) );

            real tSDFScaling = aFieldParameterList.get< real >( "image_sdf_scaling" );
            real tSDFShift   = aFieldParameterList.get< real >( "image_sdf_shift" );
            real tSDFDefault = aFieldParameterList.get< real >( "image_sdf_default" );
            bool tsDFInterp  = aFieldParameterList.get< bool >( "image_sdf_interpolate" );

            return std::make_shared< Image_SDF_Geometry >(
                    tImageFileName,
                    tDomainDimensions,
                    tDomainOffset,
                    tSDFScaling,
                    tSDFShift,
                    tSDFDefault,
                    tsDFInterp );
        }
        else if ( tFieldType == "user_defined" )
        {
            // Check if library is given
            MORIS_ERROR( aLibrary != nullptr, "Library must be given in order to create a user-defined geometry." );

            // Get sensitivity function if needed
            std::string tSensitivityFunctionName = aFieldParameterList.get< std::string >( "sensitivity_function_name" );
            Sensitivity_Function tSensitivityFunction =
                    ( tSensitivityFunctionName.empty() ? nullptr : aLibrary->load_function< Sensitivity_Function >( tSensitivityFunctionName ) );

            // Create user-defined geometry
            return std::make_shared< User_Defined_Field >(
                    aLibrary->load_function< Field_Function >( aFieldParameterList.get< std::string >( "field_function_name" ) ),
                    tSensitivityFunction,
                    aADVs,
                    tVariableIndices,
                    tADVIndices,
                    tConstants,
                    tName );
        }
        else if ( tFieldType == "voxel" and aFieldDependencies( 0 ) )
        {
            return std::make_shared< Single_Grain >(
                    aFieldDependencies( 0 ),
                    aIndex );
        }
        else if ( tFieldType == "voxel" )
        {
            // Get voxel-specific info
            std::string      tVoxelFieldName    = aFieldParameterList.get< std::string >( "voxel_field_file" );
            Matrix< DDRMat > tDomainDimensions  = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "domain_dimensions" ) );
            Matrix< DDRMat > tDomainOffset      = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "domain_offset" ) );
            Matrix< DDRMat > tGrainIdToValueMap = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "grain_id_value_map" ) );

            return std::make_shared< Voxel_Input >(
                    tVoxelFieldName,
                    tDomainDimensions,
                    tDomainOffset,
                    tGrainIdToValueMap );
        }
        else if ( tFieldType == "swiss_cheese_slice" )
        {
            // Check for definition
            uint tNumXHoles      = (uint)aFieldParameterList.get< sint >( "number_of_x_holes" );
            uint tNumYHoles      = (uint)aFieldParameterList.get< sint >( "number_of_y_holes" );
            real tTargetXSpacing = aFieldParameterList.get< real >( "target_x_spacing" );
            real tTargetYSpacing = aFieldParameterList.get< real >( "target_y_spacing" );

            MORIS_ERROR( ( tNumXHoles > 1 && tNumYHoles > 1 ) || ( tTargetXSpacing && tTargetYSpacing ),
                    "In a swiss cheese parameter list, you must specify either a number of holes > 1 %s",
                    "or a target spacing in each direction.\n" );

            // FIXME re-enable swiss cheese from geometry side
//            if ( tNumXHoles )
//            {
//                return std::make_shared< Swiss_Cheese_Slice >(
//                        aFieldParameterList.get< real >( "left_bound" ),
//                        aFieldParameterList.get< real >( "right_bound" ),
//                        aFieldParameterList.get< real >( "bottom_bound" ),
//                        aFieldParameterList.get< real >( "top_bound" ),
//                        tNumXHoles,
//                        tNumYHoles,
//                        aFieldParameterList.get< real >( "hole_x_semidiameter" ),
//                        aFieldParameterList.get< real >( "hole_y_semidiameter" ),
//                        aFieldParameterList.get< real >( "superellipse_exponent" ),
//                        aFieldParameterList.get< real >( "superellipse_scaling" ),
//                        aFieldParameterList.get< real >( "superellipse_regularization" ),
//                        aFieldParameterList.get< real >( "superellipse_shift" ),
//                        aFieldParameterList.get< real >( "row_offset" ) );
//            }
//            else
//            {
//                return std::make_shared< Swiss_Cheese_Slice >(
//                        aFieldParameterList.get< real >( "left_bound" ),
//                        aFieldParameterList.get< real >( "right_bound" ),
//                        aFieldParameterList.get< real >( "bottom_bound" ),
//                        aFieldParameterList.get< real >( "top_bound" ),
//                        tTargetXSpacing,
//                        tTargetYSpacing,
//                        aFieldParameterList.get< real >( "hole_x_semidiameter" ),
//                        aFieldParameterList.get< real >( "hole_y_semidiameter" ),
//                        aFieldParameterList.get< real >( "superellipse_exponent" ),
//                        aFieldParameterList.get< real >( "superellipse_scaling" ),
//                        aFieldParameterList.get< real >( "superellipse_regularization" ),
//                        aFieldParameterList.get< real >( "superellipse_shift" ),
//                        aFieldParameterList.get< real >( "row_offset" ),
//                        aFieldParameterList.get< bool >( "allow_less_than_target_spacing" ) );
//            }
            return nullptr;
        }
        else
        {
            MORIS_ERROR( false, "%s is not recognized as a valid Geometry type in fn_GEN_create_geometry.", tFieldType.c_str() );
            return nullptr;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
