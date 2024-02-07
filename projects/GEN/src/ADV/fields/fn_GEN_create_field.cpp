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

#include "cl_GEN_Constant_Field.hpp"
#include "cl_GEN_Scaled_Field.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Superellipse.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Superellipsoid.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Combined_Fields.hpp"
#include "cl_GEN_Field_Array_Factory.hpp"
#include "cl_GEN_Mesh_Field.hpp"
#include "cl_GEN_Signed_Distance_Field.hpp"
#include "cl_GEN_Image_Signed_Distance_Field.hpp"

namespace moris::gen
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
        // Field type
        std::string tFieldType = aFieldParameterList.get< std::string >( "field_type" );

        // ADV inputs
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
                    tADVIndices( tIndex )      = tIndex;
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

        // Shared pointer to field
        std::shared_ptr< Field > tField;

        // Build field
        if ( tFieldType == "constant" )
        {
            tField = std::make_shared< Constant_Field >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "scaled_field" )
        {
            tField = std::make_shared< Scaled_Field >( aFieldDependencies( 0 ), aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "circle" )
        {
            tField = std::make_shared< Circle >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "superellipse" )
        {
            tField = std::make_shared< Superellipse >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "sphere" )
        {
            tField = std::make_shared< Sphere >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "superellipsoid" )
        {
            tField = std::make_shared< Superellipsoid >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "plane" )
        {
            tField = std::make_shared< Plane >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
        }
        else if ( tFieldType == "combined_fields" )
        {
            bool tUseMinimum = true;
            if ( tConstants.length() == 1 )
            {
                tUseMinimum = tConstants( 0 );
            }
            tField = std::make_shared< Combined_Fields >( aFieldDependencies, tUseMinimum, tName );
        }
        else if ( tFieldType == "nodal_field" )
        {
            MORIS_ERROR( aMTKMesh != nullptr, "Mesh is a null ptr for nodal field." );

            tField = std::make_shared< gen::Mesh_Field >( aMTKMesh, tName, mtk::EntityRank::NODE );
        }
        else if ( tFieldType == "nodal_field_from_file" )
        {
            std::string tFileName    = aFieldParameterList.get< std::string >( "file_name" );
            std::string tFieldName   = aFieldParameterList.get< std::string >( "field_name" );
            std::string tFieldFormat = aFieldParameterList.get< std::string >( "file_format" );
            real        tOffset      = aFieldParameterList.get< real >( "offset" );

            tField = std::make_shared< gen::Mesh_Field >(
                    aMTKMesh,
                    tFileName,
                    tFieldName,
                    tFieldFormat,
                    tOffset,
                    mtk::EntityRank::NODE );
        }
        else if ( tFieldType == "sdf_field" )
        {
            std::string  tObjectPath   = aFieldParameterList.get< std::string >( "sdf_object_path" );
            Cell< real > tObjectOffset = string_to_cell< real >( aFieldParameterList.get< std::string >( "sdf_object_offset" ) );
            real         tSDFShift     = aFieldParameterList.get< real >( "sdf_shift" );

            tField = std::make_shared< gen::Signed_Distance_Field >(
                    tObjectPath,
                    tObjectOffset,
                    tSDFShift );
        }
        else if ( tFieldType == "image_sdf" )
        {
            // Get SDF-specific info
            std::string      tImageFileName    = aFieldParameterList.get< std::string >( "image_file" );
            Matrix< DDRMat > tDomainDimensions = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "image_dimensions" ) );
            Matrix< DDRMat > tDomainOffset     = string_to_mat< DDRMat >( aFieldParameterList.get< std::string >( "image_offset" ) );

            real tSDFScaling = aFieldParameterList.get< real >( "image_sdf_scaling" );
            real tSDFShift   = aFieldParameterList.get< real >( "image_sdf_shift" );
            real tSDFDefault = aFieldParameterList.get< real >( "image_sdf_default" );
            bool tsDFInterp  = aFieldParameterList.get< bool >( "image_sdf_interpolate" );

            tField = std::make_shared< Image_Signed_Distance_Field >(
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
            MORIS_ERROR( aLibrary != nullptr, "Library must be given in order to create a user-defined field." );

            // Get sensitivity function if needed
            std::string          tSensitivityFunctionName = aFieldParameterList.get< std::string >( "sensitivity_function_name" );
            Sensitivity_Function tSensitivityFunction =
                    ( tSensitivityFunctionName.empty() ? nullptr : aLibrary->load_function< Sensitivity_Function >( tSensitivityFunctionName ) );

            // Create user-defined field
            tField = std::make_shared< User_Defined_Field >(
                    aLibrary->load_function< Field_Function >( aFieldParameterList.get< std::string >( "field_function_name" ) ),
                    tSensitivityFunction,
                    aADVs,
                    tVariableIndices,
                    tADVIndices,
                    tConstants,
                    tName );
        }
        else
        {
            MORIS_ERROR( false, "%s is not recognized as a valid field type in fn_GEN_create_field().", tFieldType.c_str() );
            tField = nullptr;
        }

        // Check for definition of array parameters
        if ( aFieldParameterList.exists( "number_of_fields_x" ) )
        {
            // Return field array
            Field_Array_Factory tFieldArrayFactory( aFieldParameterList );
            return tFieldArrayFactory.create_field_array(
                    tField,
                    aFieldParameterList.get< bool >( "minimum" ) );
        }
        else
        {
            // Return single field
            return tField;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
