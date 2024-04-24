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
#include "GEN_Data_Types.hpp"
#include "fn_Parsing_Tools.hpp"

#include "cl_GEN_Constant_Field.hpp"
#include "cl_GEN_Scaled_Field.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Superellipse.hpp"
#include "cl_GEN_Line.hpp"
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
            Parameter_List                     aFieldParameterList,
            Vector< real >&                    aADVs,
            Vector< std::shared_ptr< Field > > aFieldDependencies,
            std::shared_ptr< Library_IO >      aLibrary,
            mtk::Mesh*                         aMTKMesh )
    {
        // Field type
        auto tFieldType = aFieldParameterList.get< gen::Field_Type >( "field_type" );

        // ADV inputs
        Vector< uint > tVariableIndices = aFieldParameterList.get< Vector< uint > >( "field_variable_indices" );
        Vector< uint > tADVIndices = aFieldParameterList.get< Vector< uint > >( "adv_indices" );
        Vector< real > tConstants = aFieldParameterList.get< Vector< real > >( "constant_parameters" );

        // Name of the field
        std::string tName = aFieldParameterList.get< std::string >( "name" );

        // Shared pointer to field
        std::shared_ptr< Field > tField;

        // Build field
        switch ( tFieldType )
        {
            case gen::Field_Type::NONE:
                // TODO Right now this isn't possible. In the future, make this an assignable field
                break;
            case gen::Field_Type::CONSTANT:
            {
                tField = std::make_shared< Constant_Field >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::LINE:
            {
                tField = std::make_shared< Line >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::CIRCLE:
            {
                tField = std::make_shared< Circle >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::SUPERELLIPSE:
            {
                tField = std::make_shared< Superellipse >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::PLANE:
            {
                tField = std::make_shared< Plane >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::SPHERE:
            {
                tField = std::make_shared< Sphere >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::SUPERELLIPSOID:
            {
                tField = std::make_shared< Superellipsoid >( aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::SCALED_FIELD:
            {
                tField = std::make_shared< Scaled_Field >( aFieldDependencies( 0 ), aADVs, tVariableIndices, tADVIndices, tConstants, tName );
                break;
            }
            case gen::Field_Type::COMBINED_FIELDS:
            {
                bool tUseMinimum = true;
                if ( tConstants.size() == 1 )
                {
                    tUseMinimum = tConstants( 0 );
                }
                tField = std::make_shared< Combined_Fields >( aFieldDependencies, tUseMinimum, tName );
                break;
            }
            case gen::Field_Type::NODAL:
            {
                MORIS_ERROR( aMTKMesh != nullptr, "Mesh is a null ptr for nodal field." );
                tField = std::make_shared< gen::Mesh_Field >( aMTKMesh, tName, mtk::EntityRank::NODE );
                break;
            }
            case gen::Field_Type::NODAL_FROM_FILE:
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
                break;
            }
            case gen::Field_Type::SIGNED_DISTANCE_OBJECT:
            {
                std::string  tObjectPath   = aFieldParameterList.get< std::string >( "sdf_object_path" );
                Vector< real > tObjectOffset = string_to_cell< real >( aFieldParameterList.get< std::string >( "sdf_object_offset" ) );
                real         tSDFShift     = aFieldParameterList.get< real >( "sdf_shift" );

                tField = std::make_shared< gen::Signed_Distance_Field >(
                        tObjectPath,
                        tObjectOffset,
                        tSDFShift );
                break;
            }
            case gen::Field_Type::SIGNED_DISTANCE_IMAGE:
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
                break;
            }
            case gen::Field_Type::USER_DEFINED:
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
