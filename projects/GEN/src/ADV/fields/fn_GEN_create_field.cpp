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
            const Parameter_List&                aFieldParameterList,
            ADV_Manager&                         aADVManager,
            Vector< std::shared_ptr< Field > >   aFieldDependencies,
            const std::shared_ptr< Library_IO >& aLibrary,
            mtk::Mesh*                           aMTKMesh )
    {
        // Field type
        auto tFieldType = aFieldParameterList.get< Field_Type >( "field_type" );

        // Name of the field
        std::string tName = aFieldParameterList.get< std::string >( "name" );

        // Shared pointer to field
        std::shared_ptr< Field > tField;

        // Build field
        switch ( tFieldType )
        {
            case Field_Type::NONE:
            {
                MORIS_ERROR( false, "A field must be created with a field type." );
                break;
            }
            case Field_Type::CONSTANT:
            {
                ADV tConstant = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "constant" ) );
                tField = std::make_shared< Constant_Field >( tConstant, tName );
                break;
            }
            case Field_Type::LINE:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tNormalX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "normal_x" ) );
                ADV tNormalY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "normal_y" ) );
                tField = std::make_shared< Line >( tCenterX, tCenterY, tNormalX, tNormalY, tName );
                break;
            }
            case Field_Type::CIRCLE:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tRadius = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "radius" ) );
                tField = std::make_shared< Circle >( tCenterX, tCenterY, tRadius, tName );
                break;
            }
            case Field_Type::SUPERELLIPSE:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tSemidiameterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "semidiameter_x" ) );
                ADV tSemidiameterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "semidiameter_y" ) );
                real tExponent = aFieldParameterList.get< real >( "exponent" );
                tField = std::make_shared< Superellipse >( tCenterX, tCenterY, tSemidiameterX, tSemidiameterY, tExponent, 1.0, 0.0 );
                break;
            }
            case Field_Type::PLANE:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tCenterZ = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_z" ) );
                ADV tNormalX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "normal_x" ) );
                ADV tNormalY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "normal_y" ) );
                ADV tNormalZ = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "normal_z" ) );
                tField = std::make_shared< Plane >( tCenterX, tCenterY, tCenterZ, tNormalX, tNormalY, tNormalZ, tName );
                break;
            }
            case Field_Type::SPHERE:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tCenterZ = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_z" ) );
                ADV tRadius = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "radius" ) );
                tField = std::make_shared< Sphere >( tCenterX, tCenterY, tCenterZ, tRadius, tName );
                break;
            }
            case Field_Type::SUPERELLIPSOID:
            {
                ADV tCenterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_x" ) );
                ADV tCenterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_y" ) );
                ADV tCenterZ = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "center_z" ) );
                ADV tSemidiameterX = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "semidiameter_x" ) );
                ADV tSemidiameterY = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "semidiameter_y" ) );
                ADV tSemidiameterZ = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "semidiameter_z" ) );
                real tExponent = aFieldParameterList.get< real >( "exponent" );
                tField = std::make_shared< Superellipsoid >( tCenterX, tCenterY, tCenterZ, tSemidiameterX, tSemidiameterY, tSemidiameterZ, tExponent, tName );
                break;
            }
            case Field_Type::SCALED_FIELD:
            {
                ADV tScalingFactor = aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( "scaling_factor" ) );
                tField = std::make_shared< Scaled_Field >( aFieldDependencies( 0 ), tScalingFactor, tName );
                break;
            }
            case Field_Type::COMBINED_FIELDS:
            {
                bool tUseMinimum = aFieldParameterList.get< bool >( "use_minimum" );
                tField = std::make_shared< Combined_Fields >( aFieldDependencies, tUseMinimum, tName );
                break;
            }
            case Field_Type::NODAL_FROM_FILE:
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
            case Field_Type::SIGNED_DISTANCE_OBJECT:
            {
                auto tObjectPath   = aFieldParameterList.get< std::string >( "sdf_object_path" );
                auto tObjectOffset = aFieldParameterList.get< Vector< real > >( "sdf_object_offset" );
                auto tSDFShift     = aFieldParameterList.get< real >( "sdf_shift" );

                tField = std::make_shared< gen::Signed_Distance_Field >(
                        tObjectPath,
                        tObjectOffset,
                        tSDFShift );
                break;
            }
            case Field_Type::SIGNED_DISTANCE_IMAGE:
            {
                // Get SDF-specific info
                auto tImageFileName    = aFieldParameterList.get< std::string >( "image_file" );
                auto tDomainDimensions = aFieldParameterList.get< Vector< real > >( "image_dimensions" );
                auto tDomainOffset     = aFieldParameterList.get< Vector< real > >( "image_offset" );
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
            case Field_Type::USER_DEFINED:
            {
                // Check if library is given
                MORIS_ERROR( aLibrary != nullptr, "Library must be given in order to create a user-defined field." );

                // get field function
                auto tFieldFunction = aLibrary->load_function< Field_Function >( aFieldParameterList.get< std::string >( "field_function_name" ) );

                // Get sensitivity function if needed
                Sensitivity_Function tSensitivityFunction = nullptr;
                std::string tSensitivityFunctionName = aFieldParameterList.get< std::string >( "sensitivity_function_name" );
                if ( not tSensitivityFunctionName.empty() )
                {
                    tSensitivityFunction = aLibrary->load_function< Sensitivity_Function >( tSensitivityFunctionName );
                }

                // Loop over parameters to create ADVs
                Vector< ADV > tADVs;
                for ( const auto& iParameter : aFieldParameterList )
                {
                    // Determine if parameter is design variable
                    if ( iParameter.get_parameter().index() == variant_index< Design_Variable >() )
                    {
                        // Get design variable from parameter list
                        tADVs.push_back( aADVManager.create_adv( aFieldParameterList.get< Design_Variable >( iParameter.get_name() ) ) );
                    }
                }

                // Create user-defined field
                tField = std::make_shared< User_Defined_Field >(
                        tFieldFunction,
                        tSensitivityFunction,
                        tADVs,
                        tName );
                break;
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

    Vector< char > get_active_parameter_ids( const Parameter_List& aFieldParameterList )
    {
        // Vector of parameter IDs
        Vector< char > tParameterIDs;

        // Loop over all parameters
        for ( const auto& iParameter : aFieldParameterList )
        {
            // Determine if parameter is design variable
            if ( iParameter.get_parameter().index() == variant_index< Design_Variable >() )
            {
                // Get design variable from parameter list
                auto tDesignVariable = aFieldParameterList.get< Design_Variable >( iParameter.get_name() );

                // Add ID if not constant
                if ( not tDesignVariable.is_constant() )
                {
                    tParameterIDs.push_back( tDesignVariable.get_id() );
                }
            }
        }

        // Return resulting vector
        return tParameterIDs;
    }

    //--------------------------------------------------------------------------------------------------------------
}
