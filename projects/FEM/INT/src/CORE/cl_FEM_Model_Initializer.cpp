/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Model_Initializer.cpp
 *
 */

#include "cl_FEM_Model_Initializer.hpp"

namespace moris::fem
{
    void Model_Initializer::initialize()
    {
        this->create_properties();
        this->create_fields();
        this->create_material_models();
        this->create_constitutive_models();
        this->create_stabilization_parameters();
        this->create_iwgs();
        this->create_iqis();
        this->create_set_info();
        this->print_physics_model();
    }

    //----------------------------------------------------------------

    void Model_Initializer::create_properties()
    {
        Vector< Parameter_List > tPropParameterList = mParameterList( 0 );
        uint                     tNumProps          = tPropParameterList.size();

        mProperties.resize( tNumProps, nullptr );

        // loop over the parameter lists
        for ( uint iProp = 0; iProp < tNumProps; iProp++ )
        {
            // get property parameter list
            Parameter_List tPropParameter = tPropParameterList( iProp );
            auto           tProperty      = std::make_shared< fem::Property >();

            // get property name from parameter list
            std::string tPropertyName = tPropParameter.get< std::string >( "property_name" );
            tProperty->set_name( tPropertyName );

            // fill property map
            mPropertyMap[ tPropertyName ] = iProp;

            // set dof dependencies
            auto tDofTypes = this->property_to_vec_of_vec( tPropParameter, "dof_dependencies", mMSIDofTypeMap );
            tProperty->set_dof_type_list( tDofTypes );

            // set dv dependencies
            auto tDvTypes = property_to_vec_of_vec( tPropParameter, "dv_dependencies", mMSIDvTypeMap );
            tProperty->set_dv_type_list( tDvTypes );

            // set field dependencies
            auto tFieldTypes = property_to_vec_of_vec( tPropParameter, "field_dependencies", mFieldTypeMap );
            tProperty->set_field_type_list( tFieldTypes );

            // set function parameters
            auto tFuncParameters = string_to_cell_mat_2< DDRMat >( tPropParameter.get< std::string >( "function_parameters" ) );
            tProperty->set_parameters( tFuncParameters );

            // set value function for property
            std::string  tValFuncName = tPropParameter.get< std::string >( "value_function" );
            FEM_Function tValFunction = nullptr;
            if ( tValFuncName.size() > 1 )
            {
                tValFunction = mLibrary->load_function< FEM_Function >( tValFuncName );
                tProperty->set_val_function( tValFunction );
            }

            // set dof derivative function for property
            Vector< fem::PropertyFunc > tDofDerFunctions = load_library_property_functions( tPropParameter, "dof_derivative_functions" );
            tProperty->set_dof_derivative_functions( tDofDerFunctions );

            // set dv derivative function for property
            Vector< fem::PropertyFunc > tDvDerFunctions = load_library_property_functions( tPropParameter, "dv_derivative_functions" );
            tProperty->set_dv_derivative_functions( tDvDerFunctions );

            // set space derivative function for property
            Vector< fem::PropertyFunc > tSpaceDerFunctions = load_library_property_functions( tPropParameter, "space_derivative_functions" );
            tProperty->set_space_der_functions( tSpaceDerFunctions );

            mProperties( iProp ) = tProperty;
        }
    }

    //----------------------------------------------------------------

    void Model_Initializer::create_fields()
    {
        Vector< Parameter_List > tFieldParameterList = mParameterList( 6 );
        sint                     tNumFields          = tFieldParameterList.size();

        mFields.resize( tNumFields, nullptr );

        // loop over the parameter lists
        for ( sint iField = 0; iField < tNumFields; iField++ )
        {
            // get property parameter list
            Parameter_List tFieldParameter = tFieldParameterList( iField );

            moris::map< std::string, mtk::Field_Entity_Type > tFieldEntityTypeMap = mtk::get_field_entity_type_map();

            mtk::Field_Entity_Type        tFieldEntityType = tFieldEntityTypeMap.find( tFieldParameter.get< std::string >( "field_entity_type" ) );
            std::shared_ptr< fem::Field > tField           = std::make_shared< fem::Field >( *mMeshPair, tFieldEntityType );

            std::string tFieldName = tFieldParameter.get< std::string >( "field_name" );
            tField->set_label( tFieldName );

            // fill property map
            mFieldMap[ tFieldName ] = iField;

            // set field type
            moris::map< std::string, mtk::Field_Type > tFieldTypeMap = mtk::get_field_type_map();

            // set field type
            Vector< mtk::Field_Type > tFieldTypes = string_to_cell< mtk::Field_Type >( tFieldParameter.get< std::string >( "field_type" ), tFieldTypeMap );
            tField->set_field_type( tFieldTypes );

            mFieldTypes.resize( std::max( static_cast< uint >( tFieldTypes( 0 ) ) + 1, (uint)mFieldTypes.size() ), -1 );
            mFieldTypes( static_cast< uint >( tFieldTypes( 0 ) ) ) = iField;

            MORIS_ERROR( ( tFieldParameter.get< std::string >( "field_create_from_file" ).empty() ) or ( tFieldParameter.get< std::string >( "IQI_Name" ).empty() ),
                    "FEM_Model::create_fields(); Field must be either created based on IQI or read from file." );

            MORIS_ERROR( not( ( not tFieldParameter.get< std::string >( "field_create_from_file" ).empty() ) and ( not tFieldParameter.get< std::string >( "IQI_Name" ).empty() ) ),
                    "FEM_Model::create_fields(); Field must be either created based on IQI or read from file." );

            if ( not tFieldParameter.get< std::string >( "field_create_from_file" ).empty() )
            {
                tField->set_field_from_file(
                        tFieldParameter.get< std::string >( "field_create_from_file" ),
                        tFieldParameter.get< sint >( "field_file_time_index" ),
                        tFieldParameter.get< sint >( "field_file_field_index" ) );
            }

            if ( not tFieldParameter.get< std::string >( "IQI_Name" ).empty() )
            {
                tField->set_IQI_name( tFieldParameter.get< std::string >( "IQI_Name" ) );
            }

            if ( not tFieldParameter.get< std::string >( "field_output_to_file" ).empty() )
            {
                tField->set_field_to_file( tFieldParameter.get< std::string >( "field_output_to_file" ) );
            }

            mFields( iField ) = tField;
        }
    }

    //----------------------------------------------------------------

    void Model_Initializer::print_physics_model()
    {
        Parameter_List tComputationParameterList = this->mParameterList( 5 )( 0 );
        bool           tPrintPhysics             = tComputationParameterList.get< bool >( "print_physics_model" );
        if ( tPrintPhysics && par_rank() == 0 )
        {
            std::cout << "Set info \n";
            for ( auto &tSetInfo : mSetInfo )
            {
                std::cout << "%-------------------------------------------------\n";
                tSetInfo.print_names();
                std::cout << "%-------------------------------------------------\n";
            }
        }
    }

    //----------------------------------------------------------------

    Vector< fem::PropertyFunc >
    Model_Initializer::load_library_property_functions(
            Parameter_List const &aParameterList,
            std::string const    &aPropertyName )
    {
        auto tFuncNames    = string_to_cell< std::string >( aParameterList.get< std::string >( aPropertyName ) );
        uint tNumFunctions = tFuncNames.size();

        Vector< fem::PropertyFunc > tPropertyFunctions( tNumFunctions, nullptr );

        for ( uint iFunc = 0; iFunc < tNumFunctions; iFunc++ )
        {
            if ( tFuncNames( iFunc ).size() > 1 )
            {
                tPropertyFunctions( iFunc ) = mLibrary->load_function< FEM_Function >( tFuncNames( iFunc ) );
            }
        }
        return tPropertyFunctions;
    }

    //----------------------------------------------------------------

    void Model_Initializer::check_and_set_ghost_set_names( std::string &aMeshSetName, MSI::Dof_Type aDofType )
    {
        // check whether new ghost sets should be used
        if ( mUseNewGhostSets )
        {
            if ( aMeshSetName.find( "ghost_p" ) != std::string::npos )
            {
                // find the phase index
                size_t tPos = aMeshSetName.find( 'p' );
                MORIS_ERROR(
                        tPos != std::string::npos,
                        "Model_Initializer::check_and_set_ghost_set_names() - Phase index not found in ghost set name." );
                aMeshSetName.erase( 0, tPos + 1 );

                moris_index tBsplineMeshIndex = mDofTypeToBsplineMeshIndex.find( aDofType )->second;

                aMeshSetName = "ghost_B" + std::to_string( tBsplineMeshIndex ) + "_p" + aMeshSetName;
            }
        }
    }
}    // namespace moris::fem
