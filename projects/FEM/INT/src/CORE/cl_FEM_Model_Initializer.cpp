//
// Created by frank on 12/14/23.
//

#include "cl_FEM_Model_Initializer.hpp"
#include <memory>

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

    void Model_Initializer::create_properties()
    {
        Vector< ParameterList > tPropParameterList = mParameterList( 0 );

        // loop over the parameter lists
        for ( auto const &tPropParameter : tPropParameterList )
        {
            // get property parameter list
            std::shared_ptr< Property > tProperty = std::make_shared< fem::Property >();

            // get property name from parameter list
            std::string tPropertyName = tPropParameter.get< std::string >( "property_name" );
            tProperty->set_name( tPropertyName );
            mProperties[ tPropertyName ] = tProperty;

            // set dof dependencies
            auto tDofTypes = this->parameter_to_vec_of_vec( tPropParameter, "dof_dependencies", mMSIDofTypeMap );
            tProperty->set_dof_type_list( tDofTypes );

            // set dv dependencies
            auto tDvTypes = parameter_to_vec_of_vec( tPropParameter, "dv_dependencies", mMSIDvTypeMap );
            tProperty->set_dv_type_list( tDvTypes );

            // set field dependencies
            auto tFieldTypes = parameter_to_vec_of_vec( tPropParameter, "field_dependencies", mMTKFieldTypeMap );
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
        }
    }

    void Model_Initializer::create_fields()
    {
        Vector< ParameterList > tFieldParameterList = mParameterList( FEM_PARAMETER_FIELD_INDEX );

        // loop over the parameter lists
        for ( auto const &tFieldParameter : tFieldParameterList )
        {
            moris::map< std::string, mtk::Field_Entity_Type > tFieldEntityTypeMap = mtk::get_field_entity_type_map();
            mtk::Field_Entity_Type                            tFieldEntityType    = tFieldEntityTypeMap.find( tFieldParameter.get< std::string >( "field_entity_type" ) );
            std::shared_ptr< fem::Field >                     tField              = std::make_shared< fem::Field >( *mMeshPair, tFieldEntityType );

            std::string tFieldName = tFieldParameter.get< std::string >( "field_name" );
            tField->set_label( tFieldName );
            mFields[ tFieldName ] = tField;

            // set field type
            moris::map< std::string, mtk::Field_Type > tFieldTypeMap =
                    mtk::get_field_type_map();

            // set field type
            Vector< mtk::Field_Type > tFieldTypes = string_to_cell< mtk::Field_Type >( tFieldParameter.get< std::string >( "field_type" ), tFieldTypeMap );
            tField->set_field_type( tFieldTypes );
            mFieldTypeToName[ static_cast< uint >( tFieldTypes( 0 ) ) ] = tFieldName;

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
        }
    }

    void Model_Initializer::print_physics_model()
    {
        ParameterList tComputationParameterList = this->mParameterList( 5 )( 0 );
        bool          tPrintPhysics             = tComputationParameterList.get< bool >( "print_physics_model" );
        if ( tPrintPhysics && par_rank() == 0 )
        {
            std::cout << "Set info \n";
            for ( auto &[ _, tSetInfo ] : mSetInfo )
            {
                std::cout << "%-------------------------------------------------\n";
                tSetInfo.print_names();
                std::cout << "%-------------------------------------------------\n";
            }
        }
    }


    Vector< fem::PropertyFunc >
    Model_Initializer::load_library_property_functions( ParameterList const &aParameterList, std::string const &aPropertyName )
    {
        auto                        tFuncNames    = string_to_cell< std::string >( aParameterList.get< std::string >( aPropertyName ) );
        uint                        tNumFunctions = tFuncNames.size();
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

    Vector< std::pair< std::shared_ptr< Property >, std::string > > Model_Initializer::read_properties( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const
    {
        std::string const                     tKey             = get_leader_follower_key( "properties", aLeaderFollower );    // returns "properties", "leader_properties", or "follower_properties"
        Vector< Vector< std::string > > const tPropertyStrings = string_to_cell_of_cell< std::string >( aParameterList.get< std::string >( tKey ) );

        Vector< std::pair< std::shared_ptr< Property >, std::string > > tPropertyNamesPairVec;
        tPropertyNamesPairVec.reserve( tPropertyStrings.size() );

        for ( uint iProp = 0; iProp < tPropertyStrings.size(); iProp++ )
        {
            std::string tPropertyName   = tPropertyStrings( iProp )( 0 );
            std::string tPropertyString = tPropertyStrings( iProp )( 1 );
            ensure_existing_parameter( mProperties, tPropertyName, tKey );
            tPropertyNamesPairVec.emplace_back( mProperties.at( tPropertyName ), tPropertyString );
        }
        return tPropertyNamesPairVec;
    }

    std::string Model_Initializer::get_leader_follower_key( std::string aKey, mtk::Leader_Follower const aLeaderFollower ) const
    {
        switch ( aLeaderFollower )
        {
            case mtk::Leader_Follower::LEADER:
                return "leader_" + aKey;
            case mtk::Leader_Follower::FOLLOWER:
                return "follower_" + aKey;
            default:
                return aKey;
        }
    }

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
