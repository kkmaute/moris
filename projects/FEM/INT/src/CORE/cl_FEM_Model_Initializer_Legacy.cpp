//
// Created by frank on 12/14/23.
//

#include "cl_FEM_Model_Initializer_Legacy.hpp"
#include "cl_Cell.hpp"
#include "cl_Param_List.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"

namespace moris::fem
{
    void Model_Initializer_Legacy::create_material_models()
    {
        std::map< std::string, uint > tMMMap;

        // create a material model factory
        MM_Factory tMMFactory;

        // get the MM parameter list
        moris::Cell< ParameterList > tMMParameterList = mParameterList( 7 );

        // get number of material models
        uint tNumMMs = tMMParameterList.size();

        // create a list of MMs
        mMaterialModels.resize( tNumMMs, nullptr );

        // loop over the parameter lists for MM
        for ( uint iMM = 0; iMM < tNumMMs; iMM++ )
        {
            // get the material type from parameter list
            fem::Material_Type tMMType =
                    static_cast< fem::Material_Type >( tMMParameterList( iMM ).get< uint >( "material_type" ) );

            // create a material model pointer
            mMaterialModels( iMM ) = tMMFactory.create_MM( tMMType );

            // set MM name
            mMaterialModels( iMM )->set_name( tMMParameterList( iMM ).get< std::string >( "material_name" ) );

            // fill MM map
            tMMMap[ tMMParameterList( iMM ).get< std::string >( "material_name" ) ] = iMM;

            // set MM space dimension
            mMaterialModels( iMM )->set_space_dim( mSpatialDimension );

            // set MM dof dependencies
            moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes = string_to_cell_of_cell< moris::MSI::Dof_Type >(
                    std::get< 0 >( tMMParameterList( iMM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),

                    mMSIDofTypeMap );
            moris::Cell< std::string > tDofTypeNames;
            string_to_cell(
                    std::get< 1 >( tMMParameterList( iMM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypeNames );
            mMaterialModels( iMM )->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set MM properties
            moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
            string_to_cell_of_cell(
                    tMMParameterList( iMM ).get< std::string >( "properties" ),
                    tPropertyNamesPair );
            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                // if property name is in the property map
                if ( mPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != mPropertyMap.end() )
                {
                    // get property index
                    uint tPropertyIndex = mPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                    // set property for MM
                    mMaterialModels( iMM )->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ) );
                }
                else
                {
                    // error message for unknown property
                    MORIS_ERROR( false,
                            "Model_Initializer::create_MMs - Unknown aPropertyString : %s \n",
                            tPropertyNamesPair( iProp )( 0 ).c_str() );
                }
            }
            // set local properties
            mMaterialModels( iMM )->set_local_properties();
        }
    }

    void Model_Initializer_Legacy::create_constitutive_models()
    {
        // create a constitutive model factory
        CM_Factory tCMFactory;

        // get the CM parameter list
        moris::Cell< ParameterList > tCMParameterList = mParameterList( 1 );

        // get number of constitutive models
        uint tNumCMs = tCMParameterList.size();

        // create a list of CMs
        mConstitutiveModels.resize( tNumCMs, nullptr );

        // loop over the parameter lists for CM
        for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
        {
            // get the constitutive type from parameter list
            fem::Constitutive_Type tCMType =
                    static_cast< fem::Constitutive_Type >( tCMParameterList( iCM ).get< uint >( "constitutive_type" ) );

            // create a constitutive model pointer
            mConstitutiveModels( iCM ) = tCMFactory.create_CM( tCMType );

            // set CM name
            mConstitutiveModels( iCM )->set_name( tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) );

            // fill CM map
            mConstitutiveModelMap[ tCMParameterList( iCM ).get< std::string >( "constitutive_name" ) ] = iCM;

            // set CM model type
            fem::Model_Type tCMModelType =
                    static_cast< fem::Model_Type >( tCMParameterList( iCM ).get< uint >( "model_type" ) );
            if ( tCMModelType != fem::Model_Type::UNDEFINED )
            {
                mConstitutiveModels( iCM )->set_model_type( tCMModelType );
            }

            // set CM space dimension
            mConstitutiveModels( iCM )->set_space_dim( mSpatialDimension );

            // set CM dof dependencies
            moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
            string_to_cell_of_cell(
                    std::get< 0 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypes,
                    mMSIDofTypeMap );
            moris::Cell< std::string > tDofTypeNames;
            string_to_cell(
                    std::get< 1 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypeNames );
            mConstitutiveModels( iCM )->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set CM dv dependencies
            moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
            string_to_cell_of_cell(
                    std::get< 0 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                    tDvTypes,
                    mMSIDvTypeMap );
            moris::Cell< std::string > tDvTypeNames;
            string_to_cell(
                    std::get< 1 >( tCMParameterList( iCM ).get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                    tDvTypeNames );
            mConstitutiveModels( iCM )->set_dv_type_list( tDvTypes, tDvTypeNames );

            // set CM properties
            moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
            string_to_cell_of_cell(
                    tCMParameterList( iCM ).get< std::string >( "properties" ),
                    tPropertyNamesPair );
            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                // if property name is in the property map
                if ( mPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != mPropertyMap.end() )
                {
                    // get property index
                    uint tPropertyIndex = mPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                    // set property for CM
                    mConstitutiveModels( iCM )->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ) );
                }
                else
                {
                    // error message for unknown property
                    MORIS_ERROR( false,
                            "Model_Initializer::create_CMs - Unknown aPropertyString : %s \n",
                            tPropertyNamesPair( iProp )( 0 ).c_str() );
                }
            }
            // set local properties
            mConstitutiveModels( iCM )->set_local_properties();

            // set material model
            moris::Cell< moris::Cell< std::string > > tMMNamesPair;
            string_to_cell_of_cell(
                    tCMParameterList( iCM ).get< std::string >( "material_model" ),
                    tMMNamesPair );

            for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
            {
                // if MM name is in the CM map
                if ( mMaterialModelMap.find( tMMNamesPair( iMM )( 0 ) ) != mMaterialModelMap.end() )
                {
                    // get MM index
                    uint tMMIndex = mMaterialModelMap[ tMMNamesPair( iMM )( 0 ) ];

                    // set MM for IWG
                    mConstitutiveModels( iCM )->set_material_model(
                            mMaterialModels( tMMIndex ),
                            tMMNamesPair( iMM )( 1 ) );
                }
                else
                {
                    // error message unknown MM
                    MORIS_ERROR( false,
                            "Model_Initializer::create_CMs - Unknown aMMString: %s \n",
                            tMMNamesPair( iMM )( 0 ).c_str() );
                }
            }
        }
    }

    void Model_Initializer_Legacy::create_stabilization_parameters()
    {
        // create a stabilization parameter factory
        SP_Factory tSPFactory;

        // get the SP parameter list
        moris::Cell< ParameterList > tSPParameterList = mParameterList( 2 );

        // get the number of stabilization parameters
        uint tNumSPs = tSPParameterList.size();

        // set size for the list of stabilization parameter pointer
        mStabilizationParameters.resize( tNumSPs, nullptr );

        // loop over the parameter list
        for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
        {
            // get the SP parameter
            ParameterList tSPParameter = tSPParameterList( iSP );

            // get the stabilization type from parameter list
            fem::Stabilization_Type tSPType =
                    static_cast< fem::Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );

            // create a stabilization parameter pointer
            mStabilizationParameters( iSP ) = tSPFactory.create_SP( tSPType );

            // set name
            mStabilizationParameters( iSP )->set_name( tSPParameter.get< std::string >( "stabilization_name" ) );

            // set SP space dimension
            mStabilizationParameters( iSP )->set_space_dim( mSpatialDimension );

            // fill stabilization map
            mStabilizationParameterMap[ tSPParameter.get< std::string >( "stabilization_name" ) ] = iSP;

            // set parameters
            moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tSPParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );
            mStabilizationParameters( iSP )->set_parameters( tFuncParameters );

            // init string for leader or follower
            std::string          tIsLeaderString = "leader";
            mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

            // loop on leader and follower
            for ( uint iLeader = 0; iLeader <= mStabilizationParameters( iSP )->get_has_follower(); iLeader++ )
            {
                // if follower
                if ( iLeader )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                        tDofTypes,
                        mMSIDofTypeMap );
                moris::Cell< std::string > tDofTypeNames;
                string_to_cell( std::get< 1 >(
                                        tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                        tDofTypeNames );
                mStabilizationParameters( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsLeader );

                // set dv dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                        tDvTypes,
                        mMSIDvTypeMap );
                moris::Cell< std::string > tDvTypeNames;
                string_to_cell(
                        std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                        tDvTypeNames );
                mStabilizationParameters( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsLeader );

                // set leader properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tSPParameter.get< std::string >( tIsLeaderString + "_properties" ),
                        tPropertyNamesPair );

                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get the property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check for unknown property
                    MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                            "Model_Initializer::create_stabilization_parameters_without_phase - Unknown leader aPropertyString : %s \n",
                            tPropertyName.c_str() );

                    // get property index
                    uint tPropertyIndex = mPropertyMap[ tPropertyName ];

                    // set property for CM
                    mStabilizationParameters( iSP )->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ),
                            tIsLeader );
                }

                // set constitutive models
                moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tSPParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                // loop over the CM names
                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // get the CM name
                    std::string tCMName = tCMNamesPair( iCM )( 0 );

                    // check for unknown CM
                    MORIS_ERROR( mConstitutiveModelMap.find( tCMName ) != mConstitutiveModelMap.end(),
                            "Model_Initializer::create_stabilization_parameters_without_phase - Unknown leader aCMString: %s \n",
                            tCMName.c_str() );

                    // get CM index
                    uint tCMIndex = mConstitutiveModelMap[ tCMName ];

                    // set CM for SP
                    mStabilizationParameters( iSP )->set_constitutive_model(
                            mConstitutiveModels( tCMIndex ),
                            tCMNamesPair( iCM )( 1 ) );
                }
            }
        }
    }

    void Model_Initializer_Legacy::create_iwgs()
    {
        IWG_Factory                  tIWGFactory;
        moris::Cell< ParameterList > tIWGParameterList = mParameterList( 3 );
        uint const                   tNumIWGs          = tIWGParameterList.size();
        mIWGs.resize( tNumIWGs, nullptr );    // list of IWG pointers

        for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
        {
            // get the treated IWG parameter list
            ParameterList const tIWGParameter = tIWGParameterList( iIWG );

            fem::IWG_Type const tIWGType =
                    static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

            auto tIWG = tIWGFactory.create_IWG( tIWGType );

            std::string const tIWGName = tIWGParameter.get< std::string >( "IWG_name" );
            tIWG->set_name( tIWGName );

            // fill IWG map
            mIWGMap[ tIWGName ] = iIWG;

            // get function parameters
            auto tFuncParameters = string_to_cell_mat_2< DDRMat >( tIWGParameter.get< std::string >( "function_parameters" ) );
            tIWG->set_parameters( tFuncParameters );

            // get the ghost order from parameter list
            uint const tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );
            tIWG->set_interpolation_order( tGhostOrder );

            // set residual dof type
            set_iwg_residual_dof_type( tIWGParameter, tIWG );

            // loop over leader and follower and set the appropriate properties
            moris::Cell< mtk::Leader_Follower > const tLeaderFollower{ mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER };
            for ( auto const &tLeaderFollowerType : tLeaderFollower )
            {
                set_iwg_dof_dependencies( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_dv_dependencies( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_field_types( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_properties( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_material_models( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_constitutive_models( tIWGParameter, tIWG, tLeaderFollowerType );
            }
            set_iwg_stabilization_parameters( tIWGParameter, tIWG );
            mIWGs( iIWG ) = tIWG;
        }
    }

    void Model_Initializer_Legacy::set_iwg_residual_dof_type( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG ) const
    {
        std::string const tDofResidualString = aIWGParameter.get< std::string >( "dof_residual" );
        auto              tResDofTypes       = string_to_cell_of_cell< MSI::Dof_Type >( tDofResidualString, mMSIDofTypeMap );
        aIWG->set_residual_dof_type( tResDofTypes );
    }

    void Model_Initializer_Legacy::set_iwg_stabilization_parameters( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG )
    {
        auto tSPNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( "stabilization_parameters" ) );

        for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
        {
            // if CM name is in the CM map
            if ( mStabilizationParameterMap.find( tSPNamesPair( iSP )( 0 ) ) != mStabilizationParameterMap.end() )
            {
                // get SP index
                uint const tSPIndex = mStabilizationParameterMap[ tSPNamesPair( iSP )( 0 ) ];

                // set SP for IWG
                aIWG->set_stabilization_parameter(
                        mStabilizationParameters( tSPIndex ),
                        tSPNamesPair( iSP )( 1 ) );
            }
            else
            {
                // error message unknown SP
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown aSPString: %s \n",
                        tSPNamesPair( iSP )( 0 ).c_str() );
            }
        }
    }

    void Model_Initializer_Legacy::set_iwg_constitutive_models( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tCMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_constitutive_models" ) );

        for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
        {
            // if CM name is in the CM map
            if ( mConstitutiveModelMap.find( tCMNamesPair( iCM )( 0 ) ) != mConstitutiveModelMap.end() )
            {
                // get CM index
                uint const tCMIndex = mConstitutiveModelMap[ tCMNamesPair( iCM )( 0 ) ];

                // set CM for IWG
                aIWG->set_constitutive_model(
                        mConstitutiveModels( tCMIndex ),
                        tCMNamesPair( iCM )( 1 ),
                        aLeaderFollowerType );
            }
            else
            {
                // error message unknown CM
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aCMString: %s \n",
                        tPrefix.c_str(),
                        tCMNamesPair( iCM )( 0 ).c_str() );
            }
        }
    }

    void Model_Initializer_Legacy::set_iwg_material_models( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tMMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_material_model" ) );

        for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
        {
            // if MM name is in the CM map
            if ( mMaterialModelMap.find( tMMNamesPair( iMM )( 0 ) ) != mMaterialModelMap.end() )
            {
                // get MM index
                uint const tMMIndex = mMaterialModelMap[ tMMNamesPair( iMM )( 0 ) ];

                // set MM for IWG
                aIWG->set_material_model(
                        mMaterialModels( tMMIndex ),
                        tMMNamesPair( iMM )( 1 ),
                        aLeaderFollowerType );
            }
            else
            {
                // error message unknown MM
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aMMString: %s \n",
                        tPrefix.c_str(),
                        tMMNamesPair( iMM )( 0 ).c_str() );
            }
        }
    }

    void Model_Initializer_Legacy::set_iwg_properties( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix = mtk::get_leader_follower_string( aLeaderFollowerType );

        auto tPropertyNamesPair = string_to_cell_of_cell< std::string >(
                aIWGParameter.get< std::string >( tPrefix + "_properties" ) );

        for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
        {
            // if property name is in the property map
            if ( mPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != mPropertyMap.end() )
            {
                // get property index
                uint const tPropertyIndex = mPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                // set property for IWG
                aIWG->set_property(
                        mProperties( tPropertyIndex ),
                        tPropertyNamesPair( iProp )( 1 ),
                        aLeaderFollowerType );
            }
            else
            {
                // create error message unknown property
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aPropertyString: %s \n",
                        tPrefix.c_str(),
                        tPropertyNamesPair( iProp )( 0 ).c_str() );
            }
        }
    }

    void Model_Initializer_Legacy::set_iwg_field_types( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix     = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tFieldTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_field_types", mFieldTypeMap );
        aIWG->set_field_type_list( tFieldTypes, aLeaderFollowerType );
    }

    void Model_Initializer_Legacy::set_iwg_dof_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower aLeaderFollowerType ) const
    {
        // get the prefix of the property based on the leader or follower type
        std::string const tPrefix   = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tDofTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_dof_dependencies", mMSIDofTypeMap );
        aIWG->set_dof_type_list( tDofTypes, aLeaderFollowerType );
    }

    void Model_Initializer_Legacy::set_iwg_dv_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG, mtk::Leader_Follower const &aLeaderFollowerType ) const
    {
        // get the prefix of the property based on the leader or follower type
        std::string const tPrefix  = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tDvTypes = property_to_cell_of_cell( aIWGParameter, tPrefix + "_dv_dependencies", mMSIDvTypeMap );
        aIWG->set_dv_type_list( tDvTypes, aLeaderFollowerType );
    }

    // TODO: This method should be refactored: Too long and repeating code!
    void Model_Initializer_Legacy::create_iqis()
    {
        // create an IQI factory
        IQI_Factory tIQIFactory;

        // get the IQI parameter list
        moris::Cell< ParameterList > tIQIParameterList = mParameterList( 4 );

        // get number of IQIs
        uint tNumIQIs = tIQIParameterList.size();

        // set size for list of IQI pointers
        mIQIs.resize( tNumIQIs, nullptr );

        // loop over the parameter lists
        for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
        {
            // get the treated IQI parameter list
            ParameterList tIQIParameter = tIQIParameterList( iIQI );

            // get name from parameter list
            std::string tIQIName = tIQIParameter.get< std::string >( "IQI_name" );

            // get the IQI type from parameter list
            fem::IQI_Type tIQIType =
                    static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );

            // get the treated IQI bulk type
            fem::Element_Type tIQIBulkType =
                    static_cast< fem::Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );

            // create an IQI pointer
            mIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

            // set name
            mIQIs( iIQI )->set_name( tIQIName );

            mIQIs( iIQI )->set_bulk_type( tIQIBulkType );

            // fill IQI map
            mIQIMap[ tIQIName ] = iIQI;

            // get the treated IQI quantity dof type
            moris::Cell< moris::MSI::Dof_Type > tQuantityDofTypes;
            string_to_cell(
                    tIQIParameter.get< std::string >( "dof_quantity" ),
                    tQuantityDofTypes,
                    mMSIDofTypeMap );
            mIQIs( iIQI )->set_quantity_dof_type( tQuantityDofTypes );

            // set index for vectorial field
            mIQIs( iIQI )->set_output_type_index(
                    tIQIParameter.get< moris::sint >( "vectorial_field_index" ) );

            // set function parameters
            moris::Cell< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tIQIParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );
            mIQIs( iIQI )->set_parameters( tFuncParameters );

            // init string for leader or follower
            std::string          tIsLeaderString = "leader";
            mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

            // loop on leader and follower
            for ( uint iLeader = 0; iLeader <= 1; iLeader++ )
            {
                // if follower
                if ( iLeader )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // set dof dependencies
                moris::Cell< moris::Cell< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_dof_dependencies" ),
                        tDofTypes,
                        mMSIDofTypeMap );
                mIQIs( iIQI )->set_dof_type_list( tDofTypes, tIsLeader );

                // set dv dependencies
                moris::Cell< moris::Cell< PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_dv_dependencies" ),
                        tDvTypes,
                        mMSIDvTypeMap );
                mIQIs( iIQI )->set_dv_type_list( tDvTypes, tIsLeader );

                // set field types
                moris::Cell< moris::Cell< moris::mtk::Field_Type > > tFieldTypes;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_field_types" ),
                        tFieldTypes,
                        mFieldTypeMap );
                mIQIs( iIQI )->set_field_type_list( tFieldTypes, tIsLeader );

                // set properties
                moris::Cell< moris::Cell< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_properties" ),
                        tPropertyNamesPair );

                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // if property name is in the property map
                    if ( mPropertyMap.find( tPropertyNamesPair( iProp )( 0 ) ) != mPropertyMap.end() )
                    {
                        // get property index
                        uint tPropertyIndex = mPropertyMap[ tPropertyNamesPair( iProp )( 0 ) ];

                        // set property for IWG
                        mIQIs( iIQI )->set_property(
                                mProperties( tPropertyIndex ),
                                tPropertyNamesPair( iProp )( 1 ),
                                tIsLeader );
                    }
                    else
                    {
                        // error message unknown property
                        MORIS_ERROR( false,
                                "Model_Initializer::create_IQIs_without_phase - Unknown %s aPropertyString: %s \n",
                                tIsLeaderString.c_str(),
                                tPropertyNamesPair( iProp )( 0 ).c_str() );
                    }
                }

                // set constitutive models
                moris::Cell< moris::Cell< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // if CM name is in the CM map
                    if ( mConstitutiveModelMap.find( tCMNamesPair( iCM )( 0 ) ) != mConstitutiveModelMap.end() )
                    {
                        // get CM index
                        uint tCMIndex = mConstitutiveModelMap[ tCMNamesPair( iCM )( 0 ) ];

                        // set CM for IQI
                        mIQIs( iIQI )->set_constitutive_model(
                                mConstitutiveModels( tCMIndex ),
                                tCMNamesPair( iCM )( 1 ),
                                tIsLeader );
                    }
                    else
                    {
                        // error message unknown CM
                        MORIS_ERROR( false,
                                "Model_Initializer::create_IQIs_without_phase - Unknown %s aCMString: %s \n",
                                tIsLeaderString.c_str(),
                                tCMNamesPair( iCM )( 0 ).c_str() );
                    }
                }
            }

            // set stabilization parameters
            moris::Cell< moris::Cell< std::string > > tSPNamesPair;
            string_to_cell_of_cell(
                    tIQIParameter.get< std::string >( "stabilization_parameters" ),
                    tSPNamesPair );

            for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
            {
                // if SP name is in the SP map
                if ( mStabilizationParameterMap.find( tSPNamesPair( iSP )( 0 ) ) != mStabilizationParameterMap.end() )
                {
                    // get SP index
                    uint tSPIndex = mStabilizationParameterMap[ tSPNamesPair( iSP )( 0 ) ];

                    // set SP for IQI
                    mIQIs( iIQI )->set_stabilization_parameter(
                            mStabilizationParameters( tSPIndex ),
                            tSPNamesPair( iSP )( 1 ) );
                }
                else
                {
                    // error message unknown SP
                    MORIS_ERROR( false,
                            "Model_Initializer::create_IQIs_without_phase - Unknown aSPString: %s \n",
                            tSPNamesPair( iSP )( 0 ).c_str() );
                }
            }

            // debug - uncomment if needed to list IQIs the code actually sees
            // mIQIs( iIQI )->print_names();
        }
    }

    void Model_Initializer_Legacy::create_set_info()
    {
        // get fem computation type parameter list
        ParameterList const tComputationParameterList = mParameterList( 5 )( 0 );

        bool const tIsAnalyticalSA = tComputationParameterList.get< bool >( "is_analytical_sensitivity" );
        auto const tFDSchemeForSA  = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme" ) );
        real const tFDPerturbation = tComputationParameterList.get< real >( "finite_difference_perturbation_size" );

        // create a map of the set
        std::map< std::tuple< std::string, bool, bool >, uint > tMeshToFemSetIndex;
        this->create_fem_set_info_from_iwgs( tIsAnalyticalSA, tFDSchemeForSA, tFDPerturbation, tMeshToFemSetIndex );
        this->create_fem_set_info_from_iqis( tIsAnalyticalSA, tFDSchemeForSA, tFDPerturbation, tMeshToFemSetIndex );
    }

    // TODO: This method has a similar logic as create_fem_set_info_from_iqis. Proper refactoring could merge the two methods.
    void Model_Initializer_Legacy::create_fem_set_info_from_iwgs(
            bool                                                     aIsAnalyticalSA,
            FDScheme_Type const                                     &aFDSchemeForSA,
            real const                                               aFDPerturbation,
            std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )
    {
        moris::Cell< ParameterList > tIWGParameterLists = this->mParameterList( 3 );
        for ( uint iIWG = 0; iIWG < tIWGParameterLists.size(); iIWG++ )
        {
            ParameterList const          &tIWGParameterList = tIWGParameterLists( iIWG );
            std::shared_ptr< IWG > const &tIWG              = this->mIWGs( iIWG );

            // get the time continuity and time boundary flags from the IWG parameter list to uniquely identify the fem sets
            bool const tTimeContinuity = tIWGParameterList.get< bool >( "time_continuity" );
            bool const tTimeBoundary   = tIWGParameterList.get< bool >( "time_boundary" );

            // get a representative DoF type
            MSI::Dof_Type const tFirstResidualDofType = tIWG->get_residual_dof_type()( 0 )( 0 );

            // loop over the mesh set names
            auto tMeshSetNames = string_to_cell< std::string >( tIWGParameterList.get< std::string >( "mesh_set_names" ) );
            for ( auto &tMeshSetName : tMeshSetNames )
            {
                // check for ghost set names and select correct B-spline mesh automatically when new ghost sets need to be used
                this->check_and_set_ghost_set_names( tMeshSetName, tFirstResidualDofType );

                // create a tuple with the mesh set name, time continuity and time boundary flags that will be used to
                // uniquely identify the fem set
                auto tMeshTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );

                // check if the mesh set name already in map
                if ( aMeshToFemSetIndex.find( tMeshTuple ) == aMeshToFemSetIndex.end() )
                {
                    // if the set did not yet exist, create a new one and add keep track of the index
                    aMeshToFemSetIndex[ tMeshTuple ] = this->mSetInfo.size();

                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );
                    aSetUserInfo.set_time_continuity( tTimeContinuity );
                    aSetUserInfo.set_time_boundary( tTimeBoundary );
                    aSetUserInfo.set_is_analytical_sensitivity_analysis( aIsAnalyticalSA );
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( aFDSchemeForSA );
                    aSetUserInfo.set_finite_difference_perturbation_size( aFDPerturbation );
                    aSetUserInfo.add_IWG( tIWG );
                    this->mSetInfo.push_back( aSetUserInfo );
                }
                else
                {
                    // if the fem set already exister, only add the IWG
                    this->mSetInfo( aMeshToFemSetIndex[ tMeshTuple ] ).add_IWG( tIWG );
                }
            }
        }
    }

    // TODO: This method has a similar logic as create_fem_set_info_from_iwgs. Proper refactoring could merge the two methods.
    void Model_Initializer_Legacy::create_fem_set_info_from_iqis(
            bool const                                               aIsAnalyticalSA,
            FDScheme_Type const                                     &aFDSchemeForSA,
            real const                                               aFDPerturbation,
            std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )
    {
        moris::Cell< ParameterList > tIQIParameterLists = this->mParameterList( 4 );
        for ( uint iIQI = 0; iIQI < tIQIParameterLists.size(); iIQI++ )
        {
            ParameterList const          &tIQIParameterList = tIQIParameterLists( iIQI );
            std::shared_ptr< IQI > const &tIQI              = this->mIQIs( iIQI );

            bool tTimeContinuity = tIQIParameterList.get< bool >( "time_continuity" );
            bool tTimeBoundary   = tIQIParameterList.get< bool >( "time_boundary" );

            // loop over the mesh set names
            auto tMeshSetNames = string_to_cell< std::string >( tIQIParameterList.get< std::string >( "mesh_set_names" ) );
            for ( auto &tMeshSetName : tMeshSetNames )
            {
                auto tMeshTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );

                // if the mesh set name not in map
                if ( aMeshToFemSetIndex.find( tMeshTuple ) == aMeshToFemSetIndex.end() )
                {
                    // add the mesh set name map
                    aMeshToFemSetIndex[ tMeshTuple ] = this->mSetInfo.size();

                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );
                    aSetUserInfo.set_time_continuity( tTimeContinuity );
                    aSetUserInfo.set_time_boundary( tTimeBoundary );
                    aSetUserInfo.set_is_analytical_sensitivity_analysis( aIsAnalyticalSA );
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( aFDSchemeForSA );
                    aSetUserInfo.set_finite_difference_perturbation_size( aFDPerturbation );
                    aSetUserInfo.add_IQI( tIQI );

                    // add it to the list of fem set info
                    this->mSetInfo.push_back( aSetUserInfo );
                }
                else
                {
                    // set the IQI
                    this->mSetInfo( aMeshToFemSetIndex[ tMeshTuple ] ).add_IQI( tIQI );
                }
            }
        }
    }


}    // namespace moris::fem