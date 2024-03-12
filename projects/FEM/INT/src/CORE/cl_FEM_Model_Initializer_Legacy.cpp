//
// Created by frank on 12/14/23.
//

#include "cl_FEM_Model_Initializer_Legacy.hpp"
#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_Vector.hpp"
#include "cl_Param_List.hpp"
#include "GEN_Data_Types.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MSI_Dof_Type_Enums.hpp"
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
        Vector< ParameterList > tMMParameterList = mParameterList( 7 );

        // loop over the parameter lists for MM
        for ( auto const &tMMParameter : tMMParameterList )
        {
            // get the material type from parameter list
            fem::Material_Type                tMMType = static_cast< fem::Material_Type >( tMMParameter.get< uint >( "material_type" ) );
            std::shared_ptr< Material_Model > tMM     = tMMFactory.create_MM( tMMType );

            std::string tMMName = tMMParameter.get< std::string >( "material_name" );
            tMM->set_name( tMMName );
            mMaterialModels[ tMMName ] = tMM;

            // set MM space dimension
            tMM->set_space_dim( mSpatialDimension );

            // set MM dof dependencies
            Vector< Vector< moris::MSI::Dof_Type > > tDofTypes = string_to_cell_of_cell< moris::MSI::Dof_Type >(
                    std::get< 0 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),

                    mMSIDofTypeMap );
            Vector< std::string > tDofTypeNames;
            string_to_cell(
                    std::get< 1 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypeNames );
            tMM->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set MM properties
            for ( auto const &[ tProperty, tPropertyString ] : read_properties( tMMParameter, mtk::Leader_Follower::UNDEFINED ) )
            {
                tMM->set_property( tProperty, tPropertyString );
            }
            // set local properties
            tMM->set_local_properties();
        }
    }

    void Model_Initializer_Legacy::create_constitutive_models()
    {
        // create a constitutive model factory
        CM_Factory tCMFactory;

        // get the CM parameter list
        Vector< ParameterList > tCMParameterList = mParameterList( 1 );

        // loop over the parameter lists for CM
        for ( auto const &tCMParameter : tCMParameterList )
        {
            // get the constitutive type from parameter list
            fem::Constitutive_Type tCMType =
                    static_cast< fem::Constitutive_Type >( tCMParameter.get< uint >( "constitutive_type" ) );

            // create a constitutive model pointer
            std::shared_ptr< Constitutive_Model > tCM     = tCMFactory.create_CM( tCMType );
            std::string const                     tCMName = tCMParameter.get< std::string >( "constitutive_name" );
            tCM->set_name( tCMName );
            mConstitutiveModels[ tCMName ] = tCM;

            // set CM model type
            fem::Model_Type tCMModelType =
                    static_cast< fem::Model_Type >( tCMParameter.get< uint >( "model_type" ) );
            if ( tCMModelType != fem::Model_Type::UNDEFINED )
            {
                tCM->set_model_type( tCMModelType );
            }

            // set CM space dimension
            tCM->set_space_dim( mSpatialDimension );

            // set CM dof dependencies
            Vector< Vector< moris::MSI::Dof_Type > > tDofTypes;
            string_to_cell_of_cell(
                    std::get< 0 >( tCMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypes,
                    mMSIDofTypeMap );
            Vector< std::string > tDofTypeNames;
            string_to_cell(
                    std::get< 1 >( tCMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypeNames );
            tCM->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set CM dv dependencies
            Vector< Vector< gen::PDV_Type > > tDvTypes;
            string_to_cell_of_cell(
                    std::get< 0 >( tCMParameter.get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                    tDvTypes,
                    mMSIDvTypeMap );
            Vector< std::string > tDvTypeNames;
            string_to_cell(
                    std::get< 1 >( tCMParameter.get< std::pair< std::string, std::string > >( "dv_dependencies" ) ),
                    tDvTypeNames );
            tCM->set_dv_type_list( tDvTypes, tDvTypeNames );

            // set CM properties
            for ( auto const &[ tProperty, tPropertyString ] : read_properties( tCMParameter, mtk::Leader_Follower::UNDEFINED ) )
            {
                tCM->set_property( tProperty, tPropertyString );
            }
            // set local properties
            tCM->set_local_properties();

            // set material model
            Vector< Vector< std::string > > tMMNamesPair;
            string_to_cell_of_cell(
                    tCMParameter.get< std::string >( "material_model" ),
                    tMMNamesPair );

            for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
            {
                std::string const tMMName = tMMNamesPair( iMM )( 0 );
                std::string const tMMType = tMMNamesPair( iMM )( 1 );
                // if MM name is in the CM map
                if ( mMaterialModels.find( tMMName ) != mMaterialModels.end() )
                {
                    // set MM for IWG
                    tCM->set_material_model(
                            mMaterialModels[ tMMName ],
                            tMMType );
                }
                else
                {
                    // error message unknown MM
                    MORIS_ERROR( false,
                            "Model_Initializer::create_CMs - Unknown aMMString: %s \n",
                            tMMName.c_str() );
                }
            }
        }
    }

    void Model_Initializer_Legacy::create_stabilization_parameters()
    {
        // create a stabilization parameter factory
        SP_Factory tSPFactory;

        // get the SP parameter list
        Vector< ParameterList > tSPParameterList = mParameterList( 2 );

        // loop over the parameter list
        for ( auto const &tSPParameter : tSPParameterList )
        {
            // get the stabilization type from parameter list
            fem::Stabilization_Type tSPType =
                    static_cast< fem::Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );
            std::shared_ptr< Stabilization_Parameter > tSP     = tSPFactory.create_SP( tSPType );
            std::string const                          tSPName = tSPParameter.get< std::string >( "stabilization_name" );
            tSP->set_name( tSPName );
            mStabilizationParameters[ tSPName ] = tSP;

            // set SP space dimension
            tSP->set_space_dim( mSpatialDimension );

            // set parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tSPParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );
            tSP->set_parameters( tFuncParameters );

            // init string for leader or follower
            std::string          tIsLeaderString = "leader";
            mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

            // loop on leader and follower
            for ( uint iLeader = 0; iLeader <= tSP->get_has_follower(); iLeader++ )
            {
                // if follower
                if ( iLeader )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // set dof dependencies
                Vector< Vector< moris::MSI::Dof_Type > > tDofTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                        tDofTypes,
                        mMSIDofTypeMap );
                Vector< std::string > tDofTypeNames;
                string_to_cell( std::get< 1 >(
                                        tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dof_dependencies" ) ),
                        tDofTypeNames );
                tSP->set_dof_type_list( tDofTypes, tDofTypeNames, tIsLeader );

                // set dv dependencies
                Vector< Vector< gen::PDV_Type > > tDvTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                        tDvTypes,
                        mMSIDvTypeMap );
                Vector< std::string > tDvTypeNames;
                string_to_cell(
                        std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( tIsLeaderString + "_dv_dependencies" ) ),
                        tDvTypeNames );
                tSP->set_dv_type_list( tDvTypes, tDvTypeNames, tIsLeader );

                // set leader properties
                for ( auto const &[ tProperty, tPropertyString ] : read_properties( tSPParameter, tIsLeader ) )
                {
                    tSP->set_property( tProperty, tPropertyString, tIsLeader );
                }

                // set constitutive models
                Vector< Vector< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tSPParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                // loop over the CM names
                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // get the CM name
                    std::string const tCMName = tCMNamesPair( iCM )( 0 );
                    std::string const tCMType = tCMNamesPair( iCM )( 1 );
                    // check for unknown CM
                    MORIS_ERROR( mConstitutiveModels.find( tCMName ) != mConstitutiveModels.end(),
                            "Model_Initializer::create_stabilization_parameters_without_phase - Unknown leader aCMString: %s \n",
                            tCMName.c_str() );

                    // set CM for SP
                    tSP->set_constitutive_model( mConstitutiveModels[ tCMName ], tCMType );
                }
            }
        }
    }

    void Model_Initializer_Legacy::create_iwgs()
    {
        IWG_Factory             tIWGFactory;
        Vector< ParameterList > tIWGParameterList = mParameterList( 3 );

        for ( auto const &tIWGParameter : tIWGParameterList )
        {
            fem::IWG_Type const tIWGType =
                    static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

            std::shared_ptr< IWG > tIWG = tIWGFactory.create_IWG( tIWGType );

            std::string const tIWGName = tIWGParameter.get< std::string >( "IWG_name" );
            tIWG->set_name( tIWGName );
            mIWGs[ tIWGName ] = tIWG;

            tIWG->set_is_fd_jacobian( not tIWGParameter.get< bool >( "analytical_jacobian" ) ); // warning: this is a negation

            // get function parameters
            auto tFuncParameters = string_to_cell_mat_2< DDRMat >( tIWGParameter.get< std::string >( "function_parameters" ) );
            tIWG->set_parameters( tFuncParameters );

            // get the ghost order from parameter list
            uint const tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );
            tIWG->set_interpolation_order( tGhostOrder );

            // set residual dof type
            set_iwg_residual_dof_type( tIWGParameter, tIWG );

            // loop over leader and follower and set the appropriate properties
            Vector< mtk::Leader_Follower > const tLeaderFollower{ mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER };
            for ( auto const &tLeaderFollowerType : tLeaderFollower )
            {
                set_iwg_iqi_dof_dependencies( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_iqi_dv_dependencies( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_iqi_field_types( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_iqi_properties( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_iqi_material_models( tIWGParameter, tIWG, tLeaderFollowerType );
                set_iwg_iqi_constitutive_models( tIWGParameter, tIWG, tLeaderFollowerType );
            }
            set_iwg_iqi_stabilization_parameters( tIWGParameter, tIWG );
        }
    }

    // TODO: This method should be refactored: Too long and repeating code!
    void Model_Initializer_Legacy::create_iqis()
    {
        // create an IQI factory
        IQI_Factory tIQIFactory;

        // get the IQI parameter list
        Vector< ParameterList > tIQIParameterList = mParameterList( 4 );

        // loop over the parameter lists
        for ( auto const &tIQIParameter : tIQIParameterList )
        {
            fem::IQI_Type          tIQIType = static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );
            std::shared_ptr< IQI > tIQI     = tIQIFactory.create_IQI( tIQIType );

            std::string tIQIName = tIQIParameter.get< std::string >( "IQI_name" );
            tIQI->set_name( tIQIName );
            mIQIs[ tIQIName ] = tIQI;

            fem::Element_Type tIQIBulkType = static_cast< fem::Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );
            tIQI->set_bulk_type( tIQIBulkType );

            // get the treated IQI quantity dof type
            Vector< moris::MSI::Dof_Type > tQuantityDofTypes;
            string_to_cell(
                    tIQIParameter.get< std::string >( "dof_quantity" ),
                    tQuantityDofTypes,
                    mMSIDofTypeMap );
            tIQI->set_quantity_dof_type( tQuantityDofTypes );

            tIQI->set_output_type_index( tIQIParameter.get< moris::sint >( "vectorial_field_index" ) );

            // set function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tIQIParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );
            tIQI->set_parameters( tFuncParameters );

            tIQI->set_normalization_type( tIQIParameter.get< std::string >( "normalization" ) );

            // loop over leader and follower and set the appropriate properties
            Vector< mtk::Leader_Follower > const tLeaderFollower{ mtk::Leader_Follower::LEADER, mtk::Leader_Follower::FOLLOWER };
            for ( auto const &tLeaderFollowerType : tLeaderFollower )
            {
                set_iwg_iqi_dof_dependencies( tIQIParameter, tIQI, tLeaderFollowerType );
                set_iwg_iqi_dv_dependencies( tIQIParameter, tIQI, tLeaderFollowerType );
                set_iwg_iqi_field_types( tIQIParameter, tIQI, tLeaderFollowerType );
                set_iwg_iqi_properties( tIQIParameter, tIQI, tLeaderFollowerType );
                set_iwg_iqi_constitutive_models( tIQIParameter, tIQI, tLeaderFollowerType );
            }
            set_iwg_iqi_stabilization_parameters( tIQIParameter, tIQI );
        }
    }

    void Model_Initializer_Legacy::set_iwg_residual_dof_type( ParameterList const &aIWGParameter, std::shared_ptr< IWG > &aIWG ) const
    {
        std::string const tDofResidualString = aIWGParameter.get< std::string >( "dof_residual" );
        auto              tResDofTypes       = string_to_cell_of_cell< MSI::Dof_Type >( tDofResidualString, mMSIDofTypeMap );
        aIWG->set_residual_dof_type( tResDofTypes );
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_stabilization_parameters( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi )
    {
        auto tSPNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( "stabilization_parameters" ) );

        for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
        {
            std::string const tSPName = tSPNamesPair( iSP )( 0 );
            std::string const tSPType = tSPNamesPair( iSP )( 1 );
            // if CM name is in the CM map
            if ( mStabilizationParameters.find( tSPName ) != mStabilizationParameters.end() )
            {
                // set SP for IWG
                aIwgIqi->set_stabilization_parameter( mStabilizationParameters[ tSPName ], tSPType );
            }
            else
            {
                // error message unknown SP
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown aSPString: %s \n",
                        tSPName.c_str() );
            }
        }
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_constitutive_models( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tCMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_constitutive_models" ) );

        for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
        {
            std::string const tCMName = tCMNamesPair( iCM )( 0 );
            std::string const tCMType = tCMNamesPair( iCM )( 1 );

            // if CM name is in the CM map
            if ( mConstitutiveModels.find( tCMName ) != mConstitutiveModels.end() )
            {
                // set CM for IWG
                aIwgIqi->set_constitutive_model( mConstitutiveModels[ tCMName ], tCMType, aLeaderFollowerType );
            }
            else
            {
                // error message unknown CM
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aCMString: %s \n",
                        tPrefix.c_str(),
                        tCMName.c_str() );
            }
        }
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_material_models( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix      = mtk::get_leader_follower_string( aLeaderFollowerType );
        auto              tMMNamesPair = string_to_cell_of_cell< std::string >( aIWGParameter.get< std::string >( tPrefix + "_material_model" ) );

        for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
        {
            std::string const tMMName = tMMNamesPair( iMM )( 0 );
            std::string const tMMType = tMMNamesPair( iMM )( 1 );
            // if MM name is in the CM map
            if ( mMaterialModels.find( tMMName ) != mMaterialModels.end() )
            {
                // set MM for IWG
                aIwgIqi->set_material_model(
                        mMaterialModels[ tMMName ],
                        tMMType,
                        aLeaderFollowerType );
            }
            else
            {
                // error message unknown MM
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aMMString: %s \n",
                        tPrefix.c_str(),
                        tMMName.c_str() );
            }
        }
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_properties( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower const &aLeaderFollowerType )
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const tPrefix = mtk::get_leader_follower_string( aLeaderFollowerType );

        auto tPropertyNamesPair = string_to_cell_of_cell< std::string >(
                aIWGParameter.get< std::string >( tPrefix + "_properties" ) );

        for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
        {
            std::string const tPropertyName = tPropertyNamesPair( iProp )( 0 );
            std::string const tPropertyType = tPropertyNamesPair( iProp )( 1 );

            // if property name is in the property map
            if ( mProperties.find( tPropertyName ) != mProperties.end() )
            {
                // set property for IWG
                aIwgIqi->set_property( mProperties[ tPropertyName ], tPropertyType, aLeaderFollowerType );
            }
            else
            {
                // create error message unknown property
                MORIS_ERROR( false,
                        "FEM_Model::create_IWGs_without_phase - Unknown %s aPropertyString: %s \n",
                        tPrefix.c_str(),
                        tPropertyName.c_str() );
            }
        }
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_field_types( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower const &aLeaderFollowerType ) const
    {
        // get the prefix of the property name based on the leader or follower type (either "leader" or "follower")
        std::string const                          tPrefix     = mtk::get_leader_follower_string( aLeaderFollowerType );
        Vector< Vector< moris::mtk::Field_Type > > tFieldTypes = parameter_to_vec_of_vec( aIWGParameter, tPrefix + "_field_types", mMTKFieldTypeMap );
        aIwgIqi->set_field_type_list( tFieldTypes, aLeaderFollowerType );
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_dof_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower aLeaderFollowerType ) const
    {
        // get the prefix of the property based on the leader or follower type
        std::string const                        tPrefix   = mtk::get_leader_follower_string( aLeaderFollowerType );
        Vector< Vector< moris::MSI::Dof_Type > > tDofTypes = parameter_to_vec_of_vec( aIWGParameter, tPrefix + "_dof_dependencies", mMSIDofTypeMap );
        aIwgIqi->set_dof_type_list( tDofTypes, aLeaderFollowerType );
    }

    template<typename T>
    void Model_Initializer_Legacy::set_iwg_iqi_dv_dependencies( ParameterList const &aIWGParameter, std::shared_ptr< T > &aIwgIqi, mtk::Leader_Follower const &aLeaderFollowerType ) const
    {
        // get the prefix of the property based on the leader or follower type
        std::string const                        tPrefix  = mtk::get_leader_follower_string( aLeaderFollowerType );
        Vector< Vector< moris::gen::PDV_Type > > tDvTypes = parameter_to_vec_of_vec( aIWGParameter, tPrefix + "_dv_dependencies", mMSIDvTypeMap );
        aIwgIqi->set_dv_type_list( tDvTypes, aLeaderFollowerType );
    }

    void Model_Initializer_Legacy::create_set_info()
    {
        // create a map of the set
        std::map< std::tuple< std::string, bool, bool >, uint > tMeshToFemSetIndex;
        this->create_fem_set_info_from_iwgs( tMeshToFemSetIndex );
        this->create_fem_set_info_from_iqis( tMeshToFemSetIndex );
    }

    // TODO: This method has a similar logic as create_fem_set_info_from_iqis. Proper refactoring could merge the two methods.
    void Model_Initializer_Legacy::create_fem_set_info_from_iwgs( std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )
    {
        ParameterList const tComputationParameterList = mParameterList( 5 )( 0 );

        // forward analysis
        bool const         tIsAnalyticalFA   = tComputationParameterList.get< bool >( "is_analytical_forward" );
        fem::FDScheme_Type tFDSchemeForFA    = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme_forward" ) );
        real               tFDPerturbationFA = tComputationParameterList.get< real >( "finite_difference_perturbation_size_forward" );

        // sensitivity analysis
        bool const tIsAnalyticalSA = tComputationParameterList.get< bool >( "is_analytical_sensitivity" );
        auto const tFDSchemeForSA  = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme" ) );
        real const tFDPerturbation = tComputationParameterList.get< real >( "finite_difference_perturbation_size" );

        fem::Perturbation_Type tPerturbationStrategy = static_cast< fem::Perturbation_Type >( tComputationParameterList.get< uint >( "finite_difference_perturbation_strategy" ) );

        Vector< ParameterList > tIWGParameterLists = this->mParameterList( 3 );
        for ( auto const &tIWGParameterList : tIWGParameterLists )
        {
            std::string                   tIWGName = tIWGParameterList.get< std::string >( "IWG_name" );
            std::shared_ptr< IWG > const &tIWG     = this->mIWGs[ tIWGName ];

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

                auto tSetTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );
                // check if the mesh set name already in map
                if ( mSetInfo.find( tSetTuple ) == mSetInfo.end() )
                {
                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );
                    aSetUserInfo.set_time_continuity( tTimeContinuity );
                    aSetUserInfo.set_time_boundary( tTimeBoundary );

                    aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );
                    aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );
                    aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                    aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                    aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                    aSetUserInfo.add_IWG( tIWG );
                    mSetInfo[ tSetTuple ] = aSetUserInfo;
                }
                else
                {
                    // add the current iwg to the existing fem set info
                    mSetInfo[ tSetTuple ].add_IWG( tIWG );
                }
            }
        }
    }

    // TODO: This method has a similar logic as create_fem_set_info_from_iwgs. Proper refactoring could merge the two methods.
    void Model_Initializer_Legacy::create_fem_set_info_from_iqis( std::map< std::tuple< std::string, bool, bool >, uint > &aMeshToFemSetIndex )
    {
        ParameterList const tComputationParameterList = mParameterList( 5 )( 0 );

        // forward analysis
        bool const         tIsAnalyticalFA   = tComputationParameterList.get< bool >( "is_analytical_forward" );
        fem::FDScheme_Type tFDSchemeForFA    = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme_forward" ) );
        real               tFDPerturbationFA = tComputationParameterList.get< real >( "finite_difference_perturbation_size_forward" );

        // sensitivity analysis
        bool const tIsAnalyticalSA = tComputationParameterList.get< bool >( "is_analytical_sensitivity" );
        auto const tFDSchemeForSA  = static_cast< fem::FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme" ) );
        real const tFDPerturbation = tComputationParameterList.get< real >( "finite_difference_perturbation_size" );

        fem::Perturbation_Type tPerturbationStrategy = static_cast< fem::Perturbation_Type >( tComputationParameterList.get< uint >( "finite_difference_perturbation_strategy" ) );

        Vector< ParameterList > tIQIParameterLists = this->mParameterList( 4 );
        for ( auto const &tIQIParameter : tIQIParameterLists )
        {
            std::string                   tIQIName = tIQIParameter.get< std::string >( "IQI_name" );
            std::shared_ptr< IQI > const &tIQI     = this->mIQIs[ tIQIName ];

            bool tTimeContinuity = tIQIParameter.get< bool >( "time_continuity" );
            bool tTimeBoundary   = tIQIParameter.get< bool >( "time_boundary" );

            // loop over the mesh set names
            auto tMeshSetNames = string_to_cell< std::string >( tIQIParameter.get< std::string >( "mesh_set_names" ) );
            for ( auto &tMeshSetName : tMeshSetNames )
            {
                auto tSetTuple = std::make_tuple( tMeshSetName, tTimeContinuity, tTimeBoundary );
                if ( mSetInfo.find( tSetTuple ) == mSetInfo.end() )    // check if the mesh set name already in map
                {
                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );
                    aSetUserInfo.set_time_continuity( tTimeContinuity );
                    aSetUserInfo.set_time_boundary( tTimeBoundary );

                    aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );
                    aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );
                    aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                    aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbation );

                    aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                    aSetUserInfo.add_IQI( tIQI );
                    mSetInfo[ tSetTuple ] = aSetUserInfo;    // add it to the list of fem set info
                }
                else
                {
                    mSetInfo[ tSetTuple ].add_IQI( tIQI );    // add the current iqi to the existing fem set info
                }
            }
        }
    }


}    // namespace moris::fem