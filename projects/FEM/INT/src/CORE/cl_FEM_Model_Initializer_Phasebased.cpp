//
// Created by frank on 12/14/23.
//

#include "cl_FEM_Model_Initializer_Phasebased.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"

namespace moris::fem
{
    void Model_Initializer_Phasebased::initialize()
    {
        this->create_phases();
        Model_Initializer::initialize();
    }

    void Model_Initializer_Phasebased::create_phases()
    {
        // get the phase parameter list
        Vector< ParameterList > tPhaseParameterList = mParameterList( 7 );

        // get number of phases
        uint tNumPhases = tPhaseParameterList.size();

        // resize the list of phase user info
        mPhaseInfo.resize( tNumPhases );

        // loop over the parameter lists
        for ( uint iPhase = 0; iPhase < tNumPhases; iPhase++ )
        {
            // get the treated phase parameter list
            ParameterList tPhaseParameter = tPhaseParameterList( iPhase );

            // get the phase name from parameter list
            std::string tPhaseName =
                    tPhaseParameter.get< std::string >( "phase_name" );

            // set phase name to phase
            mPhaseInfo( iPhase ).set_phase_name( tPhaseName );

            // get the phase index from parameter list
            moris::Matrix< moris::IndexMat > tPhaseIndices;
            string_to_mat( tPhaseParameter.get< std::string >( "phase_indices" ), tPhaseIndices );

            // set phase mesh indices to phase
            mPhaseInfo( iPhase ).set_phase_indices( tPhaseIndices );

            // fill phase map
            mPhaseMap[ tPhaseName ] = iPhase;
        }
    }

    void Model_Initializer_Phasebased::create_material_models()
    {
        // create a constitutive model factory
        MM_Factory tMMFactory;

        // get the MM parameter list
        Vector< ParameterList > tMMParameterList = mParameterList( 8 );

        // get number of constitutive models
        uint tNumMMs = tMMParameterList.size();

        // loop over the parameter lists for MM
        for ( uint iMM = 0; iMM < tNumMMs; iMM++ )
        {
            // get the treated MM parameter list
            ParameterList tMMParameter = tMMParameterList( iMM );

            // get the constitutive type from parameter list
            fem::Material_Type tMMType =
                    static_cast< fem::Material_Type >( tMMParameter.get< uint >( "material_type" ) );

            // get the constitutive model name from parameter list
            std::string tMMName =
                    tMMParameter.get< std::string >( "material_name" );

            // get the phase from parameter list
            std::string tPhaseName =
                    tMMParameter.get< std::string >( "phase_name" );

            // create a constitutive model pointer
            std::shared_ptr< fem::Material_Model > tMM =
                    tMMFactory.create_MM( tMMType );

            // set MM name
            tMM->set_name( tMMName );

            // set MM space dimension
            tMM->set_space_dim( mSpatialDimension );

            // set MM dof dependencies
            Vector< Vector< moris::MSI::Dof_Type > > tDofTypes;
            string_to_cell_of_cell(
                    std::get< 0 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypes,
                    mMSIDofTypeMap );
            Vector< std::string > tDofTypeNames;
            string_to_cell(
                    std::get< 1 >( tMMParameter.get< std::pair< std::string, std::string > >( "dof_dependencies" ) ),
                    tDofTypeNames );
            tMM->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set MM properties
            Vector< Vector< std::string > > tPropertyNamesPair;
            string_to_cell_of_cell(
                    tMMParameter.get< std::string >( "properties" ),
                    tPropertyNamesPair );
            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                // get the property name
                std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                // check if property in the map
                MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                        "Model_Initializer_Phasebased::create_MMs - Unknown aPropertyString : %s \n",
                        tPropertyName.c_str() );

                // get property index
                uint tPropertyIndex = mPropertyMap[ tPropertyName ];

                // set property for MM
                tMM->set_property(
                        mProperties( tPropertyIndex ),
                        tPropertyNamesPair( iProp )( 1 ) );
            }

            // set local properties
            tMM->set_local_properties();

            // check the phase exist
            MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                    "Model_Initializer_Phasebased::create_material_models_without_phase - Unknown tPhaseName : %s \n",
                    tPhaseName.c_str() );

            // set MM to corresponding phase
            mPhaseInfo( mPhaseMap[ tPhaseName ] ).set_MM( tMM );
        }
    }

    void Model_Initializer_Phasebased::create_constitutive_models()
    {
        // create a constitutive model factory
        CM_Factory tCMFactory;

        // get the CM parameter list
        Vector< ParameterList > tCMParameterList = mParameterList( 1 );

        // get number of constitutive models
        uint tNumCMs = tCMParameterList.size();

        // loop over the parameter lists for CM
        for ( uint iCM = 0; iCM < tNumCMs; iCM++ )
        {
            // get the treated CM parameter list
            ParameterList tCMParameter = tCMParameterList( iCM );

            // get the constitutive type from parameter list
            fem::Constitutive_Type tCMType =
                    static_cast< fem::Constitutive_Type >( tCMParameter.get< uint >( "constitutive_type" ) );

            // get the constitutive model name from parameter list
            std::string tCMName =
                    tCMParameter.get< std::string >( "constitutive_name" );

            // get the phase from parameter list
            std::string tPhaseName =
                    tCMParameter.get< std::string >( "phase_name" );

            // check for unknown phase
            MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                    "Model_Initializer_Phasebased::create_CMs - Unknown phase name: %s \n",
                    tPhaseName.c_str() );

            // get the model type
            fem::Model_Type tCMModelType =
                    static_cast< fem::Model_Type >( tCMParameter.get< uint >( "model_type" ) );

            // create a constitutive model pointer
            std::shared_ptr< fem::Constitutive_Model > tCM =
                    tCMFactory.create_CM( tCMType );

            // set CM name
            tCM->set_name( tCMName );

            // set CM model type. must come before "set_space_dim"
            // fixme: currently cannot set a plane type and tensor type at the same time from an input file
            if ( tCMModelType != fem::Model_Type::UNDEFINED )
            {
                tCM->set_model_type( tCMModelType );
            }

            // set CM space dimension
            tCM->set_space_dim( mSpatialDimension );

            // set function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tCMParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );
            tCM->set_parameters( tFuncParameters );

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

            // set CM material model
            Vector< Vector< std::string > > tMMNamesPair;
            string_to_cell_of_cell(
                    tCMParameter.get< std::string >( "material_model" ),
                    tMMNamesPair );
            MORIS_ERROR( tMMNamesPair.size() <= 1, "Model_Initializer_Phasebased::create_CMs() - Only one material model per CM allowed." );

            // loop over Material Model names
            for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
            {
                // get the material name
                std::string tMaterialName = tMMNamesPair( iMM )( 0 );

                // get phase index
                uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                // get MM from phase
                std::shared_ptr< fem::Material_Model > tMM =
                        mPhaseInfo( tPhaseIndex ).get_MM_by_name( tMaterialName );

                // set material for CM
                tCM->set_material_model(
                        tMM,
                        tMMNamesPair( iMM )( 1 ) );
            }

            // set CM properties
            Vector< Vector< std::string > > tPropertyNamesPair;
            string_to_cell_of_cell(
                    tCMParameter.get< std::string >( "properties" ),
                    tPropertyNamesPair );
            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                // get the property name
                std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                // check if property in the map
                MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                        "Model_Initializer_Phasebased::create_CMs - Unknown aPropertyString : %s \n",
                        tPropertyName.c_str() );

                // get property index
                uint tPropertyIndex = mPropertyMap[ tPropertyName ];

                // set property for CM
                tCM->set_property(
                        mProperties( tPropertyIndex ),
                        tPropertyNamesPair( iProp )( 1 ) );
            }

            // set local properties
            tCM->set_local_properties();

            // check the phase exist
            MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                    "Model_Initializer_Phasebased::create_constitutive_models_without_phase - Unknown tPhaseName : %s \n",
                    tPhaseName.c_str() );

            // set CM to corresponding phase
            mPhaseInfo( mPhaseMap[ tPhaseName ] ).set_CM( tCM );
        }
    }

    void Model_Initializer_Phasebased::create_stabilization_parameters()
    {
        // create a stabilization parameter factory
        SP_Factory tSPFactory;

        // get the SP parameter list
        Vector< ParameterList > tSPParameterList = mParameterList( 2 );

        // get the number of stabilization parameters
        uint tNumSPs = tSPParameterList.size();

        // set size for the list of stabilization parameter pointer
        mStabilizationParameters.resize( tNumSPs, nullptr );

        // loop over the parameter list
        for ( uint iSP = 0; iSP < tNumSPs; iSP++ )
        {
            // get the stabilization parameters
            ParameterList tSPParameter = tSPParameterList( iSP );

            // get the stabilization parameter name
            std::string tSPName = tSPParameter.get< std::string >( "stabilization_name" );

            // get the stabilization type from parameter list
            fem::Stabilization_Type tSPType =
                    static_cast< fem::Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );

            // create a stabilization parameter pointer
            mStabilizationParameters( iSP ) = tSPFactory.create_SP( tSPType );

            // set name
            mStabilizationParameters( iSP )->set_name( tSPName );

            // set SP space dimension
            mStabilizationParameters( iSP )->set_space_dim( mSpatialDimension );

            // fill stabilization map
            mStabilizationParameterMap[ tSPName ] = iSP;

            // set parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
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
                if ( iLeader != 0u )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // get the treated IWG phase
                std::string tPhaseName =
                        tSPParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - Unknown phase name: %s \n",
                        tPhaseName.c_str() );

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
                mStabilizationParameters( iSP )->set_dof_type_list( tDofTypes, tDofTypeNames, tIsLeader );

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
                mStabilizationParameters( iSP )->set_dv_type_list( tDvTypes, tDvTypeNames, tIsLeader );

                // set leader properties
                Vector< Vector< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tSPParameter.get< std::string >( tIsLeaderString + "_properties" ),
                        tPropertyNamesPair );

                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get the property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check for unknown property
                    MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                            "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - Unknown %s aPropertyString : %s \n",
                            tIsLeaderString.c_str(),
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
                Vector< Vector< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tSPParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                // loop over CM names
                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // get the CM name
                    std::string tCMName = tCMNamesPair( iCM )( 0 );

                    // get CM from phase
                    std::shared_ptr< fem::Constitutive_Model > tCM =
                            mPhaseInfo( mPhaseMap[ tPhaseName ] ).get_CM_by_name( tCMName );

                    // set CM for SP
                    mStabilizationParameters( iSP )->set_constitutive_model(
                            tCM,
                            tCMNamesPair( iCM )( 1 ) );
                }

                // get the cluster measures specifications
                Vector< Vector< std::string > > tClusterMeasureTypes;
                string_to_cell_of_cell(
                        std::get< 0 >( tSPParameter.get< std::pair< std::string, std::string > >( "cluster_measures" ) ),
                        tClusterMeasureTypes );

                // get the cluster measures names
                Vector< std::string > tClusterMeasureNames;
                string_to_cell( std::get< 1 >( tSPParameter.get< std::pair< std::string, std::string > >( "cluster_measures" ) ),
                        tClusterMeasureNames );

                // build a cell of tuples describing the cluster measures specifications
                Vector< std::tuple<
                        fem::Measure_Type,
                        mtk::Primary_Void,
                        mtk::Leader_Follower > >
                        tClusterMeasureTuples( tClusterMeasureNames.size() );

                // get fem::Measure_Type, mtk::Primary_Void and mtk::Leader_Follower map
                // to convert string to enums
                moris::map< std::string, fem::Measure_Type >    tFemMeasureMap = fem::get_measure_type_map();
                moris::map< std::string, mtk::Primary_Void >    tMtkPrimaryMap = mtk::get_primary_type_map();
                moris::map< std::string, mtk::Leader_Follower > tMtkLeaderMap  = mtk::get_leader_type_map();

                // loop over cluster measures names
                for ( uint iCMEA = 0; iCMEA < tClusterMeasureNames.size(); iCMEA++ )
                {
                    // check that measure type is member of map
                    MORIS_ERROR( tFemMeasureMap.key_exists( tClusterMeasureTypes( iCMEA )( 0 ) ),
                            "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - key does not exist: %s",
                            tClusterMeasureTypes( iCMEA )( 0 ).c_str() );

                    // get fem measure type from map
                    fem::Measure_Type tFemMeasureType = tFemMeasureMap.find( tClusterMeasureTypes( iCMEA )( 0 ) );

                    // check that primary type is member of map
                    MORIS_ERROR( tMtkPrimaryMap.key_exists( tClusterMeasureTypes( iCMEA )( 1 ) ),
                            "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - key does not exist: %s",
                            tClusterMeasureTypes( iCMEA )( 1 ).c_str() );

                    // get mtk primary type from map
                    mtk::Primary_Void tMtkPrimaryType = tMtkPrimaryMap.find( tClusterMeasureTypes( iCMEA )( 1 ) );

                    // check that leader type is member of map
                    MORIS_ERROR( tMtkLeaderMap.key_exists( tClusterMeasureTypes( iCMEA )( 2 ) ),
                            "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - key does not exist: %s",
                            tClusterMeasureTypes( iCMEA )( 2 ).c_str() );

                    // get mtk leader type from map
                    mtk::Leader_Follower tMtkLeaderType = tMtkLeaderMap.find( tClusterMeasureTypes( iCMEA )( 2 ) );

                    // build the cluster measure specification tuple and set it in cell of tuples
                    tClusterMeasureTuples( iCMEA ) = std::make_tuple( tFemMeasureType, tMtkPrimaryType, tMtkLeaderType );
                }

                // set the cell of cluster measure specification tuples to the SP
                mStabilizationParameters( iSP )->set_cluster_measure_type_list(
                        tClusterMeasureTuples,
                        tClusterMeasureNames );
            }
        }
    }

    void Model_Initializer_Phasebased::create_iwgs()
    {
        // create an IWG factory
        IWG_Factory tIWGFactory;

        // get the IWG parameter list
        Vector< ParameterList > tIWGParameterList = mParameterList( 3 );

        // get number of IWGs
        uint tNumIWGs = tIWGParameterList.size();

        // create a list of IWG pointers
        mIWGs.resize( tNumIWGs, nullptr );

        // loop over the parameter lists
        for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
        {
            // get the treated IWG parameter list
            ParameterList tIWGParameter = tIWGParameterList( iIWG );

            // get the treated IWG name
            std::string tIWGName = tIWGParameter.get< std::string >( "IWG_name" );

            // get the treated IWG type
            fem::IWG_Type tIWGType =
                    static_cast< fem::IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );

            // get the ghost order from parameter list
            uint tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );

            // get the treated IWG residual dof type
            Vector< Vector< moris::MSI::Dof_Type > > tResDofTypes;
            string_to_cell_of_cell(
                    tIWGParameter.get< std::string >( "dof_residual" ),
                    tResDofTypes,
                    mMSIDofTypeMap );

            // get function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tIWGParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );

            // get the treated IWG bulk type
            fem::Element_Type tIWGBulkType =
                    static_cast< fem::Element_Type >( tIWGParameter.get< uint >( "IWG_bulk_type" ) );

            // set flag for leader/follower
            bool tLeaderFollower = ( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET ) || ( tIWGBulkType == fem::Element_Type::NONCONFORMAL_SIDESET );
            // create an IWG pointer
            mIWGs( iIWG ) = tIWGFactory.create_IWG( tIWGType );

            // set name
            mIWGs( iIWG )->set_name( tIWGName );

            // set interpolation order
            mIWGs( iIWG )->set_interpolation_order( tGhostOrder );

            // set residual dof type
            mIWGs( iIWG )->set_residual_dof_type( tResDofTypes );

            // set bulk type
            mIWGs( iIWG )->set_bulk_type( tIWGBulkType );

            // set constant parameters
            mIWGs( iIWG )->set_parameters( tFuncParameters );

            mIWGs( iIWG )->set_is_fd_jacobian( not tIWGParameter.get< bool >( "analytical_jacobian" ) );    // warning: this is a negation

            // initialize string for leader or follower
            std::string          tIsLeaderString = "leader";
            mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

            // loop on leader and follower
            for ( uint iLeader = 0; iLeader <= tLeaderFollower; iLeader++ )
            {
                // if follower
                if ( iLeader )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // get the treated IWG phase
                std::string tPhaseName =
                        tIWGParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "Unknown phase name: %s \n",
                        tPhaseName.c_str() );

                // set phase name
                mIWGs( iIWG )->set_phase_name( tPhaseName, tIsLeader );

                // get phase index
                uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                // get dof type list from phase - ignore double-sided side sets
                if ( !tLeaderFollower )
                {
                    mPhaseInfo( tPhaseIndex ).add_dof_type_to_list( tResDofTypes );
                }

                // set properties
                Vector< Vector< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( tIsLeaderString + "_properties" ),
                        tPropertyNamesPair );

                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check for unknown property
                    MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                            "Unknown %s aPropertyString: %s \n",
                            tIsLeaderString.c_str(),
                            tPropertyName.c_str() );

                    // get property index
                    uint tPropertyIndex = mPropertyMap[ tPropertyName ];

                    // set property for IWG
                    mIWGs( iIWG )->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ),
                            tIsLeader );
                }

                // set material model
                Vector< Vector< std::string > > tMMNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( tIsLeaderString + "_material_model" ),
                        tMMNamesPair );
                MORIS_ERROR( tMMNamesPair.size() <= 1, "Model_Initializer_Phasebased::create_iwgs() - Only one material model per CM allowed." );

                // loop over material model
                for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
                {
                    // get the MM name
                    std::string tMMName = tMMNamesPair( iMM )( 0 );

                    // get MM from phase
                    std::shared_ptr< fem::Material_Model > tMM =
                            mPhaseInfo( tPhaseIndex ).get_MM_by_name( tMMName );

                    // set CM for IWG
                    mIWGs( iIWG )->set_material_model(
                            tMM,
                            tMMNamesPair( iMM )( 1 ),
                            tIsLeader );
                }

                // set constitutive models
                Vector< Vector< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tIWGParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                // loop over constitutive models
                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // get the CM name
                    std::string tCMName = tCMNamesPair( iCM )( 0 );

                    // get CM from phase
                    std::shared_ptr< fem::Constitutive_Model > tCM =
                            mPhaseInfo( tPhaseIndex ).get_CM_by_name( tCMName );

                    // set CM for IWG
                    mIWGs( iIWG )->set_constitutive_model(
                            tCM,
                            tCMNamesPair( iCM )( 1 ),
                            tIsLeader );
                }
            }

            // set stabilization parameters
            Vector< Vector< std::string > > tSPNamesPair;
            string_to_cell_of_cell(
                    tIWGParameter.get< std::string >( "stabilization_parameters" ),
                    tSPNamesPair );

            // loop over SP names
            for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
            {
                // get the SP name
                std::string tSPName = tSPNamesPair( iSP )( 0 );

                // check for unknown SP
                MORIS_ERROR( mStabilizationParameterMap.find( tSPName ) != mStabilizationParameterMap.end(),
                        "Unknown aSPString: %s \n",
                        tSPName.c_str() );

                // get SP index
                uint tSPIndex = mStabilizationParameterMap[ tSPNamesPair( iSP )( 0 ) ];

                // set SP for IWG
                mIWGs( iIWG )->set_stabilization_parameter(
                        mStabilizationParameters( tSPIndex ),
                        tSPNamesPair( iSP )( 1 ) );
            }
        }

        // loop over the parameter lists to set dof dependencies
        for ( uint iIWG = 0; iIWG < tNumIWGs; iIWG++ )
        {
            // get the treated IWG parameter list
            ParameterList tIWGParameter = tIWGParameterList( iIWG );

            // get the IWG bulk type
            fem::Element_Type tIWGBulkType = mIWGs( iIWG )->get_bulk_type();

            // get the IWG leader phase name
            std::string tPhaseName =
                    mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::LEADER );

            // check for unknown phase
            MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                    "Unknown phase name: %s \n",
                    tPhaseName.c_str() );

            // get the phase index
            uint tPhaseIndex = mPhaseMap[ tPhaseName ];

            // get dof type list from leader phase
            const Vector< Vector< MSI::Dof_Type > > &tLeaderDofTypes =
                    mPhaseInfo( tPhaseIndex ).get_dof_type_list();

            // get dof type list from leader phase
            const Vector< Vector< gen::PDV_Type > > &tLeaderPdvTypes =
                    mPhaseInfo( tPhaseIndex ).get_dv_type_list();

            // set leader dof dependencies
            mIWGs( iIWG )->set_dof_type_list( tLeaderDofTypes, mtk::Leader_Follower::LEADER );

            // set leader dv dependencies
            mIWGs( iIWG )->set_dv_type_list( tLeaderPdvTypes, mtk::Leader_Follower::LEADER );

            if ( tIWGBulkType == fem::Element_Type::DOUBLE_SIDESET )
            {
                // get the IWG follower phase name
                std::string tFollowerPhaseName =
                        mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tFollowerPhaseName ) != mPhaseMap.end(),
                        "Unknown phase name: %s \n",
                        tFollowerPhaseName.c_str() );

                // get CM index
                uint tFollowerPhaseIndex = mPhaseMap[ tFollowerPhaseName ];

                // get dof type list from phase
                const Vector< Vector< MSI::Dof_Type > > &tFollowerDofTypes =
                        mPhaseInfo( tFollowerPhaseIndex ).get_dof_type_list();

                // get pdv type list from phase
                const Vector< Vector< gen::PDV_Type > > &tFollowerPdvTypes =
                        mPhaseInfo( tFollowerPhaseIndex ).get_dv_type_list();

                // set follower dof dependencies
                mIWGs( iIWG )->set_dof_type_list( tFollowerDofTypes, mtk::Leader_Follower::FOLLOWER );

                // set follower dv dependencies
                mIWGs( iIWG )->set_dv_type_list( tFollowerPdvTypes, mtk::Leader_Follower::FOLLOWER );
            }
        }
    }

    void Model_Initializer_Phasebased::create_iqis()
    {
        // create an IQI factory
        IQI_Factory tIQIFactory;

        // get the IQI parameter list
        Vector< ParameterList > tIQIParameterList = mParameterList( 4 );

        // get number of IQIs
        uint tNumIQIs = tIQIParameterList.size();

        // set size for list of IQI pointers
        mIQIs.resize( tNumIQIs, nullptr );

        // loop over the parameter lists
        for ( uint iIQI = 0; iIQI < tNumIQIs; iIQI++ )
        {
            // get the treated IQI parameter list
            ParameterList tIQIParameter = tIQIParameterList( iIQI );

            // get the treated IQI name from parameter list
            std::string tIQIName =
                    tIQIParameter.get< std::string >( "IQI_name" );

            // get the IQI type from parameter list
            fem::IQI_Type tIQIType =
                    static_cast< fem::IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );

            // get the quantity dof type from parameter list
            Vector< moris::MSI::Dof_Type > tQuantityDofTypes;
            string_to_cell(
                    tIQIParameter.get< std::string >( "dof_quantity" ),
                    tQuantityDofTypes,
                    mMSIDofTypeMap );

            // get the field index from parameter list
            sint tIQIFieldIndex =
                    tIQIParameter.get< moris::sint >( "vectorial_field_index" );

            // set function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters;
            string_to_cell_mat_2(
                    tIQIParameter.get< std::string >( "function_parameters" ),
                    tFuncParameters );

            // get the treated IQI bulk type
            fem::Element_Type tIQIBulkType =
                    static_cast< fem::Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );

            // set bool to true if double sideset
            bool tLeaderFollower = ( tIQIBulkType == fem::Element_Type::DOUBLE_SIDESET ) || ( tIQIBulkType == fem::Element_Type::NONCONFORMAL_SIDESET );

            // create an IQI pointer
            mIQIs( iIQI ) = tIQIFactory.create_IQI( tIQIType );

            // set name
            mIQIs( iIQI )->set_name( tIQIName );

            // set quantity dof type
            mIQIs( iIQI )->set_quantity_dof_type( tQuantityDofTypes );

            // set index for vectorial field
            mIQIs( iIQI )->set_output_type_index( tIQIFieldIndex );

            // set bulk type
            mIQIs( iIQI )->set_bulk_type( tIQIBulkType );

            // set constant parameters
            mIQIs( iIQI )->set_parameters( tFuncParameters );

            // init string for leader or follower
            std::string          tIsLeaderString = "leader";
            mtk::Leader_Follower tIsLeader       = mtk::Leader_Follower::LEADER;

            // loop on leader and follower
            for ( uint iLeader = 0; iLeader <= tLeaderFollower; iLeader++ )
            {
                // if follower
                if ( iLeader )
                {
                    // reset string for follower
                    tIsLeaderString = "follower";
                    tIsLeader       = mtk::Leader_Follower::FOLLOWER;
                }

                // get the treated IWG phase
                std::string tPhaseName = tIQIParameter.get< std::string >( tIsLeaderString + "_phase_name" );

                // check for unknown phase
                MORIS_ERROR( mPhaseMap.find( tPhaseName ) != mPhaseMap.end(),
                        "Model_Initializer_Phasebased::create_IQIs_without_phase - %s: Unknown %s phase name: %s.",
                        tIQIName.c_str(),
                        tIsLeaderString.c_str(),
                        tPhaseName.c_str() );

                // set phase name
                mIQIs( iIQI )->set_phase_name( tPhaseName, tIsLeader );

                // get the phase index
                uint tPhaseIndex = mPhaseMap[ tPhaseName ];

                // get dof type list from phase
                const Vector< Vector< MSI::Dof_Type > > &tDofTypes =
                        mPhaseInfo( tPhaseIndex ).get_dof_type_list();

                // get dof type list from phase
                const Vector< Vector< gen::PDV_Type > > &tDvTypes =
                        mPhaseInfo( tPhaseIndex ).get_dv_type_list();

                // set leader dof dependencies
                mIQIs( iIQI )->set_dof_type_list( tDofTypes );

                // set leader dv dependencies
                mIQIs( iIQI )->set_dv_type_list( tDvTypes );

                // set leader properties
                Vector< Vector< std::string > > tPropertyNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_properties" ),
                        tPropertyNamesPair );

                for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
                {
                    // get the property name
                    std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );

                    // check for unknown property
                    MORIS_ERROR( mPropertyMap.find( tPropertyName ) != mPropertyMap.end(),
                            "Model_Initializer_Phasebased::create_IQIs_without_phase - Unknown %s aPropertyString: %s \n",
                            tIsLeaderString.c_str(),
                            tPropertyName.c_str() );

                    // get property index
                    uint tPropertyIndex = mPropertyMap[ tPropertyName ];

                    // set property for IWG
                    mIQIs( iIQI )->set_property(
                            mProperties( tPropertyIndex ),
                            tPropertyNamesPair( iProp )( 1 ),
                            tIsLeader );
                }

                // set leader constitutive models
                Vector< Vector< std::string > > tCMNamesPair;
                string_to_cell_of_cell(
                        tIQIParameter.get< std::string >( tIsLeaderString + "_constitutive_models" ),
                        tCMNamesPair );

                for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
                {
                    // get CM name
                    std::string tCMName = tCMNamesPair( iCM )( 0 );

                    // get CM from phase
                    std::shared_ptr< fem::Constitutive_Model > tCM =
                            mPhaseInfo( tPhaseIndex ).get_CM_by_name( tCMName );

                    // set CM for IQI
                    mIQIs( iIQI )->set_constitutive_model(
                            tCM,
                            tCMNamesPair( iCM )( 1 ),
                            tIsLeader );
                }
            }

            // set stabilization parameters
            Vector< Vector< std::string > > tSPNamesPair;
            string_to_cell_of_cell(
                    tIQIParameter.get< std::string >( "stabilization_parameters" ),
                    tSPNamesPair );

            for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
            {
                // get the SP name
                std::string tSPName = tSPNamesPair( iSP )( 0 );

                // check for unknown SP
                MORIS_ERROR( mStabilizationParameterMap.find( tSPName ) != mStabilizationParameterMap.end(),
                        "Model_Initializer_Phasebased::create_IQIs_without_phase - Unknown aSPString: %s \n",
                        tSPName.c_str() );

                // get SP index
                uint tSPIndex = mStabilizationParameterMap[ tSPName ];

                // set SP for IWG
                mIQIs( iIQI )->set_stabilization_parameter(
                        mStabilizationParameters( tSPIndex ),
                        tSPNamesPair( iSP )( 1 ) );
            }
        }
    }

    void Model_Initializer_Phasebased::create_set_info()
    {
        // init number of fem sets to be created
        uint tNumFEMSets = 0;

        // get the IWG and IQI parameter lists
        Vector< ParameterList > tIWGParameterList = mParameterList( 3 );
        Vector< ParameterList > tIQIParameterList = mParameterList( 4 );

        // get fem computation type parameter list
        ParameterList tComputationParameterList = mParameterList( 5 )( 0 );

        // bool true for analytical forward analysis, false for finite difference
        // decide if dRdu and dQIdu are computed by A/FD
        bool tIsAnalyticalFA =
                tComputationParameterList.get< bool >( "is_analytical_forward" );

        // get enum for FD scheme for forward analysis
        fem::FDScheme_Type tFDSchemeForFA = static_cast< fem::FDScheme_Type >(
                tComputationParameterList.get< uint >( "finite_difference_scheme_forward" ) );

        // get perturbation size for FD for forward analysis
        real tFDPerturbationFA = tComputationParameterList.get< real >(
                "finite_difference_perturbation_size_forward" );

        // get bool for analytical/finite difference for sensitivity analysis
        // decide if dRdp and dQIdp are computed by A/FD
        bool tIsAnalyticalSA =
                tComputationParameterList.get< bool >( "is_analytical_sensitivity" );

        // get enum for FD scheme for sensitivity analysis
        fem::FDScheme_Type tFDSchemeForSA = static_cast< fem::FDScheme_Type >(
                tComputationParameterList.get< uint >( "finite_difference_scheme" ) );

        // get perturbation size for FD for sensitivity analysis
        real tFDPerturbationSA = tComputationParameterList.get< real >(
                "finite_difference_perturbation_size" );

        // get enum for perturbation strategy for finite difference
        fem::Perturbation_Type tPerturbationStrategy = static_cast< fem::Perturbation_Type >(
                tComputationParameterList.get< uint >( "finite_difference_perturbation_strategy" ) );

        mtk::Integration_Order tIntegrationOrder = static_cast< mtk::Integration_Order >( tComputationParameterList.get< uint >( "nonconformal_integration_order" ) );
        real tMaxNegativeRayLength = tComputationParameterList.get< real >( "nonconformal_max_negative_ray_length" );

        // create a map of the set
        std::map< std::tuple< std::string, bool, bool >, uint > tMeshToFemSet;

        // loop over the IWGs
        for ( uint iIWG = 0; iIWG < tIWGParameterList.size(); iIWG++ )
        {
            // get the treated IWG parameter list
            ParameterList tIWGParameter = tIWGParameterList( iIWG );

            // get the IWG bulk type
            fem::Element_Type tIWGBulkType = mIWGs( iIWG )->get_bulk_type();

            // get time continuity flag
            bool tTimeContinuity = mIWGs( iIWG )->get_time_continuity();

            // get time boundary flag
            bool tTimeBoundary = mIWGs( iIWG )->get_time_boundary();

            // get bool for ghost
            bool tIsGhost = mIWGs( iIWG )->get_ghost_flag();

            // get the IWG leader phase name
            std::string tLeaderPhaseName = mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::LEADER );

            // get the IWG follower phase name
            std::string tFollowerPhaseName = mIWGs( iIWG )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

            // get follower phase string from IWG input
            std::string tNeighborPhasesName = tIWGParameter.get< std::string >( "neighbor_phases" );

            // get ordinal string from IWG input
            std::string tOrdinalString = tIWGParameter.get< std::string >( "side_ordinals" );

            // get mesh set names for IWG
            Vector< std::string > tMeshSetNames;
            this->get_mesh_set_names(
                    tIWGBulkType,
                    tLeaderPhaseName,
                    tFollowerPhaseName,
                    tNeighborPhasesName,
                    tOrdinalString,
                    tIsGhost,
                    tMeshSetNames );

            // get a representative DoF type
            MSI::Dof_Type tFirstResidualDofType = mIWGs( iIWG )->get_residual_dof_type()( 0 )( 0 );

            // loop over the mesh set names
            for ( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
            {
                // get the name of the set currently treated
                std::string tMeshSetName = tMeshSetNames( iSetName );

                // check for ghost set names and select correct B-spline mesh automatically when new ghost sets need to be used
                this->check_and_set_ghost_set_names( tMeshSetName, tFirstResidualDofType );

                // check if the mesh set name already in map
                if ( tMeshToFemSet.find( std::make_tuple(
                             tMeshSetName,
                             tTimeContinuity,
                             tTimeBoundary ) )
                        == tMeshToFemSet.end() )
                {
                    // add the mesh set name map
                    tMeshToFemSet[ std::make_tuple(
                            tMeshSetName,
                            tTimeContinuity,
                            tTimeBoundary ) ] = tNumFEMSets++;

                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;

                    // set its mesh set name
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );

                    // set its time continuity flag
                    aSetUserInfo.set_time_continuity( tTimeContinuity );

                    // set its time boundary flag
                    aSetUserInfo.set_time_boundary( tTimeBoundary );

                    // set its forward analysis type flag
                    aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );

                    // set its FD scheme for forward analysis
                    aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );

                    // set its FD perturbation size for forward analysis
                    aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                    // set its sensitivity analysis type flag
                    aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                    // set its FD scheme for sensitivity analysis
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                    // set its FD perturbation size for sensitivity analysis
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );

                    // set its perturbation strategy for finite difference
                    aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                    // set the integration order for nonconformal elements
                    aSetUserInfo.set_integration_order( tIntegrationOrder );

                    // set the maximum negative ray length for nonconformal elements
                    aSetUserInfo.set_max_negative_ray_length( tMaxNegativeRayLength );

                    // set the IWG
                    aSetUserInfo.add_IWG( mIWGs( iIWG ) );

                    // add it to the list of fem set info
                    mSetInfo.push_back( aSetUserInfo );
                }
                else
                {
                    // set the IWG
                    mSetInfo( tMeshToFemSet[ std::make_tuple(
                                      tMeshSetName,
                                      tTimeContinuity,
                                      tTimeBoundary ) ] )
                            .add_IWG( mIWGs( iIWG ) );
                }
            }
        }

        // loop over the IQIs
        for ( uint iIQI = 0; iIQI < tIQIParameterList.size(); iIQI++ )
        {
            // get the treated IWG parameter list
            ParameterList tIQIParameter = tIQIParameterList( iIQI );

            // get the IWG bulk type
            fem::Element_Type tIQIBulkType = mIQIs( iIQI )->get_bulk_type();

            // get time continuity flag
            bool tTimeContinuity = mIQIs( iIQI )->get_time_continuity();

            // get time boundary flag
            bool tTimeBoundary = mIQIs( iIQI )->get_time_boundary();

            // get the IWG leader phase name
            std::string tLeaderPhaseName =
                    mIQIs( iIQI )->get_phase_name( mtk::Leader_Follower::LEADER );

            // get the IWG follower phase name
            std::string tFollowerPhaseName =
                    mIQIs( iIQI )->get_phase_name( mtk::Leader_Follower::FOLLOWER );

            // get follower phase string from IQI input
            std::string tFollowerPhaseString =
                    tIQIParameter.get< std::string >( "neighbor_phases" );

            // get ordinal string from IQI input
            std::string tOrdinalString =
                    tIQIParameter.get< std::string >( "side_ordinals" );

            // get mesh set names for IWG
            Vector< std::string > tMeshSetNames;
            this->get_mesh_set_names(
                    tIQIBulkType,
                    tLeaderPhaseName,
                    tFollowerPhaseName,
                    tFollowerPhaseString,
                    tOrdinalString,
                    false,
                    tMeshSetNames );

            // loop over the mesh set names
            for ( uint iSetName = 0; iSetName < tMeshSetNames.size(); iSetName++ )
            {
                // get the name of the set currently treated
                std::string tMeshSetName = tMeshSetNames( iSetName );

                // if the mesh set name not in map
                if ( tMeshToFemSet.find( std::make_tuple(
                             tMeshSetName,
                             tTimeContinuity,
                             tTimeBoundary ) )
                        == tMeshToFemSet.end() )
                {
                    // add the mesh set name map
                    tMeshToFemSet[ std::make_tuple(
                            tMeshSetName,
                            tTimeContinuity,
                            tTimeBoundary ) ] = tNumFEMSets++;

                    // create a fem set info for the mesh set
                    Set_User_Info aSetUserInfo;

                    // set its mesh set name
                    aSetUserInfo.set_mesh_set_name( tMeshSetName );

                    // set its time continuity flag
                    aSetUserInfo.set_time_continuity( tTimeContinuity );

                    // set its time boundary flag
                    aSetUserInfo.set_time_boundary( tTimeBoundary );

                    // set its forward analysis type flag
                    aSetUserInfo.set_is_analytical_forward_analysis( tIsAnalyticalFA );

                    // set its FD scheme for forward analysis
                    aSetUserInfo.set_finite_difference_scheme_for_forward_analysis( tFDSchemeForFA );

                    // set its FD perturbation size for forward analysis
                    aSetUserInfo.set_finite_difference_perturbation_size_for_forward_analysis( tFDPerturbationFA );

                    // set its sensitivity analysis type flag
                    aSetUserInfo.set_is_analytical_sensitivity_analysis( tIsAnalyticalSA );

                    // set its FD scheme for sensitivity analysis
                    aSetUserInfo.set_finite_difference_scheme_for_sensitivity_analysis( tFDSchemeForSA );

                    // set its FD perturbation size for sensitivity analysis
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );

                    // set its perturbation strategy for finite difference
                    aSetUserInfo.set_perturbation_strategy( tPerturbationStrategy );

                    // set the IQI
                    aSetUserInfo.add_IQI( mIQIs( iIQI ) );

                    // add it to the list of fem set info
                    mSetInfo.push_back( aSetUserInfo );
                }
                else
                {
                    // set the IQI
                    mSetInfo( tMeshToFemSet[ std::make_tuple(
                                      tMeshSetName,
                                      tTimeContinuity,
                                      tTimeBoundary ) ] )
                            .add_IQI( mIQIs( iIQI ) );
                }
            }
        }
    }

    void
    Model_Initializer_Phasebased::get_mesh_set_names(
            fem::Element_Type      aIWGBulkType,
            const std::string     &aLeaderPhaseName,
            const std::string     &aFollowerPhaseName,
            const std::string     &aNeighborPhaseString,
            const std::string     &aOrdinalString,
            bool                   aIsGhost,
            Vector< std::string > &aMeshSetNames )
    {
        // get the leader phase mesh index
        MORIS_ERROR( mPhaseMap.find( aLeaderPhaseName ) != mPhaseMap.end(), "Unknown leader phase name: %s \n", aLeaderPhaseName.c_str() );
        moris::Matrix< moris::IndexMat > tLeaderPhaseIndices = mPhaseInfo( mPhaseMap[ aLeaderPhaseName ] ).get_phase_indices();
        uint                             tNumLeaderIndices   = tLeaderPhaseIndices.numel();

        // switch on the element type
        switch ( aIWGBulkType )
        {
            case fem::Element_Type::BULK:
            {
                // loop over phase mesh indices
                for ( uint iMeshIndex = 0; iMeshIndex < tNumLeaderIndices; iMeshIndex++ )
                {
                    //                        // FIXME ! get mesh set names from integration mesh for index
                    //                        mMeshManager->get_integration_mesh( 0 )->
                    //                                get_block_set_names_with_color( tLeaderPhaseIndices( iMeshIndex ), aMeshSetNames );

                    // add mesh set name to list
                    aMeshSetNames.push_back(
                            "HMR_dummy_c_p" + std::to_string( tLeaderPhaseIndices( iMeshIndex ) ) );

                    // add mesh set name to list
                    aMeshSetNames.push_back(
                            "HMR_dummy_n_p" + std::to_string( tLeaderPhaseIndices( iMeshIndex ) ) );
                }
                break;
            }
            case fem::Element_Type::SIDESET:
            {
                // get neighbor phase names from string
                Vector< std::string > tNeighborPhaseNames = string_to_cell< std::string >( aNeighborPhaseString );

                // get number of neighbor phase
                uint tNumSingle = tNeighborPhaseNames.size();

                // get ordinals for boundary from string
                Matrix< DDSMat > tOrdinals    = string_to_mat< DDSMat >( aOrdinalString );
                uint             tNumBoundary = tOrdinals.numel();

                // loop over leader phase mesh indices
                for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                {
                    // get single sideset
                    for ( uint iSingle = 0; iSingle < tNumSingle; iSingle++ )
                    {
                        // get the neighbor phase name
                        std::string tNeighborPhaseName = tNeighborPhaseNames( iSingle );

                        // get the follower phase mesh index
                        moris::Matrix< moris::IndexMat > tNeighborPhaseIndices = mPhaseInfo( mPhaseMap[ tNeighborPhaseName ] ).get_phase_indices();

                        // get number of neighbor phase mesh indices
                        uint tNumNeighborIndices = tNeighborPhaseIndices.numel();

                        for ( uint iNeighborMeshIndex = 0; iNeighborMeshIndex < tNumNeighborIndices; iNeighborMeshIndex++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "iside_b0_" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) + "_b1_" + std::to_string( tNeighborPhaseIndices( iNeighborMeshIndex ) ) );
                        }
                    }

                    // get boundary sideset
                    for ( uint iBoundary = 0; iBoundary < tNumBoundary; iBoundary++ )
                    {
                        // FIXME get this info from the mesh
                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "SideSet_" + std::to_string( tOrdinals( iBoundary ) ) + "_c_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );

                        // FIXME get this info from the mesh
                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "SideSet_" + std::to_string( tOrdinals( iBoundary ) ) + "_n_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );
                    }
                }
                break;
            }
            case fem::Element_Type::DOUBLE_SIDESET:
            {
                // if ghost
                if ( aIsGhost )
                {
                    // loop over leader phase mesh indices
                    for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                    {
                        // FIXME get this info from the mesh
                        // add mesh set name to list
                        aMeshSetNames.push_back(
                                "ghost_p" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) );
                    }
                }
                // if interface
                else
                {
                    MORIS_ERROR( aLeaderPhaseName != aFollowerPhaseName,
                            "Model_Initializer_Phasebased::get_mesh_set_names - Leader and follower phases are the same, FIXME case not handled yet " );

                    // get the follower phase mesh index
                    moris::Matrix< moris::IndexMat > tFollowerPhaseIndices =
                            mPhaseInfo( mPhaseMap[ aFollowerPhaseName ] ).get_phase_indices();

                    // get number of follower phase mesh index
                    uint tNumFollowerIndices = tFollowerPhaseIndices.numel();

                    // loop over leader phase mesh indices
                    for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                    {
                        // get leader index
                        uint tLeaderPhaseIndex = tLeaderPhaseIndices( iLeaderMeshIndex );

                        // loop over follower phase mesh indices
                        for ( uint iFollowerMeshIndex = 0; iFollowerMeshIndex < tNumFollowerIndices; iFollowerMeshIndex++ )
                        {
                            // get follower index
                            uint tFollowerPhaseIndex = tFollowerPhaseIndices( iFollowerMeshIndex );

                            // if leader and follower index are different
                            if ( tLeaderPhaseIndex != tFollowerPhaseIndex )
                            {
                                // FIXME get this info from the mesh
                                // get interface name
                                aMeshSetNames.push_back(
                                        "dbl_iside_p0_" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) + "_p1_" + std::to_string( tFollowerPhaseIndices( iFollowerMeshIndex ) ) );
                            }
                        }
                    }
                }
                break;
            }
            case Element_Type::NONCONFORMAL_SIDESET:
            {
                MORIS_ERROR( !aIsGhost, "Nonconformal sideset cannot be a ghost set" );
                MORIS_ERROR( !aNeighborPhaseString.empty() && !aFollowerPhaseName.empty(), "Nonconformal sideset requires leader, follower and neighbor (void) phase names" );
                MORIS_ERROR( aLeaderPhaseName != aFollowerPhaseName, "Leader and follower phases are the same, FIXME case not handled yet " );

                Vector< std::string > tNeighborPhaseNames = string_to_cell< std::string >( aNeighborPhaseString );
                Vector< uint >        tNeighborPhaseIndices;
                for ( auto const &tNeighborPhaseName : tNeighborPhaseNames )
                {
                    MORIS_ERROR( tNeighborPhaseName != aLeaderPhaseName, "Nonconformal sideset cannot have leader phase as neighbor" );
                    MORIS_ERROR( tNeighborPhaseName != aFollowerPhaseName, "Nonconformal sideset cannot have follower phase as neighbor" );
                    MORIS_ERROR( mPhaseMap.find( tNeighborPhaseName ) != mPhaseMap.end(), "Unknown neighbor phase name: %s", tNeighborPhaseName.c_str() );

                    // append the neighbor indices into the continuous list
                    Matrix< IndexMat > tNeighborIndices = mPhaseInfo( mPhaseMap[ tNeighborPhaseName ] ).get_phase_indices();
                    for ( auto const tNeighborIndex : tNeighborIndices )
                    {
                        tNeighborPhaseIndices.push_back( tNeighborIndex );
                    }
                }

                MORIS_ERROR( mPhaseMap.find( aFollowerPhaseName ) != mPhaseMap.end(), "Unknown follower phase name: %s", aFollowerPhaseName.c_str() );
                Matrix< moris::IndexMat > tFollowerPhaseIndices = mPhaseInfo( mPhaseMap[ aFollowerPhaseName ] ).get_phase_indices();

                for ( auto const tLeaderPhaseIndex : tLeaderPhaseIndices )
                {
                    for ( auto const tFollowerPhaseIndex : tFollowerPhaseIndices )
                    {
                        for ( auto const tNeighborPhaseIndex : tNeighborPhaseIndices )
                        {
                            std::string tMeshName = "ncss";    // TODO @ff: this is a temporary solution, we need to find a better way to name the nonconformal sidesets
                            tMeshName += "|iside_b0_" + std::to_string( tLeaderPhaseIndex ) + "_b1_" + std::to_string( tNeighborPhaseIndex );
                            tMeshName += "|iside_b0_" + std::to_string( tFollowerPhaseIndex ) + "_b1_" + std::to_string( tNeighborPhaseIndex );
                            aMeshSetNames.push_back( tMeshName );
                        }
                    }
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Model_Initializer_Phasebased::get_mesh_set_names - Unknown set type" );
            }
        }
    }

    void Model_Initializer_Phasebased::print_physics_model()
    {
        ParameterList tComputationParameterList = this->mParameterList( 5 )( 0 );
        bool          tPrintPhysics             = tComputationParameterList.get< bool >( "print_physics_model" );
        if ( tPrintPhysics && par_rank() == 0 )
        {
            std::cout << "Phase info \n";
            for ( auto &tPhaseInfo : mPhaseInfo )
            {
                std::cout << "%-------------------------------------------------\n";
                tPhaseInfo.print_names();
                std::cout << "%-------------------------------------------------\n";
            }
            std::cout << " \n";
        }
        Model_Initializer::print_physics_model();
    }

}    // namespace moris::fem
