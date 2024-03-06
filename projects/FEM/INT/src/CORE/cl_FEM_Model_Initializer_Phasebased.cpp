//
// Created by frank on 12/14/23.
//

#include "cl_FEM_Model_Initializer_Phasebased.hpp"
#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_SP_Factory.hpp"
#include "cl_FEM_MM_Factory.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "cl_FEM_Set_User_Info.hpp"
#include "cl_Map.hpp"
#include "cl_Param_List.hpp"
#include "cl_Vector.hpp"
#include "fn_Parsing_Tools.hpp"

#include <cl_FEM_Cluster_Measure.hpp>
#include <tuple>
#include <set>

namespace moris::fem
{
    void Model_Initializer_Phasebased::initialize()
    {
        this->create_phases();
        Model_Initializer::initialize();
    }

    void Model_Initializer_Phasebased::create_phases()
    {
        Vector< ParameterList > tPhaseParameterList = mParameterList( FEM_PARAMETER_PHASE_INDEX );
        for ( auto const &tPhaseParameter : tPhaseParameterList )
        {
            auto const tPhaseInfo = std::make_shared< Phase_User_Info >();
            auto const tPhaseName = tPhaseParameter.get< std::string >( "phase_name" );
            tPhaseInfo->set_phase_name( tPhaseName );
            tPhaseInfo->set_phase_indices( string_to_mat< IndexMat >( tPhaseParameter.get< std::string >( "phase_indices" ) ) );
            mPhases[ tPhaseName ] = tPhaseInfo;
        }
    }

    void Model_Initializer_Phasebased::create_material_models()
    {
        // create a constitutive model factory
        MM_Factory              tMMFactory;
        Vector< ParameterList > tMMParameterList = mParameterList( FEM_PARAMETER_MM_INDEX );

        for ( auto const &tMMParameter : tMMParameterList )
        {
            // create a material model pointer
            auto const                              tMMType = static_cast< Material_Type >( tMMParameter.get< uint >( "material_type" ) );
            std::shared_ptr< Material_Model > const tMM     = tMMFactory.create_MM( tMMType );
            tMM->set_space_dim( mSpatialDimension );

            std::string const tMMName = tMMParameter.get< std::string >( "material_name" );
            tMM->set_name( tMMName );

            // set MM dof dependencies
            auto const &[ tDofTypes, tDofTypeNames ] = read_dof_dependencies( tMMParameter, mtk::Leader_Follower::UNDEFINED );
            tMM->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set MM properties
            Vector< Vector< std::string > > tPropertyNamesPair = string_to_cell_of_cell< std::string >( tMMParameter.get< std::string >( "properties" ) );
            for ( uint iProp = 0; iProp < tPropertyNamesPair.size(); iProp++ )
            {
                std::string tPropertyName = tPropertyNamesPair( iProp )( 0 );    // get the property name
                // check if property in the map
                ensure_existing_parameter( mProperties, tPropertyName, "properties" );
                tMM->set_property( mProperties[ tPropertyName ], tPropertyNamesPair( iProp )( 1 ) );
            }
            // set local properties
            tMM->set_local_properties();

            auto tPhaseName = tMMParameter.get< std::string >( "phase_name" );
            ensure_existing_parameter( mPhases, tPhaseName, "phase_name" );
            mPhases[ tPhaseName ]->set_MM( tMM );
        }
    }

    void Model_Initializer_Phasebased::create_constitutive_models()
    {
        // create a constitutive model factory
        CM_Factory tCMFactory;

        // get the CM parameter list
        Vector< ParameterList > tCMParameterList = mParameterList( FEM_PARAMETER_CM_INDEX );
        for ( auto const &tCMParameter : tCMParameterList )
        {
            // get the constitutive type from parameter list
            auto const                                  tCMType = static_cast< Constitutive_Type >( tCMParameter.get< uint >( "constitutive_type" ) );
            std::shared_ptr< Constitutive_Model > const tCM     = tCMFactory.create_CM( tCMType );

            // get the constitutive model name from parameter list
            std::string const tCMName = tCMParameter.get< std::string >( "constitutive_name" );
            tCM->set_name( tCMName );

            // set CM model type. must come before "set_space_dim"
            // fixme: currently cannot set a plane type and tensor type at the same time from an input file
            auto const tCMModelType = static_cast< Model_Type >( tCMParameter.get< uint >( "model_type" ) );
            if ( tCMModelType != Model_Type::UNDEFINED )
            {
                tCM->set_model_type( tCMModelType );
            }

            // set CM space dimension
            tCM->set_space_dim( mSpatialDimension );

            // set function parameters
            Vector< moris::Matrix< DDRMat > > const tFuncParameters = string_to_cell_mat_2< DDRMat >( tCMParameter.get< std::string >( "function_parameters" ) );
            tCM->set_parameters( tFuncParameters );

            // set CM dof dependencies
            auto const &[ tDofTypes, tDofTypeNames ] = read_dof_dependencies( tCMParameter, mtk::Leader_Follower::UNDEFINED );
            tCM->set_dof_type_list( tDofTypes, tDofTypeNames );

            // set CM dv dependencies
            auto const &[ tDvTypes, tDvTypeNames ] = read_dv_dependencies( tCMParameter, mtk::Leader_Follower::UNDEFINED );
            tCM->set_dv_type_list( tDvTypes, tDvTypeNames );

            // set CM material model
            Vector< Vector< std::string > > tMMNamesPair = string_to_cell_of_cell< std::string >( tCMParameter.get< std::string >( "material_model" ) );
            MORIS_ERROR( tMMNamesPair.size() <= 1, "Model_Initializer_Phasebased::create_CMs() - Only one material model per CM allowed." );

            // get the phase from parameter list
            std::string const tPhaseName = read_phase_name( tCMParameter, mtk::Leader_Follower::UNDEFINED );
            mPhases[ tPhaseName ]->set_CM( tCM );

            // loop over Material Model names
            for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
            {
                // get the material name
                std::string const                       tMaterialName = tMMNamesPair( iMM )( 0 );
                std::shared_ptr< Material_Model > const tMM           = mPhases[ tPhaseName ]->get_MM_by_name( tMaterialName );
                tCM->set_material_model( tMM, tMMNamesPair( iMM )( 1 ) );
            }

            // set CM properties
            for ( auto const &[ tProperty, tPropertyString ] : read_properties( tCMParameter, mtk::Leader_Follower::UNDEFINED ) )
            {
                tCM->set_property( tProperty, tPropertyString );
            }
            // set local properties
            tCM->set_local_properties();
        }
    }

    void Model_Initializer_Phasebased::create_stabilization_parameters()
    {
        // create a stabilization parameter factory
        SP_Factory tSPFactory;

        // get the SP parameter list
        Vector< ParameterList > tSPParameterList = mParameterList( FEM_PARAMETER_SP_INDEX );

        // loop over the parameter list
        for ( auto const &tSPParameter : tSPParameterList )
        {
            // get the stabilization type from parameter list
            auto const                                       tSPType = static_cast< Stabilization_Type >( tSPParameter.get< uint >( "stabilization_type" ) );
            std::shared_ptr< Stabilization_Parameter > const tSP     = tSPFactory.create_SP( tSPType );

            // set name
            std::string const tSPName = tSPParameter.get< std::string >( "stabilization_name" );
            tSP->set_name( tSPName );
            // set SP space dimension
            tSP->set_space_dim( mSpatialDimension );

            // set parameters
            Vector< moris::Matrix< DDRMat > > const tFuncParameters = string_to_cell_mat_2< DDRMat >( tSPParameter.get< std::string >( "function_parameters" ) );
            tSP->set_parameters( tFuncParameters );

            // loop on leader and follower
            Vector< mtk::Leader_Follower > tLeaderFollowerParameters = { mtk::Leader_Follower::LEADER };
            if ( tSP->get_has_follower() )
            {
                tLeaderFollowerParameters.push_back( mtk::Leader_Follower::FOLLOWER );
            }

            for ( auto const &tLeaderFollower : tLeaderFollowerParameters )
            {
                // set dof dependencies
                auto [ tDofTypes, tDofTypeNames ] = read_dof_dependencies( tSPParameter, tLeaderFollower );
                tSP->set_dof_type_list( tDofTypes, tDofTypeNames, tLeaderFollower );

                auto [ tDvTypes, tDvTypeNames ] = read_dv_dependencies( tSPParameter, tLeaderFollower );
                tSP->set_dv_type_list( tDvTypes, tDvTypeNames, tLeaderFollower );

                // set properties
                for ( auto const &[ tProperty, tPropertyString ] : read_properties( tSPParameter, tLeaderFollower ) )
                {
                    tSP->set_property( tProperty, tPropertyString, tLeaderFollower );
                }

                for ( auto const &[ tConstitutiveModel, tConstitutiveModelString ] : read_constitutive_models( tSPParameter, tLeaderFollower ) )
                {
                    tSP->set_constitutive_model( tConstitutiveModel, tConstitutiveModelString, tLeaderFollower );
                }
            }

            // get the cluster measures specifications
            auto const                      tClusterMeasuresParameter = tSPParameter.get< std::pair< std::string, std::string > >( "cluster_measures" );
            Vector< Vector< std::string > > tClusterMeasureTypes      = string_to_cell_of_cell< std::string >( tClusterMeasuresParameter.first );
            Vector< std::string >           tClusterMeasureNames      = string_to_cell< std::string >( tClusterMeasuresParameter.second );

            // build a cell of tuples describing the cluster measures specifications
            Vector< Cluster_Measure::ClusterMeasureSpecification > tClusterMeasureSpecifications;

            // get Measure_Type, mtk::Primary_Void and mtk::Leader_Follower map to convert string to enums
            moris::map< std::string, Measure_Type >         tFemMeasureMap = get_measure_type_map();
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
                Measure_Type const tFemMeasureType = tFemMeasureMap.find( tClusterMeasureTypes( iCMEA )( 0 ) );

                // check that primary type is member of map
                MORIS_ERROR( tMtkPrimaryMap.key_exists( tClusterMeasureTypes( iCMEA )( 1 ) ),
                        "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - key does not exist: %s",
                        tClusterMeasureTypes( iCMEA )( 1 ).c_str() );

                // get mtk primary type from map
                mtk::Primary_Void const tMtkPrimaryType = tMtkPrimaryMap.find( tClusterMeasureTypes( iCMEA )( 1 ) );

                // check that leader type is member of map
                MORIS_ERROR( tMtkLeaderMap.key_exists( tClusterMeasureTypes( iCMEA )( 2 ) ),
                        "Model_Initializer_Phasebased::create_stabilization_parameters_without_phase - key does not exist: %s",
                        tClusterMeasureTypes( iCMEA )( 2 ).c_str() );

                // get mtk leader type from map
                mtk::Leader_Follower const tMtkLeaderType = tMtkLeaderMap.find( tClusterMeasureTypes( iCMEA )( 2 ) );

                // build the cluster measure specification tuple and set it in cell of tuples
                tClusterMeasureSpecifications.emplace_back( tFemMeasureType, tMtkPrimaryType, tMtkLeaderType );
            }
            tSP->set_cluster_measure_type_list( tClusterMeasureSpecifications, tClusterMeasureNames );

            mStabilizationParameters[ tSPName ] = tSP;
        }
    }

    void Model_Initializer_Phasebased::create_iwgs()
    {
        // create an IWG factory
        IWG_Factory             tIWGFactory;
        Vector< ParameterList > tIWGParameterList = mParameterList( FEM_PARAMETER_IWG_INDEX );

        // loop over the parameter lists
        for ( auto const &tIWGParameter : tIWGParameterList )
        {
            auto const                   tIWGType = static_cast< IWG_Type >( tIWGParameter.get< uint >( "IWG_type" ) );
            std::shared_ptr< IWG > const tIWG     = tIWGFactory.create_IWG( tIWGType );

            // get the treated IWG name
            std::string const tIWGName = tIWGParameter.get< std::string >( "IWG_name" );
            tIWG->set_name( tIWGName );
            mIWGs[ tIWGName ] = tIWG;

            // get the ghost order from parameter list
            uint const tGhostOrder = tIWGParameter.get< uint >( "ghost_order" );
            tIWG->set_interpolation_order( tGhostOrder );

            // get the treated IWG residual dof type
            Vector< Vector< moris::MSI::Dof_Type > > tResDofTypes = string_to_cell_of_cell( tIWGParameter.get< std::string >( "dof_residual" ), mMSIDofTypeMap );
            tIWG->set_residual_dof_type( tResDofTypes );

            // get function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters = string_to_cell_mat_2< DDRMat >( tIWGParameter.get< std::string >( "function_parameters" ) );
            tIWG->set_parameters( tFuncParameters );

            // get the treated IWG bulk type
            auto const tIWGBulkType = static_cast< Element_Type >( tIWGParameter.get< uint >( "IWG_bulk_type" ) );
            tIWG->set_bulk_type( tIWGBulkType );

            // in the case of bulk or side sets, the following loop has to be done only for the leader side.
            // if the element type is a double side set or nonconformal side set, the loop is done for both leader and follower side
            Vector< mtk::Leader_Follower > tLeaderFollowerSwitch = { mtk::Leader_Follower::LEADER };
            if ( ( tIWGBulkType == Element_Type::DOUBLE_SIDESET ) || ( tIWGBulkType == Element_Type::NONCONFORMAL_SIDESET ) )
            {
                tLeaderFollowerSwitch.push_back( mtk::Leader_Follower::FOLLOWER );
            }

            // loop on leader and (if necessary) follower side
            for ( auto const &tLeaderFollower : tLeaderFollowerSwitch )
            {
                // get the treated IWG phase
                std::string tPhaseName = read_phase_name( tIWGParameter, tLeaderFollower );
                tIWG->set_phase_name( tPhaseName, tLeaderFollower );

                // get dof type list from phase - ignore double-sided side sets
                if ( tLeaderFollower == mtk::Leader_Follower::LEADER )
                {
                    mPhases[ tPhaseName ]->add_dof_type_to_list( tResDofTypes );
                }

                // set properties
                for ( auto const &[ tProperty, tPropertyString ] : read_properties( tIWGParameter, tLeaderFollower ) )
                {
                    tIWG->set_property( tProperty, tPropertyString, tLeaderFollower );
                }

                // set material model
                std::string const               tMMKey       = get_leader_follower_key( "material_model", tLeaderFollower );
                Vector< Vector< std::string > > tMMNamesPair = string_to_cell_of_cell< std::string >( tIWGParameter.get< std::string >( tMMKey ) );
                MORIS_ERROR( tMMNamesPair.size() <= 1, "Model_Initializer_Phasebased::create_IWGs_without_phase() - Only one material model per CM allowed." );

                // loop over material model
                for ( uint iMM = 0; iMM < tMMNamesPair.size(); iMM++ )
                {
                    // get the MM name
                    std::string tMMName   = tMMNamesPair( iMM )( 0 );
                    std::string tMMString = tMMNamesPair( iMM )( 1 );
                    // get MM from phase
                    std::shared_ptr< Material_Model > tMM = mPhases[ tPhaseName ]->get_MM_by_name( tMMName );

                    // set CM for IWG
                    tIWG->set_material_model( tMM, tMMString, tLeaderFollower );
                }

                for ( auto const &[ tConstitutiveModel, tConstitutiveModelString ] : read_constitutive_models( tIWGParameter, tLeaderFollower ) )
                {
                    tIWG->set_constitutive_model( tConstitutiveModel, tConstitutiveModelString, tLeaderFollower );
                }

                const Vector< Vector< MSI::Dof_Type > > &tDofTypes = mPhases[ tPhaseName ]->get_dof_type_list();
                tIWG->set_dof_type_list( tDofTypes, tLeaderFollower );

                const Vector< Vector< gen::PDV_Type > > &tPdvTypes = mPhases[ tPhaseName ]->get_dv_type_list();
                tIWG->set_dv_type_list( tPdvTypes, tLeaderFollower );
            }

            // set stabilization parameters
            for ( auto const &[ tStabilizationParameter, tStabilizationParameterString ] : read_stabilization_parameters( tIWGParameter ) )
            {
                tIWG->set_stabilization_parameter( tStabilizationParameter, tStabilizationParameterString );
            }
        }
    }

    void Model_Initializer_Phasebased::create_iqis()
    {
        // create an IQI factory
        IQI_Factory             tIQIFactory;
        Vector< ParameterList > tIQIParameterList = mParameterList( FEM_PARAMETER_IQI_INDEX );
        for ( auto const &tIQIParameter : tIQIParameterList )
        {
            // get the IQI type from parameter list
            auto const             tIQIType = static_cast< IQI_Type >( tIQIParameter.get< uint >( "IQI_type" ) );
            std::shared_ptr< IQI > tIQI     = tIQIFactory.create_IQI( tIQIType );

            std::string tIQIName = tIQIParameter.get< std::string >( "IQI_name" );
            tIQI->set_name( tIQIName );
            mIQIs[ tIQIName ] = tIQI;

            // get the quantity dof type from parameter list
            Vector< moris::MSI::Dof_Type > tQuantityDofTypes = string_to_cell( tIQIParameter.get< std::string >( "dof_quantity" ), mMSIDofTypeMap );
            tIQI->set_quantity_dof_type( tQuantityDofTypes );

            // get the field index from parameter list
            sint tIQIFieldIndex = tIQIParameter.get< moris::sint >( "vectorial_field_index" );
            tIQI->set_output_type_index( tIQIFieldIndex );

            tIQI->set_normalization_type( tIQIParameter.get< std::string >( "normalization" ) );

            // set function parameters
            Vector< moris::Matrix< DDRMat > > tFuncParameters = string_to_cell_mat_2< DDRMat >( tIQIParameter.get< std::string >( "function_parameters" ) );
            tIQI->set_parameters( tFuncParameters );

            // get the treated IQI bulk type
            auto const tIQIBulkType = static_cast< Element_Type >( tIQIParameter.get< uint >( "IQI_bulk_type" ) );
            tIQI->set_bulk_type( tIQIBulkType );

            // in the case of bulk or side sets, the following loop has to be done only for the leader side.
            // if the element type is a double side set or nonconformal side set, the loop is done for both leader and follower side
            Vector< mtk::Leader_Follower > tLeaderFollowerSwitch = { mtk::Leader_Follower::LEADER };
            if ( ( tIQIBulkType == Element_Type::DOUBLE_SIDESET ) || ( tIQIBulkType == Element_Type::NONCONFORMAL_SIDESET ) )
            {
                tLeaderFollowerSwitch.push_back( mtk::Leader_Follower::FOLLOWER );
            }

            // loop on leader and (if necessary) follower side
            for ( auto const &tLeaderFollower : tLeaderFollowerSwitch )
            {
                // get the treated IWG phase
                std::string tPhaseName = read_phase_name( tIQIParameter, tLeaderFollower );
                tIQI->set_phase_name( tPhaseName, tLeaderFollower );

                // get dof type list from phase
                const Vector< Vector< MSI::Dof_Type > > &tDofTypes = mPhases[ tPhaseName ]->get_dof_type_list();
                tIQI->set_dof_type_list( tDofTypes );

                // get dof type list from phase
                const Vector< Vector< gen::PDV_Type > > &tDvTypes = mPhases[ tPhaseName ]->get_dv_type_list();
                tIQI->set_dv_type_list( tDvTypes );

                // set properties
                for ( auto const &[ tProperty, tPropertyString ] : read_properties( tIQIParameter, tLeaderFollower ) )
                {
                    tIQI->set_property( tProperty, tPropertyString, tLeaderFollower );
                }

                // set constitutive model
                for ( auto const &[ tConstitutiveModel, tConstitutiveModelString ] : read_constitutive_models( tIQIParameter, tLeaderFollower ) )
                {
                    tIQI->set_constitutive_model( tConstitutiveModel, tConstitutiveModelString, tLeaderFollower );
                }
            }    // end loop on leader and follower

            // set stabilization parameters
            for ( auto const &[ tStabilizationParameter, tStabilizationParameterString ] : read_stabilization_parameters( tIQIParameter ) )
            {
                tIQI->set_stabilization_parameter( tStabilizationParameter, tStabilizationParameterString );
            }
        }
    }

    void Model_Initializer_Phasebased::create_set_info()
    {
        // get the IWG and IQI parameter lists
        Vector< ParameterList > tIWGParameterList         = mParameterList( FEM_PARAMETER_IWG_INDEX );
        Vector< ParameterList > tIQIParameterList         = mParameterList( FEM_PARAMETER_IQI_INDEX );
        ParameterList           tComputationParameterList = mParameterList( FEM_PARAMETER_COMPUTATION_INDEX )( 0 );

        // bool true for analytical forward analysis, false for finite difference
        // decide if dRdu and dQIdu are computed by A/FD
        bool       tIsAnalyticalFA   = tComputationParameterList.get< bool >( "is_analytical_forward" );
        auto const tFDSchemeForFA    = static_cast< FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme_forward" ) );
        auto const tFDPerturbationFA = tComputationParameterList.get< real >( "finite_difference_perturbation_size_forward" );

        // get bool for analytical/finite difference for sensitivity analysis
        // decide if dRdp and dQIdp are computed by A/FD
        auto const tIsAnalyticalSA   = tComputationParameterList.get< bool >( "is_analytical_sensitivity" );
        auto const tFDSchemeForSA    = static_cast< FDScheme_Type >( tComputationParameterList.get< uint >( "finite_difference_scheme" ) );
        auto const tFDPerturbationSA = tComputationParameterList.get< real >( "finite_difference_perturbation_size" );

        // get enum for perturbation strategy for finite difference
        auto const tPerturbationStrategy = static_cast< Perturbation_Type >( tComputationParameterList.get< uint >( "finite_difference_perturbation_strategy" ) );

        // loop over the IWGs
        for ( auto const &tIWGParameter : tIWGParameterList )
        {
            auto const tIWGName = tIWGParameter.get< std::string >( "IWG_name" );
            auto const tIWG     = mIWGs[ tIWGName ];

            // get mesh set names for IWG
            Vector< std::string > tMeshSetNames;
            this->get_mesh_set_names(
                    tIWG->get_bulk_type(),
                    tIWG->get_phase_name( mtk::Leader_Follower::LEADER ),
                    tIWG->get_phase_name( mtk::Leader_Follower::FOLLOWER ),
                    tIWGParameter.get< std::string >( "neighbor_phases" ),
                    tIWGParameter.get< std::string >( "side_ordinals" ),
                    tIWG->get_ghost_flag(),
                    tMeshSetNames );

            bool const tTimeContinuity = tIWG->get_time_continuity();
            bool const tTimeBoundary   = tIWG->get_time_boundary();

            // get a representative DoF type
            MSI::Dof_Type const tFirstResidualDofType = tIWG->get_residual_dof_type()( 0 )( 0 );

            // loop over the mesh set names
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
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );
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

        // loop over the IQIs
        for ( auto const &tIQIParameter : tIQIParameterList )
        {
            auto const tIQIName = tIQIParameter.get< std::string >( "IWG_name" );
            auto const tIQI     = mIQIs[ tIQIName ];

            Vector< std::string > tMeshSetNames;
            this->get_mesh_set_names(
                    tIQI->get_bulk_type(),
                    tIQI->get_phase_name( mtk::Leader_Follower::LEADER ),
                    tIQI->get_phase_name( mtk::Leader_Follower::FOLLOWER ),
                    tIQIParameter.get< std::string >( "neighbor_phases" ),
                    tIQIParameter.get< std::string >( "side_ordinals" ),
                    false,
                    tMeshSetNames );

            // get time continuity flag
            bool const tTimeContinuity = tIQI->get_time_continuity();
            bool       tTimeBoundary   = tIQI->get_time_boundary();

            // loop over the mesh set names
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
                    aSetUserInfo.set_finite_difference_perturbation_size( tFDPerturbationSA );
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

    void
    Model_Initializer_Phasebased::get_mesh_set_names(
            Element_Type           aIWGBulkType,
            const std::string     &aLeaderPhaseName,
            const std::string     &aFollowerPhaseName,
            const std::string     &aFollowerPhaseString,
            const std::string     &aOrdinalString,
            bool                   aIsGhost,
            Vector< std::string > &aMeshSetNames )
    {
        // get the leader phase mesh index
        moris::Matrix< moris::IndexMat > tLeaderPhaseIndices = mPhases[ aLeaderPhaseName ]->get_phase_indices();

        // get the number of leader phase mesh indices
        uint tNumLeaderIndices = tLeaderPhaseIndices.numel();

        // switch on the element type
        switch ( aIWGBulkType )
        {
            case Element_Type::BULK:
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
            case Element_Type::SIDESET:
            {
                // get neighbor phase names from string
                Vector< std::string > tFollowerPhaseNames;
                string_to_cell( aFollowerPhaseString, tFollowerPhaseNames );

                // get number of neighbor phase
                uint tNumSingle = tFollowerPhaseNames.size();

                // get ordinals for boundary from string
                Matrix< DDSMat > tOrdinals;
                string_to_mat( aOrdinalString, tOrdinals );
                uint tNumBoundary = tOrdinals.numel();

                // loop over leader phase mesh indices
                for ( uint iLeaderMeshIndex = 0; iLeaderMeshIndex < tNumLeaderIndices; iLeaderMeshIndex++ )
                {
                    // get single sideset
                    for ( uint iSingle = 0; iSingle < tNumSingle; iSingle++ )
                    {
                        // get the neighbor phase name
                        std::string tNeighborPhaseName = tFollowerPhaseNames( iSingle );

                        // get the follower phase mesh index
                        moris::Matrix< moris::IndexMat > tFollowerPhaseIndices = mPhases[ tNeighborPhaseName ]->get_phase_indices();

                        // get number of neighbor phase mesh indices
                        uint tNumNeighborIndices = tFollowerPhaseIndices.numel();

                        for ( uint iNeighborMeshIndex = 0; iNeighborMeshIndex < tNumNeighborIndices; iNeighborMeshIndex++ )
                        {
                            // FIXME get this info from the mesh
                            // add mesh set name to list
                            aMeshSetNames.push_back(
                                    "iside_b0_" + std::to_string( tLeaderPhaseIndices( iLeaderMeshIndex ) ) + "_b1_" + std::to_string( tFollowerPhaseIndices( iNeighborMeshIndex ) ) );
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
            case Element_Type::DOUBLE_SIDESET:
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
                    moris::Matrix< moris::IndexMat > tFollowerPhaseIndices = mPhases[ aFollowerPhaseName ]->get_phase_indices();

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
            default:
            {
                MORIS_ERROR( false, "Model_Initializer_Phasebased::get_mesh_set_names - Unknown set type" );
            }
        }
    }

    std::string Model_Initializer_Phasebased::read_phase_name( ParameterList const &aParameterList, mtk::Leader_Follower const aLeaderFollower ) const
    {
        std::string const tKey       = get_leader_follower_key( "phase_name", aLeaderFollower );    // returns either "phase_name" or "leader_phase_name" or "follower_phase_name"
        std::string const tPhaseName = aParameterList.get< std::string >( tKey );
        ensure_existing_parameter( mPhases, tPhaseName, tKey );
        return tPhaseName;
    }

    Vector< std::pair< std::shared_ptr< Constitutive_Model >, std::string > >
    Model_Initializer_Phasebased::read_constitutive_models( ParameterList const &aParameterList, mtk::Leader_Follower const &tLeaderFollower )
    {
        // set constitutive models
        std::string const                                                         tCMKey       = get_leader_follower_key( "constitutive_models", tLeaderFollower );
        Vector< Vector< std::string > >                                           tCMNamesPair = string_to_cell_of_cell< std::string >( aParameterList.get< std::string >( tCMKey ) );
        std::string const                                                         tPhaseName   = read_phase_name( aParameterList, tLeaderFollower );
        Vector< std::pair< std::shared_ptr< Constitutive_Model >, std::string > > tCMPairs;
        // loop over CM names
        for ( uint iCM = 0; iCM < tCMNamesPair.size(); iCM++ )
        {
            std::string const tCMName   = tCMNamesPair( iCM )( 0 );
            std::string const tCMString = tCMNamesPair( iCM )( 1 );
            // get CM from phase
            std::shared_ptr< Constitutive_Model > const tCM = mPhases[ tPhaseName ]->get_CM_by_name( tCMName );
            tCMPairs.push_back( std::make_pair( tCM, tCMString ) );
        }
        return tCMPairs;
    }
    Vector< std::pair< std::shared_ptr< Stabilization_Parameter >, std::string > >
    Model_Initializer_Phasebased::read_stabilization_parameters( ParameterList const &aParameterList )
    {
        Vector< Vector< std::string > > tSPNamesPair = string_to_cell_of_cell< std::string >( aParameterList.get< std::string >( "stabilization_parameters" ) );

        Vector< std::pair< std::shared_ptr< Stabilization_Parameter >, std::string > > tSPPairs;
        for ( uint iSP = 0; iSP < tSPNamesPair.size(); iSP++ )
        {
            std::string tSPName   = tSPNamesPair( iSP )( 0 );
            std::string tSPString = tSPNamesPair( iSP )( 1 );
            ensure_existing_parameter( mStabilizationParameters, tSPName, "stabilization_parameters" );
            tSPPairs.push_back( std::make_pair( mStabilizationParameters[ tSPName ], tSPString ) );
        }
        return tSPPairs;
    }

    void Model_Initializer_Phasebased::print_physics_model()
    {
        ParameterList tComputationParameterList = this->mParameterList( 5 )( 0 );
        bool          tPrintPhysics             = tComputationParameterList.get< bool >( "print_physics_model" );
        if ( tPrintPhysics && par_rank() == 0 )
        {
            std::cout << "Phase info \n";
            for ( auto &[ _, tPhaseInfo ] : mPhases )
            {
                std::cout << "%-------------------------------------------------\n";
                tPhaseInfo->print_names();
                std::cout << "%-------------------------------------------------\n";
            }
            std::cout << " \n";
        }
        Model_Initializer::print_physics_model();
    }

}    // namespace moris::fem
