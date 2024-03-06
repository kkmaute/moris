/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_EvaluableTerm.cpp
 *
 */
#include "cl_FEM_EvaluableTerm.hpp"
#include "cl_FEM_Stabilization_Parameter.hpp"

namespace moris::fem
{
    EvaluableSideInformation const &EvaluableTerm::get_side( Leader_Follower const aIsLeader ) const
    {
        MORIS_ASSERT(
                aIsLeader == Leader_Follower::LEADER || aIsLeader == Leader_Follower::FOLLOWER,
                "EvaluableTerm: Can only request leader or follower side information." );
        return aIsLeader == Leader_Follower::LEADER ? mLeaderSideInfo : mFollowerSideInfo;
    }

    EvaluableSideInformation &EvaluableTerm::get_side( Leader_Follower const aIsLeader )
    {
        MORIS_ASSERT(
                aIsLeader == Leader_Follower::LEADER || aIsLeader == Leader_Follower::FOLLOWER,
                "EvaluableTerm: Can only request leader or follower side information." );
        return aIsLeader == Leader_Follower::LEADER ? mLeaderSideInfo : mFollowerSideInfo;
    }

    void EvaluableTerm::set_fem_set( Set *aSet )
    {
        mSet = aSet;
        get_side( Leader_Follower::LEADER ).set_fem_set( aSet );
        get_side( Leader_Follower::FOLLOWER ).set_fem_set( aSet );
        for ( auto &[ _, tStabilizationParameter ] : mStabilizationParameter )
        {
            if ( tStabilizationParameter != nullptr )
            {
                tStabilizationParameter->set_set_pointer( aSet );
            }
        }
    }

    void EvaluableTerm::set_field_interpolator_manager(
            Field_Interpolator_Manager *aFieldInterpolatorManager,
            mtk::Leader_Follower        aIsLeader )
    {
        get_side( aIsLeader ).set_fi_manager( aFieldInterpolatorManager );
        for ( auto &[ _, tStabilizationParameter ] : mStabilizationParameter )
        {
            if ( tStabilizationParameter != nullptr )
            {
                tStabilizationParameter->set_field_interpolator_manager( aFieldInterpolatorManager, aIsLeader );
            }
        }

        // TODO @ff the set update was done in the previous implementation. Is this necessary?
        get_side( aIsLeader ).set_fem_set( mSet );
        this->set_fem_set( mSet );
    }

    void EvaluableTerm::set_field_interpolator_manager_eigen_vector( Field_Interpolator_Manager *aFieldInterpolatorManager, Leader_Follower aIsLeader )
    {
        MORIS_ASSERT( aIsLeader == Leader_Follower::LEADER, "EvaluableTerm: Can only set the eigen vector field interpolator manager for the leader side." );
        get_side( aIsLeader ).set_eigen_vector_fi_manager( aFieldInterpolatorManager );
    }

    void EvaluableTerm::set_field_interpolator_manager_previous_time( Field_Interpolator_Manager *aFieldInterpolatorManager, Leader_Follower aIsLeader )
    {
        MORIS_ASSERT( aIsLeader == Leader_Follower::LEADER, "EvaluableTerm: Can only set the previous time field interpolator manager for the leader side." );
        get_side( aIsLeader ).set_previous_fi_manager( aFieldInterpolatorManager );
    }

    void EvaluableTerm::set_normal( Matrix< DDRMat > const &aNormal )
    {
        mNormal = aNormal;
        for ( auto &[ _, tStabilizationParameter ] : mStabilizationParameter )
        {
            if ( tStabilizationParameter != nullptr )
            {
                tStabilizationParameter->set_normal( mNormal );
            }
        }
    }

    void EvaluableTerm::set_property( std::shared_ptr< Property > aProperty, std::string aPropertyString, Leader_Follower aIsLeader )
    {
        get_side( aIsLeader ).set_property( aProperty, aPropertyString );
    }

    void EvaluableTerm::set_material_model( std::shared_ptr< Material_Model > aMaterialModel, std::string aMaterialModelName, Leader_Follower aIsLeader )
    {
        get_side( aIsLeader ).set_material_model( aMaterialModel, aMaterialModelName );
    }

    void EvaluableTerm::set_constitutive_model( std::shared_ptr< Constitutive_Model > aConstitutiveModel, std::string aConstitutiveModelstring, Leader_Follower aIsLeader )
    {
        get_side( aIsLeader ).set_constitutive_model( aConstitutiveModel, aConstitutiveModelstring );
    }

    void EvaluableTerm::set_stabilization_parameter( std::shared_ptr< fem::Stabilization_Parameter > const &aStabilizationParameter, std::string aStabilizationParameterType )
    {
        ensure_valid_option( mStabilizationParameter, aStabilizationParameterType, "stabilization_parameter" );
        mStabilizationParameter[ aStabilizationParameterType ] = aStabilizationParameter;

        // set active cluster measure on IWG flag on/off
        mIsActiveClusterMeasure = mIsActiveClusterMeasure || ( aStabilizationParameter->get_cluster_measure_tuple_list().size() > 0 );
    }

    Vector< Vector< MSI::Dof_Type > > const &EvaluableTerm::get_global_dof_type_list( mtk::Leader_Follower aIsLeader )
    {
        return get_side( aIsLeader ).get_global_dof_types( mStabilizationParameter );
    }

    Vector< Vector< gen::PDV_Type > > const &EvaluableTerm::get_global_dv_type_list( mtk::Leader_Follower aIsLeader )
    {
        return get_side( aIsLeader ).get_global_dv_types( mStabilizationParameter );
    }

    Vector< Vector< mtk::Field_Type > > const &EvaluableTerm::get_global_field_type_list( mtk::Leader_Follower aIsLeader )
    {
        return get_side( aIsLeader ).get_global_field_types( mStabilizationParameter );
    }

    Vector< Vector< MSI::Dof_Type > > const &EvaluableTerm::get_requested_dof_type_list( bool aIsStaggered, mtk::Leader_Follower aIsLeader )
    {
        return get_side( aIsLeader ).get_requested_global_dof_types( aIsStaggered, mStabilizationParameter );
    }

    void EvaluableTerm::get_non_unique_dof_dv_and_field_types( Vector< Vector< MSI::Dof_Type > > &aDofTypes, Vector< Vector< gen::PDV_Type > > &aDvTypes, Vector< Vector< mtk::Field_Type > > &aFieldTypes )
    {
        aDofTypes.resize( 2 );
        aDvTypes.resize( 2 );
        aFieldTypes.resize( 2 );
        get_side( Leader_Follower::LEADER ).get_non_unique_dof_dv_and_field_types( mStabilizationParameter, aDofTypes( 0 ), aDvTypes( 0 ), aFieldTypes( 0 ) );
        get_side( Leader_Follower::FOLLOWER ).get_non_unique_dof_dv_and_field_types( mStabilizationParameter, aDofTypes( 1 ), aDvTypes( 1 ), aFieldTypes( 1 ) );
    }

    void EvaluableTerm::print_names() const
    {
        auto tPrintSpec = []( std::string aTitle, auto aMap ) {
            for ( auto &[ tTypeName, _ ] : aMap ) { std::cout << " - " << aTitle << " " << tTypeName << "\n"; }
        };
        std::cout << "----------------------------------------------------------------\n"
                  << "EvaluableTerm: " << mName << "\n";
        tPrintSpec( "Leader property:            ", get_side( Leader_Follower::LEADER ).get_properties() );
        tPrintSpec( "Follower property:          ", get_side( Leader_Follower::FOLLOWER ).get_properties() );
        tPrintSpec( "Leader material model:      ", get_side( Leader_Follower::LEADER ).get_material_models() );
        tPrintSpec( "Follower material model:    ", get_side( Leader_Follower::FOLLOWER ).get_material_models() );
        tPrintSpec( "Leader constitutive model:  ", get_side( Leader_Follower::LEADER ).get_constitutive_models() );
        tPrintSpec( "Follower constitutive model:", get_side( Leader_Follower::FOLLOWER ).get_constitutive_models() );
        tPrintSpec( "Stabilization parameter:    ", mStabilizationParameter );
        std::cout << "----------------------------------------------------------------" << std::endl;
    }

    void EvaluableTerm::reset_eval_flags()
    {
        get_side( Leader_Follower::LEADER ).reset_eval_flags();
        get_side( Leader_Follower::FOLLOWER ).reset_eval_flags();
        std::for_each( mStabilizationParameter.begin(), mStabilizationParameter.end(), []( auto &tStabilizationParameter ) {
            if ( tStabilizationParameter.second != nullptr )
            {
                tStabilizationParameter.second->reset_eval_flags();
            }
        } );
    }
}    // namespace moris::fem