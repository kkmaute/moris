/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_EvaluableTerm.hpp
 *
 */
#pragma once
#include "cl_FEM_Enums.hpp"
#include "moris_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_assert.hpp"
#include <cl_FEM_EvaluableSideInformation.hpp>
#include <cl_FEM_Field_Interpolator_Manager.hpp>
#include <cl_MSI_Dof_Type_Enums.hpp>
#include <cl_Matrix.hpp>
#include <map>
#include <string>

using moris::mtk::Leader_Follower;

namespace moris::fem
{
    class Set;
    class Cluster;
    class EvaluableTerm
    {
      public:
        [[nodiscard]] std::string get_name() const { return mName; }
        void                      set_name( std::string const & aName ) { mName = aName; }

        [[nodiscard]] Set* get_fem_set() const { return mSet; };
        virtual void       set_fem_set( Set* mSet );

        [[nodiscard]] Cluster* get_fem_cluster() const { return mCluster; }
        void                   set_fem_cluster( Cluster* aCluster ) { mCluster = aCluster; }

        [[nodiscard]] bool get_time_continuity() const { return mTimeContinuity; }
        void               set_time_continuity( bool const aTimeContinuity ) { mTimeContinuity = aTimeContinuity; }

        [[nodiscard]] bool get_time_boundary() const { return mTimeBoundary; }
        void               set_time_boundary( bool const aTimeBoundary ) { mTimeBoundary = aTimeBoundary; }

        [[nodiscard]] Element_Type get_bulk_type() const { return mBulkType; }
        void                       set_bulk_type( Element_Type const aBulkType ) { mBulkType = aBulkType; }

        [[nodiscard]] Matrix< DDRMat > get_normal() const { return mNormal; }
        void                           set_normal( Matrix< DDRMat > const & aNormal );

        [[nodiscard]] Vector< Matrix< DDRMat > > get_parameters() const { return mParameters; }
        void                                     set_parameters( Vector< Matrix< DDRMat > > const & aParameters ) { mParameters = aParameters; }

        bool is_active_cluster_measure() const { return mIsActiveClusterMeasure; }

        void print_names() const;

        Field_Interpolator_Manager* get_field_interpolator_manager( Leader_Follower aIsLeader = Leader_Follower::LEADER ) const { return get_side( aIsLeader ).get_field_interpolator_manager(); }
        Field_Interpolator_Manager* get_leader_field_interpolator_manager() const { return get_field_interpolator_manager( Leader_Follower::LEADER ); }
        Field_Interpolator_Manager* get_follower_field_interpolator_manager() const { return get_field_interpolator_manager( Leader_Follower::FOLLOWER ); }
        void                        set_field_interpolator_manager( Field_Interpolator_Manager* aFieldInterpolatorManager, Leader_Follower aIsLeader );
        void                        set_field_interpolator_manager_eigen_vector( Field_Interpolator_Manager* aFieldInterpolatorManager, Leader_Follower aIsLeader = Leader_Follower::LEADER );
        void                        set_field_interpolator_manager_previous_time( Field_Interpolator_Manager* aFieldInterpolatorManager, Leader_Follower aIsLeader = Leader_Follower::LEADER );

        void                      set_phase_name( std::string aPhaseName, Leader_Follower aIsLeader = Leader_Follower::LEADER ) { get_side( aIsLeader ).set_phase_name( aPhaseName ); };
        [[nodiscard]] std::string get_phase_name( Leader_Follower aIsLeader = Leader_Follower::LEADER ) const { return get_side( aIsLeader ).get_phase_name(); };

        void                                            set_dof_type_list( Vector< Vector< MSI::Dof_Type > > const & aDofTypes, Leader_Follower aIsLeader = Leader_Follower::LEADER ) { get_side( aIsLeader ).set_dof_types( aDofTypes ); }
        [[nodiscard]] Vector< Vector< MSI::Dof_Type > > get_dof_type_list( Leader_Follower aIsLeader = Leader_Follower::LEADER ) const { return get_side( aIsLeader ).get_dof_types(); }
        Vector< Vector< MSI::Dof_Type > > const &       get_global_dof_type_list( mtk::Leader_Follower aIsLeader );
        Vector< Vector< MSI::Dof_Type > > const &       get_requested_dof_type_list( bool aIsStaggered, mtk::Leader_Follower aIsLeader );

        void                                            set_dv_type_list( Vector< Vector< gen::PDV_Type > > const & aDvTypes, Leader_Follower aIsLeader = Leader_Follower::LEADER ) { get_side( aIsLeader ).set_dv_types( aDvTypes ); }
        [[nodiscard]] Vector< Vector< gen::PDV_Type > > get_dv_type_list( Leader_Follower aIsLeader = Leader_Follower::LEADER ) const { return get_side( aIsLeader ).get_dv_types(); }
        Vector< Vector< gen::PDV_Type > > const &       get_global_dv_type_list( mtk::Leader_Follower aIsLeader );

        void                                              set_field_type_list( Vector< Vector< mtk::Field_Type > > const & aFieldTypes, Leader_Follower aIsLeader = Leader_Follower::LEADER ) { get_side( aIsLeader ).set_field_types( aFieldTypes ); }
        [[nodiscard]] Vector< Vector< mtk::Field_Type > > get_field_type_list( Leader_Follower aIsLeader = Leader_Follower::LEADER ) const { return get_side( aIsLeader ).get_field_types(); }
        Vector< Vector< mtk::Field_Type > > const &       get_global_field_type_list( mtk::Leader_Follower aIsLeader );

        void get_non_unique_dof_dv_and_field_types( Vector< Vector< MSI::Dof_Type > >& aDofTypes, Vector< Vector< gen::PDV_Type > >& aDvTypes, Vector< Vector< mtk::Field_Type > >& aFieldTypes );

        void                                                         set_property( std::shared_ptr< Property > aProperty, std::string aPropertyString, Leader_Follower aIsLeader );
        std::map< std::string, std::shared_ptr< Property > > const & get_properties( Leader_Follower aIsLeader ) { return get_side( aIsLeader ).get_properties(); };

        void                                                               set_material_model( std::shared_ptr< Material_Model > aMaterialModel, std::string aMaterialModelName, Leader_Follower aIsLeader );
        std::map< std::string, std::shared_ptr< Material_Model > > const & get_material_models( Leader_Follower aIsLeader ) { return get_side( aIsLeader ).get_material_models(); };

        void                                                                   set_constitutive_models( std::shared_ptr< Constitutive_Model > aConstitutiveModel, std::string aConstitutiveModelstring, Leader_Follower aIsLeader );
        std::map< std::string, std::shared_ptr< Constitutive_Model > > const & get_constitutive_models( Leader_Follower aIsLeader ) { return get_side( aIsLeader ).get_constitutive_models(); };

        void                                                                             set_stabilization_parameter( std::shared_ptr< fem::Stabilization_Parameter > const & aStabilizationParameter, std::string aStabilizationParameterString );
        std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > const & get_stabilization_parameters() const { return mStabilizationParameter; }

        virtual void reset_eval_flags();

      private:
        [[nodiscard]] EvaluableSideInformation get_side( Leader_Follower const aIsLeader ) const;

      protected:
        Set*     mSet     = nullptr;
        Cluster* mCluster = nullptr;

      private:
        std::string      mName;
        Matrix< DDRMat > mNormal;

        EvaluableSideInformation mLeaderSideInfo{ mtk::Leader_Follower::LEADER };
        EvaluableSideInformation mFollowerSideInfo{ mtk::Leader_Follower::FOLLOWER };

        Vector< Matrix< DDRMat > >                                               mParameters;
        std::map< std::string, std::shared_ptr< fem::Stabilization_Parameter > > mStabilizationParameter;

        bool mGlobalFieldBuild = true;

        // bool for time continuity and boundary
        bool mTimeContinuity = false;
        bool mTimeBoundary   = false;

        // active cluster measure on IWG flag
        bool mIsActiveClusterMeasure = false;

        // element type
        Element_Type mBulkType = Element_Type::BULK;
    };
}    // namespace moris::fem
