/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI.cpp
 *
 */

#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_norm.hpp"
#include "fn_min.hpp"
#include "fn_max.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void
        IQI::print_names()
        {
            std::cout << "----------" << std::endl;
            std::cout << "IQI: " << mName << std::endl;

            // properties
            for ( uint iProp = 0; iProp < mLeaderProp.size(); iProp++ )
            {
                if ( mLeaderProp( iProp ) != nullptr )
                {
                    std::cout << "Leader property: " << mLeaderProp( iProp )->get_name() << std::endl;
                }
            }
            for ( uint iProp = 0; iProp < mFollowerProp.size(); iProp++ )
            {
                if ( mFollowerProp( iProp ) != nullptr )
                {
                    std::cout << "Follower property:  " << mFollowerProp( iProp )->get_name() << std::endl;
                }
            }

            // CM
            for ( uint iCM = 0; iCM < mLeaderCM.size(); iCM++ )
            {
                if ( mLeaderCM( iCM ) != nullptr )
                {
                    std::cout << "Leader CM: " << mLeaderCM( iCM )->get_name() << std::endl;
                }
            }
            for ( uint iCM = 0; iCM < mFollowerCM.size(); iCM++ )
            {
                if ( mFollowerCM( iCM ) != nullptr )
                {
                    std::cout << "Follower CM:  " << mFollowerCM( iCM )->get_name() << std::endl;
                }
            }

            // SP
            for ( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
            {
                if ( mStabilizationParam( iSP ) != nullptr )
                {
                    std::cout << "SP: " << mStabilizationParam( iSP )->get_name() << std::endl;
                }
            }
            std::cout << "----------" << std::endl;
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_function_pointers()
        {
            // switch on element type
            switch ( mSet->get_element_type() )
            {
                case fem::Element_Type::BULK:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_sideset;
                    break;
                }
                case fem::Element_Type::TIME_SIDESET:
                case fem::Element_Type::TIME_BOUNDARY:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_bulk;
                    break;
                }
                case fem::Element_Type::DOUBLE_SIDESET:
                case fem::Element_Type::NONCONFORMAL_SIDESET:
                {
                    m_compute_dQIdu_FD          = &IQI::select_dQIdu_FD;
                    m_compute_dQIdp_FD_material = &IQI::select_dQIdp_FD_material_double;
                    m_compute_dQIdp_FD_geometry = &IQI::select_dQIdp_FD_geometry_double;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IQI::set_function_pointers - unknown element type." );
            }

            // switch on perturbation strategy
            switch ( mSet->get_perturbation_strategy() )
            {
                case fem::Perturbation_Type::RELATIVE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_relative;
                    break;
                }
                case fem::Perturbation_Type::ABSOLUTE:
                {
                    m_build_perturbation_size = &IQI::build_perturbation_size_absolute;
                    break;
                }
                default:
                    MORIS_ERROR( false, "IQI::set_function_pointers - unknown perturbation type." );
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::reset_eval_flags()
        {
            // reset properties
            for ( const std::shared_ptr< Property >& tProp : mLeaderProp )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            for ( const std::shared_ptr< Property >& tProp : mFollowerProp )
            {
                if ( tProp != nullptr )
                {
                    tProp->reset_eval_flags();
                }
            }

            // reset constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
            {
                if ( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
            {
                if ( tCM != nullptr )
                {
                    tCM->reset_eval_flags();
                }
            }

            // reset stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    tSP->reset_eval_flags();
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_phase_name(
                std::string          aPhaseName,
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderPhaseName = aPhaseName;
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    mFollowerPhaseName = aPhaseName;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_phase_name - aIsLeader can only be leader or follower." );
                }
            }
        }

        //------------------------------------------------------------------------------

        std::string
        IQI::get_phase_name( mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    return mLeaderPhaseName;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    return mFollowerPhaseName;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::get_phase_name - aIsLeader can only be leader or follower." );
                    return mLeaderPhaseName;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_reference_value( real aReferenceValue )
        {
            if ( not mNormalized )
            {
                mReferenceValue = aReferenceValue;
                mNormalized     = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_field_interpolator_manager(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderFIManager = aFieldInterpolatorManager;
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    mFollowerFIManager = aFieldInterpolatorManager;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_field_interpolator_manager - can only be leader or follower" );
                }
            }

            // loop over the the SP
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : this->get_stabilization_parameters() )
            {
                if ( tSP != nullptr )
                {
                    // set the field interpolator manager for the SP
                    tSP->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ), aIsLeader );

                    // set the fem set pointer for the SP
                    tSP->set_set_pointer( mSet );
                }
            }

            // loop over the constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : this->get_constitutive_models( aIsLeader ) )
            {
                if ( tCM != nullptr )
                {
                    // set the field interpolator manager for the CM
                    tCM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ) );

                    // set the fem set pointe for the CM
                    tCM->set_set_pointer( mSet );
                }
            }

            // loop over the properties
            for ( const std::shared_ptr< Property >& tProp : this->get_properties( aIsLeader ) )
            {
                if ( tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ) );

                    // set the fem set pointer for the property
                    tProp->set_set_pointer( mSet );
                }
            }
        }

        //------------------------------------------------------------------------------

        Field_Interpolator_Manager*
        IQI::get_field_interpolator_manager(
                mtk::Leader_Follower aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                    return mLeaderFIManager;

                case mtk::Leader_Follower::FOLLOWER:
                    return mFollowerFIManager;

                default:
                    MORIS_ERROR( false, "IQI::get_field_inetrpolator_manager - can only be leader or follower." );
                    return mLeaderFIManager;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_dof_type_list(
                const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
                mtk::Leader_Follower                               aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderDofTypes = aDofTypes;
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    mFollowerDofTypes = aDofTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dof_type_list - can only be LEADER or FOLLOWER." );
                }
            }
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< MSI::Dof_Type > >&
        IQI::get_dof_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dof type list
                    return mLeaderDofTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dof type list
                    return mFollowerDofTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_dof_type_list - can only be leader or follower." );
                    return mLeaderDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< MSI::Dof_Type > >&
        IQI::get_global_dof_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            // if the global list was not yet built
            if ( mGlobalDofBuild )
            {
                // build global dof type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDofBuild   = false;
                mGlobalDvBuild    = false;
                mGlobalFieldBuild = false;
            }

            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dof type list
                    return mLeaderGlobalDofTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dof type list
                    return mFollowerGlobalDofTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_global_dof_type_list - can only be leader or follower." );
                    return mLeaderGlobalDofTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_dv_type_list(
                const Vector< Vector< PDV_Type > >& aDvTypes,
                mtk::Leader_Follower                          aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderDvTypes = aDvTypes;
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    mFollowerDvTypes = aDvTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dv_type_list - can only be LEADER or FOLLOWER." );
                }
            }
        }

        //------------------------------------------------------------------------------

        Vector< Vector< PDV_Type > >&
        IQI::get_dv_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dv type list
                    return mLeaderDvTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dv type list
                    return mFollowerDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_dv_type_list - can only be leader or follower." );
                    return mLeaderDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< PDV_Type > >&
        IQI::get_global_dv_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            // if the global list was not yet built
            if ( mGlobalDvBuild )
            {
                // build global dv type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDvBuild    = false;
                mGlobalDofBuild   = false;
                mGlobalFieldBuild = false;
            }

            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dv type list
                    return mLeaderGlobalDvTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dv type list
                    return mFollowerGlobalDvTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_global_dv_type_list - can only be leader or follower." );
                    return mLeaderGlobalDvTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< mtk::Field_Type > >&
        IQI::get_global_field_type_list(
                mtk::Leader_Follower aIsLeader )
        {
            // if the global list was not yet built
            if ( mGlobalFieldBuild )
            {
                // build global dv type list
                this->build_global_dof_dv_and_field_type_list();

                // update build flag
                mGlobalDvBuild    = false;
                mGlobalDofBuild   = false;
                mGlobalFieldBuild = false;
            }

            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dof type list
                    return mLeaderGlobalFieldTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dof type list
                    return mFollowerGlobalFieldTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_global_field_type_list - can only be leader or follower." );
                    return mLeaderGlobalFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_field_type_list(
                const Vector< Vector< mtk::Field_Type > >& aDofTypes,
                mtk::Leader_Follower                                 aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderFieldTypes = aDofTypes;
                    break;
                }
                case mtk::Leader_Follower::FOLLOWER:
                {
                    mFollowerFieldTypes = aDofTypes;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IQI::set_dof_type_list - can only be LEADER or FOLLOWER." );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_field_interpolator_manager_eigen_vector(
                Field_Interpolator_Manager* aFieldInterpolatorManager,
                mtk::Leader_Follower        aIsLeader )
        {
            switch ( aIsLeader )
            {
                case mtk::Leader_Follower::LEADER:
                {
                    mLeaderEigenFIManager = aFieldInterpolatorManager;
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IQI::set_field_interpolator_manager_eigen_vector - can only be leader" );
                }
            }
        }

        //------------------------------------------------------------------------------

        const Vector< Vector< mtk::Field_Type > >&
        IQI::get_field_type_list(
                mtk::Leader_Follower aIsLeader ) const
        {
            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader global dof type list
                    return mLeaderFieldTypes;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower global dof type list
                    return mFollowerFieldTypes;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_dof_type_list - can only be leader or follower." );
                    return mLeaderFieldTypes;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Leader_Follower        aIsLeader )
        {
            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "IQI::set_property - IQI %s - Unknown aPropertyString: %s ",
                    mName.c_str(),
                    aPropertyString.c_str() );

            // set the property in the property pointer cell
            this->get_properties( aIsLeader )( mPropertyMap[ aPropertyString ] ) = aProperty;
        }

        //------------------------------------------------------------------------------

        Vector< std::shared_ptr< fem::Property > >&
        IQI::get_properties(
                mtk::Leader_Follower aIsLeader )
        {
            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader property pointers
                    return mLeaderProp;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower property pointers
                    return mFollowerProp;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_properties - can only be leader or follower." );
                    return mLeaderProp;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Leader_Follower                  aIsLeader )
        {
            // check that aConstitutiveString makes sense
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                    "IQI::set_constitutive_model - IQI %s - Unknown aConstitutiveString: %s ",
                    mName.c_str(),
                    aConstitutiveString.c_str() );

            // set the CM in the CM pointer cell
            this->get_constitutive_models( aIsLeader )( mConstitutiveMap[ aConstitutiveString ] ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        Vector< std::shared_ptr< fem::Constitutive_Model > >&
        IQI::get_constitutive_models(
                mtk::Leader_Follower aIsLeader )
        {
            // switch on leader/follower
            switch ( aIsLeader )
            {
                // if leader
                case mtk::Leader_Follower::LEADER:
                {
                    // return leader property pointers
                    return mLeaderCM;
                }
                // if follower
                case mtk::Leader_Follower::FOLLOWER:
                {
                    // return follower property pointers
                    return mFollowerCM;
                }
                // if none
                default:
                {
                    MORIS_ASSERT( false, "IQI::get_constitutive_models - can only be leader or follower." );
                    return mLeaderCM;
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                    "IQI::set_stabilization_parameter - IQI %s - Unknown aStabilizationString: %s ",
                    mName.c_str(),
                    aStabilizationString.c_str() );

            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( mStabilizationMap[ aStabilizationString ] ) = aStabilizationParameter;

            // set active cluster measure on IQI flag on/off
            mActiveCMEAFlag = mActiveCMEAFlag || ( aStabilizationParameter->get_cluster_measure_tuple_list().size() > 0 );
        }

        //------------------------------------------------------------------------------

        void
        IQI::get_non_unique_dof_dv_and_field_types(
                Vector< Vector< MSI::Dof_Type > >&   aDofTypes,
                Vector< Vector< PDV_Type > >&        aDvTypes,
                Vector< Vector< mtk::Field_Type > >& aFieldTypes )

        {
            // init counters for dof and dv types
            uint tLeaderDofCounter     = 0;
            uint tFollowerDofCounter   = 0;
            uint tLeaderDvCounter      = 0;
            uint tFollowerDvCounter    = 0;
            uint tLeaderFieldCounter   = 0;
            uint tFollowerFieldCounter = 0;

            // get number of direct leader dof dependencies
            for ( uint iDof = 0; iDof < mLeaderDofTypes.size(); iDof++ )
            {
                tLeaderDofCounter += mLeaderDofTypes( iDof ).size();
            }

            // get number of direct leader dv dependencies
            for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
            {
                tLeaderDvCounter += mLeaderDvTypes( iDv ).size();
            }

            // get number of direct leader field dependencies
            for ( uint iFi = 0; iFi < mLeaderFieldTypes.size(); iFi++ )
            {
                tLeaderFieldCounter += mLeaderFieldTypes( iFi ).size();
            }

            // get number of direct follower dof dependencies
            for ( uint iDof = 0; iDof < mFollowerDofTypes.size(); iDof++ )
            {
                tFollowerDofCounter += mFollowerDofTypes( iDof ).size();
            }

            // get number of direct follower dv dependencies
            for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
            {
                tFollowerDvCounter += mFollowerDvTypes( iDv ).size();
            }

            // get number of direct follower field dependencies
            for ( uint iFi = 0; iFi < mFollowerFieldTypes.size(); iFi++ )
            {
                tFollowerFieldCounter += mFollowerFieldTypes( iFi ).size();
            }

            // loop over the leader properties
            for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tLeaderDofCounter += tActiveDofTypes.size();
                    tLeaderDvCounter += tActiveDvTypes.size();
                    tLeaderFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over follower properties
            for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type lists
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counter
                    tFollowerDofCounter += tActiveDofTypes.size();
                    tFollowerDvCounter += tActiveDvTypes.size();
                    tFollowerFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over leader constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;
                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tLeaderDofCounter += tActiveDofTypes.size();
                    tLeaderDvCounter += tActiveDvTypes.size();
                    tLeaderFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over follower constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;
                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // update dof and dv counters
                    tFollowerDofCounter += tActiveDofTypes.size();
                    tFollowerDvCounter += tActiveDvTypes.size();
                    tFollowerFieldCounter += tActiveFieldTypes.size();
                }
            }

            // loop over leader stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique dof type list
                    Vector< MSI::Dof_Type > tActiveDofTypes;
                    Vector< PDV_Type >      tActiveDvTypes;
                    tSP->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // update dof and dv counters
                    tLeaderDofCounter += tActiveDofTypes.size();
                    tLeaderDvCounter += tActiveDvTypes.size();
                    tFollowerDofCounter += tActiveDofTypes.size();
                    tFollowerDvCounter += tActiveDvTypes.size();
                }
            }

            // reserve memory for dof and dv type lists
            aDofTypes.resize( 2 );
            aDvTypes.resize( 2 );
            aFieldTypes.resize( 2 );

            aDofTypes( 0 ).reserve( tLeaderDofCounter );
            aDvTypes( 0 ).reserve( tLeaderDvCounter );
            aFieldTypes( 0 ).reserve( tLeaderFieldCounter );
            aDofTypes( 1 ).reserve( tFollowerDofCounter );
            aDvTypes( 1 ).reserve( tFollowerDvCounter );
            aFieldTypes( 1 ).reserve( tFollowerFieldCounter );

            // loop over leader dof direct dependencies
            for ( uint iDof = 0; iDof < mLeaderDofTypes.size(); iDof++ )
            {
                // populate the dof list
                aDofTypes( 0 ).append( mLeaderDofTypes( iDof ) );
            }

            // loop over leader dv direct dependencies
            for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
            {
                // populate the dv list
                aDvTypes( 0 ).append( mLeaderDvTypes( iDv ) );
            }

            // loop over leader field direct dependencies
            for ( uint iFi = 0; iFi < mLeaderFieldTypes.size(); iFi++ )
            {
                // populate the dv list
                aFieldTypes( 0 ).append( mLeaderFieldTypes( iFi ) );
            }

            // loop over follower dof direct dependencies
            for ( uint iDof = 0; iDof < mFollowerDofTypes.size(); iDof++ )
            {
                // populate the dof list
                aDofTypes( 1 ).append( mFollowerDofTypes( iDof ) );
            }

            // loop over follower dv direct dependencies
            for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
            {
                // populate the dv list
                aDvTypes( 1 ).append( mFollowerDvTypes( iDv ) );
            }

            // loop over follower dv direct dependencies
            for ( uint iFi = 0; iFi < mFollowerFieldTypes.size(); iFi++ )
            {
                aFieldTypes( 1 ).append( mFollowerFieldTypes( iFi ) );
            }

            // loop over leader properties
            for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aFieldTypes( 0 ).append( tActiveFieldTypes );
                }
            }

            // loop over follower properties
            for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
            {
                if ( tProperty != nullptr )
                {
                    // get property non unique dof and dv type list
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tProperty->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                    aFieldTypes( 1 ).append( tActiveFieldTypes );
                }
            }

            // loop over the leader constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aFieldTypes( 0 ).append( tActiveFieldTypes );
                }
            }

            // loop over the follower constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
            {
                if ( tCM != nullptr )
                {
                    // get CM non unique dof and dv type lists
                    Vector< MSI::Dof_Type >   tActiveDofTypes;
                    Vector< PDV_Type >        tActiveDvTypes;
                    Vector< mtk::Field_Type > tActiveFieldTypes;

                    tCM->get_non_unique_dof_dv_and_field_types(
                            tActiveDofTypes,
                            tActiveDvTypes,
                            tActiveFieldTypes );

                    // populate the dof and dv lists
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                    aFieldTypes( 1 ).append( tActiveFieldTypes );
                }
            }

            // FIXME this is potentially problematic since it will add follower dependencies even for bulk elements
            // FIXME Ask lise about it. We could ask the set for the element type. should work for DOUBLE_SIDED.
            // FIXME Whats with time boundary
            // loop over the stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get SP non unique leader dof type list
                    Vector< MSI::Dof_Type > tActiveDofTypes;
                    Vector< PDV_Type >      tActiveDvTypes;

                    tSP->get_non_unique_dof_and_dv_types(
                            tActiveDofTypes,
                            tActiveDvTypes );

                    // populate the dof and dv lists
                    aDofTypes( 0 ).append( tActiveDofTypes );
                    aDvTypes( 0 ).append( tActiveDvTypes );
                    aDofTypes( 1 ).append( tActiveDofTypes );
                    aDvTypes( 1 ).append( tActiveDvTypes );
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI::build_global_dof_dv_and_field_type_list()
        {
            // LEADER-------------------------------------------------------
            // get number of dof and dv types on set
            uint tNumDofTypes   = mSet->get_num_unique_dof_types();
            uint tNumDvTypes    = mSet->get_num_unique_dv_types();
            uint tNumFieldTypes = mSet->get_num_unique_field_types();

            // set size for the global dof and dv type lists
            mLeaderGlobalDofTypes.reserve( tNumDofTypes );
            mLeaderGlobalDvTypes.reserve( tNumDvTypes );
            mLeaderGlobalFieldTypes.reserve( tNumFieldTypes );

            // set a size for the dof and dv checkLists
            //( used to avoid repeating a dof or a dv type)
            Matrix< DDSMat > tDofCheckList( tNumDofTypes, 1, -1 );
            Matrix< DDSMat > tDvCheckList( tNumDvTypes, 1, -1 );
            Matrix< DDSMat > tFieldCheckList( tNumFieldTypes, 1, -1 );

            // get dof type from direct dependencies
            for ( uint iDof = 0; iDof < mLeaderDofTypes.size(); iDof++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( mLeaderDofTypes( iDof )( 0 ) );    // FIXME'

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mLeaderGlobalDofTypes.push_back( mLeaderDofTypes( iDof ) );
            }

            // get dv type from direct dependencies
            for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mLeaderDvTypes( iDv )( 0 ) );    // FIXME'

                // put the dv type in the checklist
                tDvCheckList( tDvTypeIndex ) = 1;

                // put the dv type in the global type list
                mLeaderGlobalDvTypes.push_back( mLeaderDvTypes( iDv ) );
            }

            // get field type from direct dependencies
            for ( uint iFi = 0; iFi < mLeaderFieldTypes.size(); iFi++ )
            {
                // get set index for field type
                sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( mLeaderFieldTypes( iFi )( 0 ) );    // FIXME'

                // put the field type in the checklist
                tFieldCheckList( tFieldTypeIndex ) = 1;

                // put the field type in the global type list
                mLeaderGlobalFieldTypes.push_back( mLeaderFieldTypes( iFi ) );
            }

            // get dof type from leader properties
            for ( const std::shared_ptr< Property >& tProperty : mLeaderProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for property
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dof enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mLeaderGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for property
                    const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                            tProperty->get_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mLeaderGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from leader constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mLeaderGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for property
                    const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                            tCM->get_global_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mLeaderGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from leader stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tSP->get_global_dof_type_list( mtk::Leader_Follower::LEADER );

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tSP->get_global_dv_type_list( mtk::Leader_Follower::LEADER );

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mLeaderGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }
                }
            }

            // reduce size of dof and dv lists to fit unique list
            mLeaderGlobalDofTypes.shrink_to_fit();
            mLeaderGlobalDvTypes.shrink_to_fit();
            mLeaderGlobalFieldTypes.shrink_to_fit();

            // FOLLOWER--------------------------------------------------------

            // set size for the global dof type list
            mFollowerGlobalDofTypes.reserve( tNumDofTypes );
            mFollowerGlobalDvTypes.reserve( tNumDvTypes );
            mFollowerGlobalFieldTypes.reserve( tNumFieldTypes );

            // set a size for the checkList ( used to avoid repeating a dof type)
            tDofCheckList.fill( -1 );
            tDvCheckList.fill( -1 );
            tFieldCheckList.fill( -1 );

            // get dof type from follower direct dependencies
            for ( uint iDof = 0; iDof < mFollowerDofTypes.size(); iDof++ )
            {
                // get set index for dof type
                sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( mFollowerDofTypes( iDof )( 0 ) );

                // put the dof type in the checklist
                tDofCheckList( tDofTypeIndex ) = 1;

                // put the dof type in the global type list
                mFollowerGlobalDofTypes.push_back( mFollowerDofTypes( iDof ) );
            }

            // get dv type from follower direct dependencies
            for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
            {
                // get set index for dv type
                sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( mFollowerDvTypes( iDv )( 0 ) );

                // put the dv type in the checklist
                tDvCheckList( tDvTypeIndex ) = 1;

                // put the dv type in the global type list
                mFollowerGlobalDvTypes.push_back( mFollowerDvTypes( iDv ) );
            }

            // get field type from follower direct dependencies
            for ( uint iFi = 0; iFi < mFollowerFieldTypes.size(); iFi++ )
            {
                // get set index for field type
                sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( mFollowerFieldTypes( iFi )( 0 ) );

                // put the field type in the check list
                tFieldCheckList( tFieldTypeIndex ) = 1;

                // put the field type in the global type list
                mFollowerGlobalFieldTypes.push_back( mFollowerFieldTypes( iFi ) );
            }

            // get dof type from leader properties
            for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
            {
                if ( tProperty != nullptr )
                {
                    // get dof types for property
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tProperty->get_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for property
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tProperty->get_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mFollowerGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for property
                    const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                            tProperty->get_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mFollowerGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from follower constitutive models
            for ( const std::shared_ptr< Constitutive_Model >& tCM : mFollowerCM )
            {
                if ( tCM != nullptr )
                {
                    // get dof types for constitutive model
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tCM->get_global_dof_type_list();

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for constitutive model
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tCM->get_global_dv_type_list();

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mFollowerGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }

                    // get field types for property
                    const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                            tCM->get_global_field_type_list();

                    // loop on property field type
                    for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                    {
                        // get set index for field type
                        sint tFieldTypeIndex = mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

                        // if field enum not in the list
                        if ( tFieldCheckList( tFieldTypeIndex ) != 1 )
                        {
                            // put the field type in the check list
                            tFieldCheckList( tFieldTypeIndex ) = 1;

                            // put the field type in the global type list
                            mFollowerGlobalFieldTypes.push_back( tActiveFieldTypes( iFi ) );
                        }
                    }
                }
            }

            // get dof type from stabilization parameters
            for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
            {
                if ( tSP != nullptr )
                {
                    // get dof types for constitutive model
                    const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                            tSP->get_global_dof_type_list( mtk::Leader_Follower::FOLLOWER );

                    // loop on property dof type
                    for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                    {
                        // get set index for dof type
                        sint tDofTypeIndex = mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                        // if dof enum not in the list
                        if ( tDofCheckList( tDofTypeIndex ) != 1 )
                        {
                            // put the dof type in the checklist
                            tDofCheckList( tDofTypeIndex ) = 1;

                            // put the dof type in the global type list
                            mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                        }
                    }

                    // get dv types for stabilization parameter
                    const Vector< Vector< PDV_Type > >& tActiveDvTypes =
                            tSP->get_global_dv_type_list( mtk::Leader_Follower::FOLLOWER );

                    // loop on property dv type
                    for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                    {
                        // get set index for dv type
                        sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                        // if dv enum not in the list
                        if ( tDvCheckList( tDvTypeIndex ) != 1 )
                        {
                            // put the dv type in the checklist
                            tDvCheckList( tDvTypeIndex ) = 1;

                            // put the dv type in the global type list
                            mFollowerGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                        }
                    }
                }
            }

            // reduce size of dof list to fit unique list
            mFollowerGlobalDofTypes.shrink_to_fit();
            mFollowerGlobalDvTypes.shrink_to_fit();
            mFollowerGlobalFieldTypes.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        IQI::build_requested_dof_type_list()
        {
            // clear the dof lists
            mRequestedLeaderGlobalDofTypes.clear();
            mRequestedFollowerGlobalDofTypes.clear();

            // get the requested dof types
            Vector< enum MSI::Dof_Type > tRequestedDofTypes = mSet->get_requested_dof_types();

            // reserve possible max size for requested dof lists
            mRequestedLeaderGlobalDofTypes.reserve( tRequestedDofTypes.size() );
            mRequestedFollowerGlobalDofTypes.reserve( tRequestedDofTypes.size() );

            // loop over the requested dof types
            for ( auto tDofTypes : tRequestedDofTypes )
            {
                // loop over the IWG leader dof types groups
                for ( uint Ik = 0; Ik < mLeaderGlobalDofTypes.size(); Ik++ )
                {
                    // if requested dof type matches IWG leader dof type
                    if ( mLeaderGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                    {
                        // add the IWG leader dof type to the requested dof list
                        mRequestedLeaderGlobalDofTypes.push_back( mLeaderGlobalDofTypes( Ik ) );
                        break;
                    }
                }

                // loop over the IWG follower dof types groups
                for ( uint Ik = 0; Ik < mFollowerGlobalDofTypes.size(); Ik++ )
                {
                    // if requested dof type matches IWG follower dof type
                    if ( mFollowerGlobalDofTypes( Ik )( 0 ) == tDofTypes )
                    {
                        // add the IWG follower dof type to the requested dof list
                        mRequestedFollowerGlobalDofTypes.push_back( mFollowerGlobalDofTypes( Ik ) );
                        break;
                    }
                }
            }

            // reduce size for requested dof lists
            mRequestedLeaderGlobalDofTypes.shrink_to_fit();
            mRequestedFollowerGlobalDofTypes.shrink_to_fit();
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tQIIndex );

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get leader number of dof types
            uint tNumLeaderDofType = mRequestedLeaderGlobalDofTypes.size();

            // reset the QI
            mSet->get_QI()( tQIIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tQIIndex );

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tNumLeaderDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iFI );

                // get leader index for residual dof type, indices for assembly
                uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );

                // get field interpolator for dependency dof type
                Field_Interpolator* tFI =
                        mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );
                            tFI->reset_eval_flags();    // not useful

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // get follower number of dof types
            uint tNumFollowerDofType = mRequestedFollowerGlobalDofTypes.size();

            // loop over the follower dof types
            for ( uint iFI = 0; iFI < tNumFollowerDofType; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tFollowerDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
                uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDepDofIndex )( tFollowerDepDofIndex, 0 );

                // get follower dependency field interpolator
                Field_Interpolator* tFI =
                        mFollowerFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // loop over the coefficients columns
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the coefficients rows
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( 0 ) * tQI /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();
                            tFI->reset_eval_flags();    // not useful

                            // reset the QI
                            mSet->get_QI()( tQIIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the dQIdu
                            mSet->get_residual()( tQIIndex )(
                                    { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter },
                                    { 0, 0 } ) +=
                                    tFDScheme( 1 )( iPoint ) *      //
                                    mSet->get_QI()( tQIIndex ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tQIIndex ) = tQIStore;
        }

        //------------------------------------------------------------------------------

        bool
        IQI::check_dQIdu_FD(
                real               aWStar,
                real               aPerturbation,
                real               aEpsilon,
                Matrix< DDRMat >&  adQIdu,
                Matrix< DDRMat >&  adQIduFD,
                bool               aErrorPrint,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // compute dQIdu with IQI
            this->compute_dQIdu( aWStar );
            adQIdu = mSet->get_residual()( tQIIndex );

            // reset dQIdu
            mSet->get_residual()( tQIIndex ).fill( 0.0 );

            // compute dQIdu by FD
            this->compute_dQIdu_FD( aWStar, aPerturbation, aFDSchemeType );
            adQIduFD = mSet->get_residual()( tQIIndex );

            // define a boolean for check
            bool tCheckdQIdu = true;

            // check if adQIdu and adQIduFD have the same size
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_rows() == adQIduFD.n_rows() );
            tCheckdQIdu = tCheckdQIdu && ( adQIdu.n_cols() == adQIduFD.n_cols() );

            // check that matrices to compare have same size
            MORIS_ERROR(
                    ( adQIdu.n_rows() == adQIduFD.n_rows() ) &&    //
                            ( adQIdu.n_cols() == adQIduFD.n_cols() ),
                    "IQI::check_dQIdu - matrices to check do not share same dimensions." );

            // define a real for absolute difference
            real tAbsolute = 0.0;

            // define a real for relative difference
            real tRelative = 0.0;

            // loop over the rows
            for ( uint iRow = 0; iRow < adQIdu.n_rows(); iRow++ )
            {
                // loop over the columns
                for ( uint jCol = 0; jCol < adQIdu.n_cols(); jCol++ )
                {
                    // get absolute difference
                    tAbsolute = std::abs( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) );

                    // get relative difference
                    tRelative = std::abs( ( adQIduFD( iRow, jCol ) - adQIdu( iRow, jCol ) ) / adQIduFD( iRow, jCol ) );

                    // update check value
                    tCheckdQIdu = tCheckdQIdu && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                    // debug print
                    if ( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                    {
                        if ( aErrorPrint )
                        {
                            std::cout << "iRow " << iRow << " - jCol " << jCol << "\n"
                                      << std::flush;
                            std::cout << "adQIdu( iRow, jCol )    " << std::setprecision( 12 ) << adQIdu( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "adQIduFD( iRow, jCol )  " << std::setprecision( 12 ) << adQIduFD( iRow, jCol ) << "\n"
                                      << std::flush;
                            std::cout << "Absolute difference " << tAbsolute << "\n"
                                      << std::flush;
                            std::cout << "Relative difference " << tRelative << "\n"
                                      << std::flush;
                        }
                    }
                }
            }
            // return bool
            return tCheckdQIdu;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_material(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the column index to assemble in residual
            sint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( aFDSchemeType, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get the requested ip pdv types
            Vector< Vector< PDV_Type > > tRequestedPdvTypes;
            mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

            // get number of requested dv types
            uint tNumPdvType = tRequestedPdvTypes.size();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // loop over the pdv types
            for ( uint iFI = 0; iFI < tNumPdvType; iFI++ )
            {
                // get the FI for the dv type
                Field_Interpolator* tFI =
                        mLeaderFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

                // get number of leader FI bases and fields
                uint tDerNumBases  = tFI->get_number_of_space_time_bases();
                uint tDerNumFields = tFI->get_number_of_fields();

                // coefficients for dof type wrt which derivative is computed
                Matrix< DDRMat > tCoeff = tFI->get_coeff();

                // init pdv coeff counter
                uint tPdvCoeffCounter = 0;

                // loop over the pdv coefficient column
                for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
                {
                    // loop over the pdv coefficient row
                    for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                    {
                        // compute the perturbation absolute value
                        real tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // get the pdv index for assembly
                        uint tPdvAssemblyIndex = mSet->get_mat_pdv_assembly_map()( iFI )( 0, 0 ) + tPdvCoeffCounter;

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over the points for FD
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // set the perturbed coefficients to FI
                            tFI->set_coeff( tCoeffPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            this->compute_QI( aWStar );

                            // assemble the jacobian
                            mSet->get_dqidpmat()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update the pdv coefficient counter
                        tPdvCoeffCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpmat()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_material - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            return ( this->*m_build_perturbation_size )(
                    aPerturbation,
                    aCoefficientToPerturb,
                    aMaxPerturbation,
                    aTolerance );
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_relative(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_relative - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // compute the perturbation value using fraction of maximum allowable perturbation
            real tDeltaH = std::abs( aPerturbation * aMaxPerturbation );

            // compute perturbation such that it is not smaller than tolerance
            // and not larger than maximum value
            tDeltaH = std::max( std::min( tDeltaH, aMaxPerturbation ), tActualTol );

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        real
        IQI::build_perturbation_size_absolute(
                const real& aPerturbation,
                const real& aCoefficientToPerturb,
                const real& aMaxPerturbation,
                const real& aTolerance )
        {
            // check that maximum perturbation size is larger than tolerance
            MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                    "IQI::build_perturbation_size_absolute - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                    aMaxPerturbation,
                    aTolerance );

            // determine actual tolerance (only useful when above assert inactive)
            const real tActualTol = std::min( aMaxPerturbation, aTolerance );

            // check that absolute value of perturbation is not smaller than tolerance
            // and not larger than maximum value
            return std::max( std::min( std::abs( aPerturbation ), aMaxPerturbation ), tActualTol );
        }

        //------------------------------------------------------------------------------

        real
        IQI::check_ig_coordinates_inside_ip_element(
                const real&         aPerturbation,
                const real&         aCoefficientToPerturb,
                const uint&         aSpatialDirection,
                fem::FDScheme_Type& aUsedFDScheme )
        {
            // FIXME: only works for rectangular IP elements
            // FIXME: only works for forward, backward, central, not for higher as 5-point FD

            // get the IP element geometry interpolator
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

            // IP element max/min
            real tMaxIP = max( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get maximum values of coordinates of IP nodes
            real tMinIP = min( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );    // get minimum values of coordinates of IP nodes

            // get maximum possible perturbation
            real tMaxPerturb = ( tMaxIP - tMinIP ) / 3.0;

            // compute the perturbation value
            real tDeltaH = build_perturbation_size(
                    aPerturbation,
                    aCoefficientToPerturb,
                    tMaxPerturb,
                    mToleranceFD );

            // check that IG node coordinate is consistent with minimum and maximum IP coordinates
            MORIS_ASSERT(
                    tMaxIP >= aCoefficientToPerturb - tDeltaH &&    //
                            tMinIP <= aCoefficientToPerturb + tDeltaH,
                    "ERROR: IG coordinates are outside IP element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  \n",
                    aSpatialDirection,
                    tMinIP,
                    tMaxIP,
                    aCoefficientToPerturb );

            // check point location
            if ( aCoefficientToPerturb + tDeltaH >= tMaxIP )
            {
                aUsedFDScheme = fem::FDScheme_Type::POINT_1_BACKWARD;

                // check for correctness of perturbation size for backward FD
                MORIS_ASSERT( tDeltaH < aCoefficientToPerturb - tMinIP,
                        "ERROR: backward perturbation size exceed limits of interpolation element:\n"
                        "dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                        aSpatialDirection,
                        tMinIP,
                        tMaxIP,
                        aCoefficientToPerturb,
                        tMaxPerturb,
                        tDeltaH,
                        aPerturbation );
            }
            else
            {
                if ( aCoefficientToPerturb - tDeltaH <= tMinIP )
                {
                    aUsedFDScheme = fem::FDScheme_Type::POINT_1_FORWARD;

                    // check for correctness of perturbation size for backward FD
                    MORIS_ASSERT( tDeltaH < tMaxIP - aCoefficientToPerturb,
                            "ERROR: forward perturbation size exceeds limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
                else
                {
                    // check for correctness of perturbation size for central FD
                    MORIS_ASSERT(
                            tDeltaH < tMaxIP - aCoefficientToPerturb &&    //
                                    tDeltaH < aCoefficientToPerturb - tMinIP,
                            "ERROR: central perturbation size exceed limits of interpolation element: dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
                            aSpatialDirection,
                            tMinIP,
                            tMaxIP,
                            aCoefficientToPerturb,
                            tMaxPerturb,
                            tDeltaH,
                            aPerturbation );
                }
            }

            return tDeltaH;
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_bulk(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // get number of leader GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // coefficients for dv type wrt which derivative is computed
            Matrix< DDRMat > tCoeff      = tIGGI->get_space_coeff();
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );
            real tGPWeight = aWStar / tIGGI->det_J();

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );

                            Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();

                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) *                    //
                                    mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset the coefficients values
                tIGGI->set_space_coeff( tCoeff );
                tIGGI->set_space_param_coeff( tParamCoeff );
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::select_dQIdp_FD_geometry_sideset(
                moris::real                        aWStar,
                moris::real                        aPerturbation,
                fem::FDScheme_Type                 aFDSchemeType,
                Matrix< DDSMat >&                  aGeoLocalAssembly,
                Vector< Matrix< IndexMat > >& aVertexIndices )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // get the GI for the IP and IG element considered
            Geometry_Interpolator* tIPGI =
                    mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
            Geometry_Interpolator* tIGGI =
                    mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();

            // store unperturbed xyz coordinates
            Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();

            // store unperturbed local coordinates
            Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();

            // store unperturbed evaluation point
            Matrix< DDRMat > tEvaluationPoint;
            tIGGI->get_space_time( tEvaluationPoint );

            // store unperturbed evaluation point weight
            real tGPWeight = aWStar / tIGGI->det_J();

            // store unperturbed normal
            Matrix< DDRMat > tNormal;
            tIGGI->get_normal( tNormal );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // init perturbation
            real tDeltaH = 0.0;

            // get number of leader GI bases and space dimensions
            uint tDerNumBases      = tIGGI->get_number_of_space_bases();
            uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

            // loop over the spatial directions/loop on pdv type
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
            {
                // loop over the IG nodes/loop one nodes
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // get the geometry pdv assembly index
                    sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                    // if pdv is active
                    if ( tPdvAssemblyIndex != -1 )
                    {
                        // check point location and define perturbation size and FD scheme accordingly
                        fem::FDScheme_Type tUsedFDScheme = aFDSchemeType;
                        tDeltaH                          = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tCoeff( iCoeffRow, iCoeffCol ),
                                iCoeffCol,
                                tUsedFDScheme );

                        // finalize FD scheme
                        fd_scheme( tUsedFDScheme, tFDScheme );
                        uint tNumFDPoints = tFDScheme( 0 ).size();

                        // set starting point for FD
                        uint tStartPoint = 0;

                        // if backward or forward fd
                        if ( ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                                ( tUsedFDScheme == fem::FDScheme_Type::POINT_1_FORWARD ) )
                        {
                            // add unperturbed QI contribution to dQIdp
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( 0 ) * tQI( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // skip first point in FD
                            tStartPoint = 1;
                        }

                        // loop over point of FD scheme
                        for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                        {
                            // reset the perturbed coefficients
                            Matrix< DDRMat > tCoeffPert = tCoeff;

                            // perturb the coefficient
                            tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                            // setting the perturbed coefficients
                            tIGGI->set_space_coeff( tCoeffPert );

                            // update local coordinates
                            Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                            Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                            tIPGI->update_local_coordinates( tXCoords, tXiCoords );

                            Matrix< DDRMat > tParamCoeffPert = tParamCoeff;

                            tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();

                            tIGGI->set_space_param_coeff( tParamCoeffPert );

                            // set evaluation point for interpolators (FIs and GIs)
                            mSet->get_field_interpolator_manager()->    //
                                    set_space_time_from_local_IG_point( tEvaluationPoint );

                            // reset the normal
                            Matrix< DDRMat > tNormalPert;
                            tIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );

                            // reset properties, CM and SP for IWG
                            this->reset_eval_flags();

                            // reset the QI
                            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                            // compute the QI
                            real tWStarPert = tGPWeight * tIGGI->det_J();
                            this->compute_QI( tWStarPert );

                            // evaluate dQIdpGeo
                            mSet->get_dqidpgeo()( tIQIAssemblyIndex )( tPdvAssemblyIndex ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_QI()( tIQIAssemblyIndex )( 0 ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                }
                // reset xyz values
                tIGGI->set_space_coeff( tCoeff );

                // reset local coordinates
                tIGGI->set_space_param_coeff( tParamCoeff );

                // reset evaluation point
                mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                // reset normal
                this->set_normal( tNormal );
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // if active cluster measure on IQI
            if ( mActiveCMEAFlag )
            {
                // add their contribution to dQIdp
                this->add_cluster_measure_dQIdp_FD_geometry(
                        aWStar,
                        aPerturbation,
                        aFDSchemeType );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI::add_cluster_measure_dQIdp_FD_geometry(
                moris::real        aWStar,
                moris::real        aPerturbation,
                fem::FDScheme_Type aFDSchemeType )
        {
            // get the IQI index
            uint tIQIAssemblyIndex = mSet->get_QI_assembly_index( mName );

            // store QI value
            Matrix< DDRMat > tQIStore = mSet->get_QI()( tIQIAssemblyIndex );

            // reset the QI
            mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

            // compute the QI
            this->compute_QI( aWStar );

            // store QI value
            Matrix< DDRMat > tQI = mSet->get_QI()( tIQIAssemblyIndex );

            // init FD scheme
            Vector< Vector< real > > tFDScheme;

            // initialize perturbation of cluster measure
            real tDeltaCM = 0.0;

            // loop over the cluster measures
            for ( uint iCMEA = 0; iCMEA < mCluster->get_cluster_measures().size(); iCMEA++ )
            {
                // get treated cluster measure
                std::shared_ptr< Cluster_Measure >& tClusterMeasure =
                        mCluster->get_cluster_measures()( iCMEA );

                // evaluate the perturbation of cluster measure
                tDeltaCM = this->build_perturbation_size(
                        aPerturbation,
                        tClusterMeasure->val()( 0 ),
                        std::max( tClusterMeasure->val()( 0 ), mToleranceFD ),
                        mToleranceFD );

                // create FD scheme
                fd_scheme( aFDSchemeType, tFDScheme );
                uint tNumFDPoints = tFDScheme( 0 ).size();

                // set starting point for FD
                uint tStartPoint = 0;

                // if backward or forward fd
                if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                        ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                {
                    // add unperturbed QI contribution to dQIdp
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( 0 ) * tQI( 0 ) *    //
                            tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                    // skip first point in FD
                    tStartPoint = 1;
                }

                // loop over point of FD scheme
                for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                {
                    // perturb the cluster measure
                    tClusterMeasure->perturb_cluster_measure( tFDScheme( 0 )( iPoint ) * tDeltaCM );

                    // reset properties, CM and SP for IWG
                    this->reset_eval_flags();

                    // reset the QI
                    mSet->get_QI()( tIQIAssemblyIndex ).fill( 0.0 );

                    // compute the QI
                    this->compute_QI( aWStar );

                    // evaluate dQIdpGeo
                    mSet->get_dqidpgeo()( tIQIAssemblyIndex ) +=
                            tFDScheme( 1 )( iPoint ) *                                                  //
                            mSet->get_QI()( tIQIAssemblyIndex )( 0 ) * tClusterMeasure->dMEAdPDV() /    //
                            ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                    // reset cluster measures
                    mCluster->reset_cluster_measure();
                }
            }

            // reset QI value
            mSet->get_QI()( tIQIAssemblyIndex ) = tQIStore;

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_dqidpgeo()( tIQIAssemblyIndex ) ),
                    "IQI::compute_dQIdp_FD_geometry - dQIdp contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
