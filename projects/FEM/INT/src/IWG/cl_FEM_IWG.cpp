/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG.cpp
 *
 */

#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Cluster.hpp"
#include "cl_FEM_Cluster_Measure.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Model.hpp"
#include "cl_Vector.hpp"
#include "fn_max.hpp"
#include "fn_min.hpp"
#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_vectorize.hpp"
#include "fn_join_horiz.hpp"
#include "fn_unique.hpp"

#include <cl_MTK_Ray_Line_Intersection.hpp>
#include <iomanip>
#include <utility>

namespace moris::fem
{
    //------------------------------------------------------------------------------

    void
    IWG::print_names()
    {
        std::cout << "----------" << '\n';
        std::cout << "IWG: " << mName << '\n';

        // properties
        for ( uint iProp = 0; iProp < mLeaderProp.size(); iProp++ )
        {
            if ( mLeaderProp( iProp ) != nullptr )
            {
                std::cout << "Leader property: " << mLeaderProp( iProp )->get_name() << '\n';
            }
        }
        for ( uint iProp = 0; iProp < mFollowerProp.size(); iProp++ )
        {
            if ( mFollowerProp( iProp ) != nullptr )
            {
                std::cout << "Follower property:  " << mFollowerProp( iProp )->get_name() << '\n';
            }
        }

        // CM
        for ( uint iCM = 0; iCM < mLeaderCM.size(); iCM++ )
        {
            if ( mLeaderCM( iCM ) != nullptr )
            {
                std::cout << "Leader CM:       " << mLeaderCM( iCM )->get_name() << '\n';
            }
        }
        for ( uint iCM = 0; iCM < mFollowerCM.size(); iCM++ )
        {
            if ( mFollowerCM( iCM ) != nullptr )
            {
                std::cout << "Follower CM:        " << mFollowerCM( iCM )->get_name() << '\n';
            }
        }

        // SP
        for ( uint iSP = 0; iSP < mStabilizationParam.size(); iSP++ )
        {
            if ( mStabilizationParam( iSP ) != nullptr )
            {
                std::cout << "SP:              " << mStabilizationParam( iSP )->get_name() << '\n';
            }
        }
        std::cout << "----------" << '\n';
    }

    //------------------------------------------------------------------------------

    void
    IWG::reset_eval_flags()
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

        // reset material models
        for ( const std::shared_ptr< Material_Model >& tMM : mLeaderMM )
        {
            if ( tMM != nullptr )
            {
                tMM->reset_eval_flags();
            }
        }
        for ( const std::shared_ptr< Material_Model >& tMM : mFollowerMM )
        {
            if ( tMM != nullptr )
            {
                tMM->reset_eval_flags();
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

        // reset evaluation flags specific to child implementations
        this->reset_spec_eval_flags();

        // reset gap data if it exists
        if ( mGapData != nullptr )
        {
            mGapData->mEval = true;
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_function_pointers()
    {
        // switch on element type
        switch ( mSet->get_element_type() )
        {
            case fem::Element_Type::BULK:
            {
                m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_bulk;
                break;
            }
            case fem::Element_Type::SIDESET:
            {
                m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_sideset;
                break;
            }
            case fem::Element_Type::TIME_SIDESET:
            case fem::Element_Type::TIME_FINAL_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
            {
                m_compute_jacobian_FD      = &IWG::select_jacobian_FD;
                m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material;
                m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_time_sideset;
                break;
            }
            case fem::Element_Type::DOUBLE_SIDESET:
            case fem::Element_Type::NONCONFORMAL_SIDESET:
            {
                m_compute_jacobian_FD      = &IWG::select_jacobian_FD_double;
                m_compute_dRdp_FD_material = &IWG::select_dRdp_FD_material_double;
                m_compute_dRdp_FD_geometry = &IWG::select_dRdp_FD_geometry_double;
                break;
            }
            default:
                MORIS_ERROR( false, "IWG::set_function_pointers - unknown element type." );
        }

        // switch on perturbation strategy
        switch ( mSet->get_perturbation_strategy() )
        {
            case fem::Perturbation_Type::RELATIVE:
            {
                m_build_perturbation_size = &IWG::build_perturbation_size_relative;
                break;
            }
            case fem::Perturbation_Type::ABSOLUTE:
            {
                m_build_perturbation_size = &IWG::build_perturbation_size_absolute;
                break;
            }
            default:
                MORIS_ERROR( false, "IWG::set_function_pointers - unknown perturbation type." );
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_phase_name(
            const std::string&   aPhaseName,
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
                MORIS_ERROR( false, "IWG::set_phase_name - aIsLeader can only be leader or follower." );
            }
        }
    }

    //------------------------------------------------------------------------------

    std::string
    IWG::get_phase_name( mtk::Leader_Follower aIsLeader )
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
                MORIS_ERROR( false, "IWG::get_phase_name - aIsLeader can only be leader or follower." );
                return mLeaderPhaseName;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_field_interpolator_manager(
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
                MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be leader or follower" );
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

        // loop over the material models
        for ( const std::shared_ptr< Material_Model >& tMM : this->get_material_models( aIsLeader ) )
        {
            if ( tMM != nullptr )
            {
                // set the field interpolator manager for the CM
                tMM->set_field_interpolator_manager( this->get_field_interpolator_manager( aIsLeader ) );

                // set the fem set pointe for the CM
                tMM->set_set_pointer( mSet );
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
    IWG::get_field_interpolator_manager(
            mtk::Leader_Follower aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
                return mLeaderFIManager;

            case mtk::Leader_Follower::FOLLOWER:
                return mFollowerFIManager;

            default:
                MORIS_ERROR( false, "IWG::get_field_interpolator_manager - can only be leader or follower." );
                return mLeaderFIManager;
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_field_interpolator_manager_previous_time(
            Field_Interpolator_Manager* aFieldInterpolatorManager,
            mtk::Leader_Follower        aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                mLeaderPreviousFIManager = aFieldInterpolatorManager;
                break;
            }

            default:
            {
                MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be leader" );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_normal( const Matrix< DDRMat >& aNormal )
    {
        mNormal = aNormal;

        // set normal for SP
        for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
        {
            if ( tSP != nullptr )
            {
                tSP->set_normal( mNormal );
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_interpolation_order()
    {
        // if order is already set
        if ( ( mOrder != MORIS_UINT_MAX ) ) return;

        // get field interpolator manager
        Field_Interpolator_Manager* tFIManager = mSet->get_field_interpolator_manager();

        // define axuxilliary variable for storing interpolation order
        uint tPrevOrder = MORIS_UINT_MAX;

        // loop overall residual dof types
        for ( uint iType = 0; iType < mResidualDofType.size(); ++iType )
        {
            // get field interpolator for current residual dof type;
            Field_Interpolator* tFI = tFIManager->get_field_interpolators_for_type( mResidualDofType( iType )( 0 ) );

            // get residual dof type interpolation order
            mtk::Interpolation_Order tInterpOrder = tFI->get_space_interpolation_order();

            // set the interpolation order for IWG
            switch ( tInterpOrder )
            {
                case mtk::Interpolation_Order::LINEAR:
                {
                    mOrder = 1;
                    break;
                }
                case mtk::Interpolation_Order::QUADRATIC:
                {
                    mOrder = 2;
                    break;
                }
                case mtk::Interpolation_Order::CUBIC:
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG::set_interpolation_order - order not supported" );
                }
            }

            // check order is the same for all residual types
            MORIS_ERROR( iType > 0 ? tPrevOrder == mOrder : true,
                    "IWG::set_interpolation_order - Interpolation order of residual fields need to be identical among all residual fields.\n" );

            tPrevOrder = mOrder;
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_dof_type_list(
            const Vector< Vector< MSI::Dof_Type > >& aDofTypes,
            mtk::Leader_Follower                     aIsLeader )
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
                MORIS_ERROR( false, "IWG::set_dof_type_list - can only be LEADER or FOLLOWER." );
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< MSI::Dof_Type > >&
    IWG::get_dof_type_list(
            mtk::Leader_Follower aIsLeader ) const
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
                MORIS_ASSERT( false, "IWG::get_dof_type_list - can only be leader or follower." );
                return mLeaderDofTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_dv_type_list(
            const Vector< Vector< gen::PDV_Type > >& aDvTypes,
            mtk::Leader_Follower                     aIsLeader )
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
                MORIS_ERROR( false, "IWG::set_dv_type_list - can only be LEADER or FOLLOWER." );
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< gen::PDV_Type > >&
    IWG::get_dv_type_list(
            mtk::Leader_Follower aIsLeader ) const
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader global dof type list
                return mLeaderDvTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerDvTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "IWG::get_dv_type_list - can only be leader or follower." );
                return mLeaderDvTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_field_type_list(
            const Vector< Vector< mtk::Field_Type > >& aDofTypes,
            mtk::Leader_Follower                       aIsLeader )
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

    const Vector< Vector< mtk::Field_Type > >&
    IWG::get_field_type_list(
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
    IWG::set_property(
            std::shared_ptr< Property > aProperty,
            const std::string&          aPropertyString,
            mtk::Leader_Follower        aIsLeader )
    {
        // check that aPropertyString makes sense
        MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                "IWG::set_property - IWG %s - Unknown aPropertyString: %s ",
                mName.c_str(),
                aPropertyString.c_str() );

        // set the property in the property pointer cell
        this->get_properties( aIsLeader )( mPropertyMap[ aPropertyString ] ) = std::move( aProperty );
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< Property > >&
    IWG::get_properties(
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
                MORIS_ASSERT( false, "IWG::get_properties - can only be leader or follower." );
                return mLeaderProp;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_material_model(
            std::shared_ptr< Material_Model > aMaterialModel,
            const std::string&                aMaterialModelString,
            mtk::Leader_Follower              aIsLeader )
    {
        // check that aConstitutiveString makes sense
        MORIS_ERROR( mMaterialMap.find( aMaterialModelString ) != mMaterialMap.end(),
                "IWG::set_material_model - IWG %s - Unknown aMaterialModelString: %s ",
                mName.c_str(),
                aMaterialModelString.c_str() );

        // set the CM in the CM pointer cell
        this->get_material_models( aIsLeader )( mMaterialMap[ aMaterialModelString ] ) = std::move( aMaterialModel );
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< Material_Model > >&
    IWG::get_material_models(
            mtk::Leader_Follower aIsLeader )
    {
        // switch on leader/follower
        switch ( aIsLeader )
        {
            // if leader
            case mtk::Leader_Follower::LEADER:
            {
                // return leader property pointers
                return mLeaderMM;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower property pointers
                return mFollowerMM;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "IWG::get_material_models - can only be leader or follower." );
                return mLeaderMM;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_constitutive_model(
            std::shared_ptr< Constitutive_Model > aConstitutiveModel,
            const std::string&                    aConstitutiveString,
            mtk::Leader_Follower                  aIsLeader )
    {
        // check that aConstitutiveString makes sense
        MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                "IWG::set_constitutive_model - IWG %s - Unknown aConstitutiveString: %s ",
                mName.c_str(),
                aConstitutiveString.c_str() );

        // set the CM in the CM pointer cell
        this->get_constitutive_models( aIsLeader )( mConstitutiveMap[ aConstitutiveString ] ) = std::move( aConstitutiveModel );
    }

    //------------------------------------------------------------------------------

    Vector< std::shared_ptr< Constitutive_Model > >&
    IWG::get_constitutive_models(
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
                MORIS_ASSERT( false, "IWG::get_constitutive_models - can only be leader or follower." );
                return mLeaderCM;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::set_stabilization_parameter(
            const std::shared_ptr< Stabilization_Parameter >& aStabilizationParameter,
            const std::string&                                aStabilizationString )
    {
        // check that aStabilizationString makes sense
        MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(),
                "IWG::set_stabilization_parameter - IWG %s - Unknown aStabilizationString: %s ",
                mName.c_str(),
                aStabilizationString.c_str() );

        // set the stabilization parameter in the stabilization parameter cell
        this->get_stabilization_parameters()( mStabilizationMap[ aStabilizationString ] ) = aStabilizationParameter;

        // set active cluster measure on IWG flag on/off
        mActiveCMEAFlag =
                mActiveCMEAFlag || ( aStabilizationParameter->get_cluster_measure_tuple_list().size() > 0 );
    }

    //------------------------------------------------------------------------------

    void
    IWG::get_non_unique_dof_dv_and_field_types(
            Vector< Vector< MSI::Dof_Type > >&   aDofTypes,
            Vector< Vector< gen::PDV_Type > >&   aDvTypes,
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
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tProperty->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update dof and dv counters
                tLeaderDofCounter += tActiveDofTypes.size();
                tLeaderDvCounter += tActiveDvTypes.size();
                tLeaderFieldCounter += tActiveDvTypes.size();
            }
        }

        // loop over follower properties
        for ( const std::shared_ptr< Property >& tProperty : mFollowerProp )
        {
            if ( tProperty != nullptr )
            {
                // get property non unique dof and dv type lists
                // get property non unique dof and dv type list
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
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

        // loop over leader material models
        for ( const std::shared_ptr< Material_Model >& tMM : mLeaderMM )
        {
            if ( tMM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type > tActiveDofTypes;
                Vector< gen::PDV_Type > tActiveDvTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // update dof counters (DVs not part of MM yet)
                tLeaderDofCounter += tActiveDofTypes.size();
            }
        }

        // loop over follower material models
        for ( const std::shared_ptr< Material_Model >& tMM : mFollowerMM )
        {
            if ( tMM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type > tActiveDofTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // update dof and dv counters
                tFollowerDofCounter += tActiveDofTypes.size();
            }
        }

        // loop over leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
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
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tCM->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // update dof and dv counters
                tFollowerDofCounter += tActiveDofTypes.size();
                tFollowerDvCounter += tActiveDvTypes.size();
            }
        }

        // loop over leader stabilization parameters
        for ( const std::shared_ptr< Stabilization_Parameter >& tSP : mStabilizationParam )
        {
            if ( tSP != nullptr )
            {
                // get SP non unique dof type list
                Vector< MSI::Dof_Type > tActiveDofTypes;
                Vector< gen::PDV_Type > tActiveDvTypes;

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
            // populate the field list
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
                Vector< gen::PDV_Type >   tActiveDvTypes;
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
                Vector< gen::PDV_Type >   tActiveDvTypes;
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

        // loop over the leader material models
        for ( const std::shared_ptr< Material_Model >& tMM : mLeaderMM )
        {
            if ( tMM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type > tActiveDofTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // update dof counters (DVs not part of MM yet)
                tLeaderDofCounter += tActiveDofTypes.size();
            }
        }

        // loop over the follower material models
        for ( const std::shared_ptr< Material_Model >& tMM : mFollowerMM )
        {
            if ( tMM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type > tActiveDofTypes;

                tMM->get_non_unique_dof_types( tActiveDofTypes );

                // update dof counters (DVs not part of MM yet)
                tLeaderDofCounter += tActiveDofTypes.size();
            }
        }

        // loop over the leader constitutive models
        for ( const std::shared_ptr< Constitutive_Model >& tCM : mLeaderCM )
        {
            if ( tCM != nullptr )
            {
                // get CM non unique dof and dv type lists
                Vector< MSI::Dof_Type >   tActiveDofTypes;
                Vector< gen::PDV_Type >   tActiveDvTypes;
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
                Vector< gen::PDV_Type >   tActiveDvTypes;
                Vector< mtk::Field_Type > tActiveFieldTypes;

                tCM->get_non_unique_dof_dv_and_field_types(
                        tActiveDofTypes,
                        tActiveDvTypes,
                        tActiveFieldTypes );

                // populate the dof and dv lists
                aDofTypes( 1 ).append( tActiveDofTypes );
                aDvTypes( 1 ).append( tActiveDvTypes );
                aFieldTypes( 0 ).append( tActiveFieldTypes );
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
                Vector< gen::PDV_Type > tActiveDvTypes;

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
    IWG::build_global_dof_dv_and_field_type_list()
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
            sint tDofTypeIndex =
                    mSet->get_index_from_unique_dof_type_map( mLeaderDofTypes( iDof )( 0 ) );    // FIXME'

            // put the dof type in the checklist
            tDofCheckList( tDofTypeIndex ) = 1;

            // put the dof type in the global type list
            mLeaderGlobalDofTypes.push_back( mLeaderDofTypes( iDof ) );
        }

        // get dv type from direct dependencies
        for ( uint iDv = 0; iDv < mLeaderDvTypes.size(); iDv++ )
        {
            // get set index for dv type
            sint tDvTypeIndex =
                    mSet->get_index_from_unique_dv_type_map( mLeaderDvTypes( iDv )( 0 ) );    // FIXME'

            // put the dv type in the checklist
            tDvCheckList( tDvTypeIndex ) = 1;

            // put the dv type in the global type list
            mLeaderGlobalDvTypes.push_back( mLeaderDvTypes( iDv ) );
        }

        // get field type from direct dependencies
        for ( uint iFi = 0; iFi < mLeaderFieldTypes.size(); iFi++ )
        {
            // get set index for field type
            sint tFieldTypeIndex =
                    mSet->get_index_from_unique_field_type_map( mLeaderFieldTypes( iFi )( 0 ) );    // FIXME'

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
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for property
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tProperty->get_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex =
                            mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dof enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
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
                    sint tFieldTypeIndex =
                            mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

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

        // get dof type from leader material models
        for ( const std::shared_ptr< Material_Model >& tMM : mLeaderMM )
        {
            if ( tMM != nullptr )
            {
                // get dof types for material modIWG::build_global_dof_dv_and_field_listel
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                        tMM->get_global_dof_type_list();

                // loop on property dof type
                for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                {
                    // get set index for dof type
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }
                // skip loop on material model dv type - not implemented
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
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for constitutive model
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tCM->get_global_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex =
                            mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dv type in the check list
                        tDvCheckList( tDvTypeIndex ) = 1;

                        // put the dv type in the global type list
                        mLeaderGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                    }
                }

                // get field types for constitutive model
                const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                        tCM->get_global_field_type_list();

                // loop on property field type
                for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                {
                    // get set index for field type
                    sint tFieldTypeIndex =
                            mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

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
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mLeaderGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for constitutive model
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tSP->get_global_dv_type_list( mtk::Leader_Follower::LEADER );

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex =
                            mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dv type in the check list
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
            sint tDofTypeIndex =
                    mSet->get_index_from_unique_dof_type_map( mFollowerDofTypes( iDof )( 0 ) );

            // put the dof type in the check list
            tDofCheckList( tDofTypeIndex ) = 1;

            // put the dof type in the global type list
            mFollowerGlobalDofTypes.push_back( mFollowerDofTypes( iDof ) );
        }

        // get dv type from follower direct dependencies
        for ( uint iDv = 0; iDv < mFollowerDvTypes.size(); iDv++ )
        {
            // get set index for dv type
            sint tDvTypeIndex =
                    mSet->get_index_from_unique_dv_type_map( mFollowerDvTypes( iDv )( 0 ) );

            // put the dv type in the check list
            tDvCheckList( tDvTypeIndex ) = 1;

            // put the dv type in the global type list
            mFollowerGlobalDvTypes.push_back( mFollowerDvTypes( iDv ) );
        }

        // get field type from follower direct dependencies
        for ( uint iFi = 0; iFi < mFollowerFieldTypes.size(); iFi++ )
        {
            // get set index for field type
            sint tFieldTypeIndex =
                    mSet->get_index_from_unique_field_type_map( mFollowerFieldTypes( iFi )( 0 ) );

            // put the field type in the check list
            tFieldCheckList( tFieldTypeIndex ) = 1;

            // put the field type in the global type list
            mFollowerGlobalFieldTypes.push_back( mFollowerFieldTypes( iFi ) );
        }

        // get dof type from follower properties
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
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for property
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tProperty->get_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex =
                            mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dv type in the check list
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
                    sint tFieldTypeIndex =
                            mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

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

        // get dof type from follower material models
        for ( const std::shared_ptr< Material_Model >& tMM : mFollowerMM )
        {
            if ( tMM != nullptr )
            {
                // get dof types for constitutive model
                const Vector< Vector< MSI::Dof_Type > >& tActiveDofTypes =
                        tMM->get_global_dof_type_list();

                // loop on property dof type
                for ( uint iDof = 0; iDof < tActiveDofTypes.size(); iDof++ )
                {
                    // get set index for dof type
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }
                // skip loop on material dv type - not implemented
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
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for constitutive model
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tCM->get_global_dv_type_list();

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dv type in the check list
                        tDvCheckList( tDvTypeIndex ) = 1;

                        // put the dv type in the global type list
                        mFollowerGlobalDvTypes.push_back( tActiveDvTypes( iDv ) );
                    }
                }

                // get field types for constitutive model
                const Vector< Vector< mtk::Field_Type > >& tActiveFieldTypes =
                        tCM->get_field_type_list();

                // loop on constitutive model field type
                for ( uint iFi = 0; iFi < tActiveFieldTypes.size(); iFi++ )
                {
                    // get set index for field type
                    sint tFieldTypeIndex =
                            mSet->get_index_from_unique_field_type_map( tActiveFieldTypes( iFi )( 0 ) );

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
                    sint tDofTypeIndex =
                            mSet->get_index_from_unique_dof_type_map( tActiveDofTypes( iDof )( 0 ) );

                    // if dof enum not in the list
                    if ( tDofCheckList( tDofTypeIndex ) != 1 )
                    {
                        // put the dof type in the check list
                        tDofCheckList( tDofTypeIndex ) = 1;

                        // put the dof type in the global type list
                        mFollowerGlobalDofTypes.push_back( tActiveDofTypes( iDof ) );
                    }
                }

                // get dv types for stabilization parameter
                const Vector< Vector< gen::PDV_Type > >& tActiveDvTypes =
                        tSP->get_global_dv_type_list( mtk::Leader_Follower::FOLLOWER );

                // loop on property dv type
                for ( uint iDv = 0; iDv < tActiveDvTypes.size(); iDv++ )
                {
                    // get set index for dv type
                    sint tDvTypeIndex = mSet->get_index_from_unique_dv_type_map( tActiveDvTypes( iDv )( 0 ) );

                    // if dv enum not in the list
                    if ( tDvCheckList( tDvTypeIndex ) != 1 )
                    {
                        // put the dv type in the check list
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
    IWG::build_requested_dof_type_list( const bool aIsStaggered )
    {
        // clear the dof lists
        mRequestedLeaderGlobalDofTypes.clear();
        mRequestedFollowerGlobalDofTypes.clear();

        Vector< enum MSI::Dof_Type > tRequestedDofTypes;

        // if residual evaluation
        if ( aIsStaggered )
        {
            // get the requested dof types
            tRequestedDofTypes = mSet->get_secondary_dof_types();
        }
        // if Jacobian evaluation
        else
        {
            // get the requested dof types
            tRequestedDofTypes = mSet->get_requested_dof_types();
        }

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
    IWG::check_field_interpolators( mtk::Leader_Follower aIsLeader )
    {
        switch ( aIsLeader )
        {
            case mtk::Leader_Follower::LEADER:
            {
                // loop over the dof field interpolator pointers
                for ( uint iDofFI = 0; iDofFI < mRequestedLeaderGlobalDofTypes.size(); iDofFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT(
                            this->get_field_interpolator_manager( aIsLeader )->    //
                                    get_field_interpolators_for_type( mRequestedLeaderGlobalDofTypes( iDofFI )( 0 ) )
                                    != nullptr,
                            "IWG::check_field_interpolators - Leader dof FI missing. " );
                }

                // loop over the dv field interpolator pointers
                for ( uint iDvFI = 0; iDvFI < mLeaderGlobalDvTypes.size(); iDvFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT(
                            this->get_field_interpolator_manager( aIsLeader )->    //
                                    get_field_interpolators_for_type( mLeaderGlobalDvTypes( iDvFI )( 0 ) )
                                    != nullptr,
                            "IWG::check_field_interpolators - Leader dv FI missing. " );
                }
                break;
            }
            case mtk::Leader_Follower::FOLLOWER:
            {
                // loop over the dof field interpolator pointers
                for ( uint iDofFI = 0; iDofFI < mRequestedFollowerGlobalDofTypes.size(); iDofFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT(
                            this->get_field_interpolator_manager( aIsLeader )->    //
                                    get_field_interpolators_for_type( mRequestedFollowerGlobalDofTypes( iDofFI )( 0 ) )
                                    != nullptr,
                            "IWG::check_dof_field_interpolators - Follower dof FI missing. " );
                }

                // loop over the dv field interpolator pointers
                for ( uint iDvFI = 0; iDvFI < mFollowerGlobalDvTypes.size(); iDvFI++ )
                {
                    // check that the field interpolator was set
                    MORIS_ASSERT(
                            this->get_field_interpolator_manager( aIsLeader )->    //
                                    get_field_interpolators_for_type( mFollowerGlobalDvTypes( iDvFI )( 0 ) )
                                    != nullptr,
                            "IWG::check_field_interpolators - Follower dv FI missing. " );
                }
                break;
            }
            default:
            {
                MORIS_ERROR( false, "IWG::check_field_interpolators - can only be leader or follower." );
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< MSI::Dof_Type > >&
    IWG::get_global_dof_type_list(
            mtk::Leader_Follower aIsLeader )
    {
        // if the global list was not yet built
        if ( mGlobalDofBuild )
        {
            // build the stabilization parameter global dof type list
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
                break;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerGlobalDofTypes;
                break;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "IWG::get_global_dof_type_list - can only be leader or follower." );
                return mLeaderGlobalDofTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< gen::PDV_Type > >&
    IWG::get_global_dv_type_list(
            mtk::Leader_Follower aIsLeader )
    {
        // if the global list was not yet built
        if ( mGlobalDvBuild )
        {
            // build the stabilization parameter global dof type list
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
                return mLeaderGlobalDvTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global dof type list
                return mFollowerGlobalDvTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "IWG::get_global_dv_type_list - can only be leader or follower." );
                return mLeaderGlobalDvTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    const Vector< Vector< mtk::Field_Type > >&
    IWG::get_global_field_type_list(
            mtk::Leader_Follower aIsLeader )
    {
        // if the global list was not yet built
        if ( mGlobalFieldBuild )
        {
            // build the stabilization parameter global dof type list
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
                // return leader global field type list
                return mLeaderGlobalFieldTypes;
            }
            // if follower
            case mtk::Leader_Follower::FOLLOWER:
            {
                // return follower global field type list
                return mFollowerGlobalFieldTypes;
            }
            // if none
            default:
            {
                MORIS_ASSERT( false, "IWG::get_global_field_type_list - can only be leader or follower." );
                return mLeaderGlobalFieldTypes;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::select_jacobian_FD(
            real               aWStar,
            real               aPerturbation,
            fem::FDScheme_Type aFDSchemeType,
            bool               aUseAbsolutePerturbations )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // get leader index for residual dof type, indices for assembly
        sint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get leader number of dof types
        uint tLeaderNumDofTypes = mRequestedLeaderGlobalDofTypes.size();

        // reset and evaluate the residual plus
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual =
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } );

        // loop over the IWG dof types
        for ( uint iFI = 0; iFI < tLeaderNumDofTypes; iFI++ )
        {
            // init dof counter
            uint tDofCounter = 0;

            // get the dof type
            Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iFI );

            // get the index for the dof type
            sint tLeaderDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tLeaderDepDofIndex, 0 );

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
                    real tDeltaH = aPerturbation;
                    if ( !aUseAbsolutePerturbations )
                    {
                        tDeltaH = tDeltaH * tCoeff( iCoeffRow, iCoeffCol );
                    }

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed residual contribution to dRdp
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // assemble the Jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;
    }

    //------------------------------------------------------------------------------

    void
    IWG::select_jacobian_FD_double(
            real               aWStar,
            real               aPerturbation,
            fem::FDScheme_Type aFDSchemeType,
            bool               aUseAbsolutePerturbations )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // if this is a nonconformal sideset, some things will be handled differently!
        //        bool const tIsNonconformal = ( mSet->get_element_type() == Element_Type::NONCONFORMAL_SIDESET );

        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // get leader index for residual dof type, indices for assembly
        sint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get follower index for residual dof type, indices for assembly
        sint tFollowerDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        uint tFollowerResStartIndex = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
        uint tFollowerResStopIndex  = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );

        Geometry_Interpolator* tLeaderIGGI   = this->mLeaderFIManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = this->mFollowerFIManager->get_IG_geometry_interpolator();

        // reset and evaluate the residual plus
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tLeaderResidual =
                mSet->get_residual()( 0 )(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { 0, 0 } );
        Matrix< DDRMat > tFollowerResidual =
                mSet->get_residual()( 0 )(
                        { tFollowerResStartIndex, tFollowerResStopIndex },
                        { 0, 0 } );

        // get leader number of dof types
        uint tLeaderNumDofTypes = mRequestedLeaderGlobalDofTypes.size();

        // loop over the IWG dof types
        for ( uint iFI = 0; iFI < tLeaderNumDofTypes; iFI++ )
        {
            // init dof counter
            uint tDofCounter = 0;

            // get the dof type
            Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iFI );

            // bool const tIsDisplacementField = ( tDofType( 0 ) == MSI::Dof_Type::UX );

            // get the index for the dof type
            sint tLeaderDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tLeaderDepDofIndex, 0 );

            // get field interpolator for dependency dof type
            Field_Interpolator* tLeaderFI = mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

            // get number of leader FI bases and fields
            uint const tDerNumBases  = tLeaderFI->get_number_of_space_time_bases();    // coefficients for the interpolation of the field (number of shape functions)
            uint const tDerNumFields = tLeaderFI->get_number_of_fields();              // e.g. UX, UY, UZ

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tFieldCoefficients = tLeaderFI->get_coeff();

            // loop over the coefficient column
            for ( uint iField = 0; iField < tDerNumFields; iField++ )
            {
                // loop over the coefficient row
                for ( uint iBase = 0; iBase < tDerNumBases; iBase++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation;
                    if ( !aUseAbsolutePerturbations )
                    {
                        tDeltaH = tDeltaH * tFieldCoefficients( iBase, iField );
                    }

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( 0 ) * tLeaderResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // add unperturbed follower residual contribution to dRdp
                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( 0 ) * tFollowerResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over the points for FD
                    for ( uint iFDPoint = tStartPoint; iFDPoint < tNumFDPoints; iFDPoint++ )
                    {
                        // reset the perturbed coefficients
                        Matrix< DDRMat > tCoeffPert = tFieldCoefficients;

                        // perturb the coefficient
                        tCoeffPert( iBase, iField ) += tFDScheme( 0 )( iFDPoint ) * tDeltaH;

                        // set the perturbed coefficients to FI
                        tLeaderFI->set_coeff( tCoeffPert );

                        //                        // If this is a nonconformal side set and the dof type is a displacement field, we need to remap the integration point on the follower side!
                        //                        // For this, we need to use the perturbed displacement field to map the leader integration point and the follower element into the current configuration.
                        //                        // To obtain the new normal, we currently limit the implementation to line elements to directly calculate the normal vector from the deformed coordinates of the leader element.
                        //                        if ( tIsNonconformal && tIsDisplacementField )
                        //                        {
                        //                            Field_Interpolator*    tFollowerFI             = mFollowerFIManager->get_field_interpolators_for_type( tDofType( 0 ) );
                        //                            Matrix< DDRMat > const tRemappedFollowerCoords = this->remap_nonconformal_rays( tLeaderFI, tFollowerFI );
                        //
                        //                            MORIS_ERROR( std::abs( tRemappedFollowerCoords( 0 ) ) <= 1.0,
                        //                                    "IWG::select_jacobian_FD_double - remapped follower coordinates are not in [-1,1] range." );
                        //
                        //                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
                        //                        }

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // assemble leader part of the jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( iFDPoint ) *                                                              //
                                mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // assemble follower part of the jacobian
                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( iFDPoint ) *                                                                  //
                                mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tLeaderFI->set_coeff( tFieldCoefficients );
            tLeaderIGGI->reset_eval_flags_coordinates();    // reset computation of leader geometry interpolator
        }

        // get follower number of dof types
        uint tFollowerNumDofTypes = mRequestedFollowerGlobalDofTypes.size();

        // loop over the IWG dof types
        for ( uint iFI = 0; iFI < tFollowerNumDofTypes; iFI++ )
        {
            // init dof counter
            uint tDofCounter = 0;

            // get the dof type
            Vector< MSI::Dof_Type > tDofType = mRequestedFollowerGlobalDofTypes( iFI );
            //            bool const              tIsDisplacementField = ( tDofType( 0 ) == MSI::Dof_Type::UX );

            // get the index for the dof type
            sint tFollowerDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::FOLLOWER );
            uint tFollowerDepStartIndex = mSet->get_jac_dof_assembly_map()( tFollowerDofIndex )( tFollowerDepDofIndex, 0 );

            // get field interpolator for dependency dof type
            Field_Interpolator* tFollowerFI = mFollowerFIManager->get_field_interpolators_for_type( tDofType( 0 ) );

            // get number of leader FI bases and fields
            uint tDerNumBases  = tFollowerFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFollowerFI->get_number_of_fields();

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFollowerFI->get_coeff();

            // loop over the coefficient column
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over the coefficient row
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    real tDeltaH = aPerturbation;
                    if ( !aUseAbsolutePerturbations )
                    {
                        tDeltaH = tDeltaH * tCoeff( iCoeffRow, iCoeffCol );
                    }

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( 0 ) * tLeaderResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // add unperturbed follower residual contribution to dRdp
                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( 0 ) * tFollowerResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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
                        tFollowerFI->set_coeff( tCoeffPert );

                        //                        // If this is a nonconformal side set and the dof type is a displacement field, we need to remap the integration point on the follower side!
                        //                        // For this, we need to use the perturbed displacement field to map the leader integration point and the follower element into the current configuration.
                        //                        // To obtain the new normal, we currently limit the implementation to line elements to directly calculate the normal vector from the deformed coordinates of the leader element.
                        //                        if ( tIsNonconformal && tIsDisplacementField )
                        //                        {
                        //                            Field_Interpolator*    tLeaderFI               = mLeaderFIManager->get_field_interpolators_for_type( tDofType( 0 ) );
                        //                            Matrix< DDRMat > const tRemappedFollowerCoords = this->remap_nonconformal_rays( tLeaderFI, tFollowerFI );
                        //
                        //                            MORIS_ERROR( std::abs( tRemappedFollowerCoords( 0 ) ) <= 1.0,
                        //                                    "IWG::select_jacobian_FD_double - remapped follower coordinates are not in [-1,1] range." );
                        //
                        //                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
                        //                        }

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // assemble the jacobian
                        mSet->get_jacobian()(
                                { tLeaderResStartIndex, tLeaderResStopIndex },
                                { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tLeaderResStartIndex, tLeaderResStopIndex }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // assemble the jacobian
                        mSet->get_jacobian()(
                                { tFollowerResStartIndex, tFollowerResStopIndex },
                                { tFollowerDepStartIndex + tDofCounter, tFollowerDepStartIndex + tDofCounter } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                    //
                                mSet->get_residual()( 0 )( { tFollowerResStartIndex, tFollowerResStopIndex }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update dof counter
                    tDofCounter++;
                }
            }
            // reset the coefficients values
            tFollowerFI->set_coeff( tCoeff );
            tFollowerIGGI->reset_eval_flags_coordinates();    // reset computation of follower geometry
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;
    }

    //------------------------------------------------------------------------------
    Matrix< DDRMat > IWG::remap_nonconformal_rays_undeformed_geometry() const
    {
        // get geometry interpolator for leader and follower IG elements
        Geometry_Interpolator* tLeaderIGGI   = this->mLeaderFIManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = this->mFollowerFIManager->get_IG_geometry_interpolator();

        //        // reset evaluation flags for the geometry interpolators -
        //        // FIXME: not sure whether this is needed
        //        tLeaderIGGI->reset_eval_flags();
        //        tFollowerIGGI->reset_eval_flags();

        // get the integration point location
        const Matrix< DDRMat >& tIGPoint = tLeaderIGGI->valx();

        // get the normal vector at the integration point
        const Matrix< DDRMat >& tNormal = tLeaderIGGI->get_normal();

        // Get the nodes of the follower element
        const Matrix< DDRMat >& tFollowerCoordinates = tFollowerIGGI->get_space_coeff();

        // check that we are dealing with line elements
        MORIS_ASSERT( tFollowerCoordinates.n_cols() == 2 && tNormal.n_rows() == 2,
                "IWG::remap_nonconformal_rays - Nonconformal FD scheme is only implemented for line elements." );

        // Perform the mapping
        mtk::Ray_Line_Intersection tRLI( tFollowerCoordinates.n_cols() );
        tRLI.set_ray_origin( trans( tIGPoint ) );
        tRLI.set_ray_direction( tNormal );
        tRLI.set_target_origin( trans( tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.set_target_span( trans( tFollowerCoordinates.get_row( 1 ) - tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.perform_raytracing();

        Matrix< DDRMat > tFollowerSpaceTime;
        tLeaderIGGI->get_space_time( tFollowerSpaceTime );    // use leader as follower may not be initialized

        if ( tRLI.has_intersection() )
        {
            tFollowerSpaceTime( 0 ) = tRLI.get_intersection_parametric()( 0 );

            // set the space time coordinates of the follower side
            this->mFollowerFIManager->set_space_time_from_local_IG_point( tFollowerSpaceTime );
        }
        else
        {
            tFollowerSpaceTime( 0 ) = -1.0;    // set to -1 if no intersection found
        }

        return tFollowerSpaceTime;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IWG::remap_nonconformal_rays(
            Field_Interpolator* aLeaderFieldInterpolator,
            Field_Interpolator* aFollowerFieldInterpolator )
    {
        if ( mUseDeformedGeometryForGap )
        {
            if ( mUseConsistentDeformedGeometryForGap )
            {
                return this->remap_nonconformal_rays_consistent_deformed_geometry(
                        aLeaderFieldInterpolator,
                        aFollowerFieldInterpolator );
            }
            else
            {
                return this->remap_nonconformal_rays_linear_deformed_geometry(
                        aLeaderFieldInterpolator,
                        aFollowerFieldInterpolator );
            }
        }
        else
        {
            return this->remap_nonconformal_rays_undeformed_geometry();
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IWG::remap_nonconformal_rays_linear_deformed_geometry(
            Field_Interpolator* aLeaderFieldInterpolator,
            Field_Interpolator* aFollowerFieldInterpolator )
    {
        // get the space dimension
        uint const tSpaceDim     = aLeaderFieldInterpolator->get_space_dim();
        uint const tNumTimeBases = aLeaderFieldInterpolator->get_number_of_time_bases();

        // create GapData pointer
        if ( mGapData == nullptr )
        {
            mGapData = std::make_unique< GapData >();
        }

        // check if the gap data need to be evaluated
        if ( !mGapData->mEval )
        {
            Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = mGapData->mEta.matrix_data();
            return tFollowerParamPoint;
        }

        // check that we are dealing with line elements
        MORIS_ERROR( tSpaceDim == 2,
                "IWG::remap_nonconformal_rays - only implemented for line elements." );
        MORIS_ERROR( tNumTimeBases == 1,
                "IWG::remap_nonconformal_rays - only implemented for constant time slabs." );

        // get IG geometry interpolator for leader and follower
        Geometry_Interpolator* tLeaderIGGI   = this->mLeaderFIManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = this->mFollowerFIManager->get_IG_geometry_interpolator();

        // get coordinates (global and parametric) of IG nodes
        const Matrix< DDRMat >& tLeaderIGNodes   = tLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat >& tFollowerIGNodes = tFollowerIGGI->get_space_coeff();

        const Matrix< DDRMat >& tLeaderIGNodesParam   = tLeaderIGGI->get_space_param_coeff();
        const Matrix< DDRMat >& tFollowerIGNodesParam = tFollowerIGGI->get_space_param_coeff();

        // get IP geometry interpolator for leader and follower
        // Geometry_Interpolator* tLeaderIPGI   = aLeaderFieldInterpolator->get_IG_geometry_interpolator();
        // Geometry_Interpolator* tFollowerIPGI = aFollowerFieldInterpolator->get_IG_geometry_interpolator();

        // get displacements for leader and follower
        const Matrix< DDRMat >& tLeaderUHat = aLeaderFieldInterpolator->get_coeff();
        // const Matrix< DDRMat >& tFollowerUHat = aFollowerFieldInterpolator->get_coeff();

        // determine number of IP nodes and dofs (assume same for leader and follower)
        const uint tNumNodes = tLeaderUHat.n_rows();
        const uint tNumDofs  = tLeaderUHat.numel();

        const uint tNumIGNodes = tLeaderIGNodes.n_rows();

        MORIS_ASSERT( tNumDofs == tNumNodes * tSpaceDim,
                "IWG::remap_nonconformal_rays - Number of dofs does not match number of nodes." );

        // get geometry for leader quadrature point
        const Matrix< DDRMat >& tLeaderXgp  = tLeaderIGGI->valx();    // assuming here that this is linear
        const Matrix< DDRMat >& tLeaderNXgp = tLeaderIGGI->NXi();
        // const Matrix< DDRMat >& tLeaderdNXgpdxi = tLeaderIGGI->dNdXi();

        // get displacements at IG Nodes of leader
        Matrix< DDRMat > tLeaderUIGNodes( tSpaceDim, tNumIGNodes );
        Matrix< DDRMat > tLeaderNUIGNodes( tSpaceDim, tNumDofs );

        const Matrix< DDRMat > tLeaderParamSpaceTime        = aLeaderFieldInterpolator->get_space_time();
        Matrix< DDRMat >       tLeaderIGNodesParamSpaceTime = tLeaderParamSpaceTime;

        for ( uint in = 0; in < tNumIGNodes; in++ )
        {
            tLeaderIGNodesParamSpaceTime( { 0, tSpaceDim - 1 }, { 0, 0 } ) =
                    trans( tLeaderIGNodesParam( { in, in }, { 0, tSpaceDim - 1 } ) );
            aLeaderFieldInterpolator->set_space_time( tLeaderIGNodesParamSpaceTime );

            tLeaderUIGNodes.get_column( in ) = aLeaderFieldInterpolator->val().matrix_data();

            const Matrix< DDRMat > tLeaderNUIGNodesTilde = aLeaderFieldInterpolator->N().matrix_data();

            tLeaderNUIGNodes( { in, in }, { 0, tNumDofs - 1 } ) = tLeaderNUIGNodesTilde( { 0, 0 }, { 0, tNumDofs - 1 } );
        }
        aLeaderFieldInterpolator->set_space_time( tLeaderParamSpaceTime );

        // compute displacement at IG nodes of leader
        const Matrix< DDRMat > tLeaderUgp = tLeaderUIGNodes * trans( tLeaderNXgp );

        // compute difference of shape functions at IG nodes
        const Matrix< DDRMat > tLeaderDifNUIGNodes =
                tLeaderNUIGNodes( { 1, 1 }, { 0, tNumDofs - 1 } )    //
                - tLeaderNUIGNodes( { 0, 0 }, { 0, tNumDofs - 1 } );

        // compute position of IG nodes in the deformed configuration
        const Matrix< DDRMat > tLeaderDefIGNodes = trans( tLeaderIGNodes ) + tLeaderUIGNodes;

        //        const Matrix< DDRMat >& tLeaderUgp     = aLeaderFieldInterpolator->val();
        //        const Matrix< DDRMat >& tLeaderNUgp    = tLeaderIPGI->NXi();
        //        const Matrix< DDRMat >& tLeaderdNUgpdr = tLeaderIPGI->dNdXi();
        //
        //        // compute derivatives along IG element side
        //        const Matrix< DDRMat > tLeaderdXgpdxi = trans( tLeaderdNXgpdxi * tLeaderIGNodes );
        //        const Matrix< DDRMat > tLeaderRgp     = tLeaderNXgp * tLeaderIGNodesParam;
        //        const Matrix< DDRMat > tLeaderdRgpdxi = tLeaderdNXgpdxi * tLeaderIGNodesParam;
        //        const Matrix< DDRMat > tLeaderdUgpdr  = tLeaderdNUgpdr * tLeaderUHat;
        //
        //        const Matrix< DDRMat > tLeaderdUgpdxi  = trans( tLeaderdRgpdxi * tLeaderdUgpdr );
        //        const Matrix< DDRMat > tLeaderdNUgpdxi = trans( tLeaderdRgpdxi * tLeaderdNUgpdr );


        // compute outwards pointing normal vector at the leader quadrature point
        const Matrix< DDRMat > RotMat             = { { 0, 1 }, { -1, 0 } };    // FIXME: this is only valid for 2D line elements
        const Matrix< DDRMat > tLeaderNormalTilde = RotMat * ( tLeaderDefIGNodes.get_column( 1 ) - tLeaderDefIGNodes.get_column( 0 ) );

        Matrix< DDRMat > tLeaderRefNormal = RotMat * trans( tLeaderIGNodes.get_row( 1 ) - tLeaderIGNodes.get_row( 0 ) );
        tLeaderRefNormal                  = 1.0 / norm( tLeaderRefNormal ) * tLeaderRefNormal;

        const real tNtildeSquared = dot( tLeaderNormalTilde, tLeaderNormalTilde );
        const real tNtildeOm12    = 1.0 / std::sqrt( tNtildeSquared );
        const real tNtildeOm32    = std::pow( tNtildeSquared, -1.5 );
        const real tNtilde3Om52   = 3.0 * std::pow( tNtildeSquared, -2.5 );

        const Matrix< DDRMat > tLeaderNormal = tNtildeOm12 * tLeaderNormalTilde;

        // compute first derivative of the normal vector with respect to the leader dofs
        const Matrix< DDRMat > tIdentity = eye( tSpaceDim, tSpaceDim );

        const Matrix< DDRMat > tPreMultiply =
                ( tNtildeOm12 * tIdentity - tNtildeOm32 * tLeaderNormalTilde * trans( tLeaderNormalTilde ) ) * RotMat;

        //        const Matrix< DDRMat > tLeaderdNormaldU =                                                           //
        //                tPreMultiply * ( tLeaderNUIGNodes( { tSpaceDim, tSpaceDim + 1 }, { 0, tNumDofs - 1 } ) -    //
        //                                 tLeaderNUIGNodes( { 0, tSpaceDim - 1 }, { 0, tNumDofs - 1 } ) );

        // Fixme: see for better solution above but requires reordering of the dof indices
        Matrix< DDRMat > tLeaderdNormaldU( tSpaceDim, tNumDofs );
        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    tLeaderDifNUIGNodes( in ) * tPreMultiply * tIdentity;
        }

        // compute second derivative of the normal vector with respect to the leader dofs
        const Matrix< DDRMat > tLeaderNormalAux = trans( RotMat ) * tLeaderNormalTilde;

        Matrix< DDRMat > tLeaderdNormal2dU2( tSpaceDim, tNumDofs * tNumDofs );

        uint tCounter = 0;
        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                for ( uint id = 0; id < tSpaceDim; id++ )
                {
                    for ( uint jn = 0; jn < tNumNodes; jn++ )
                    {
                        for ( uint jd = 0; jd < tSpaceDim; jd++ )
                        {
                            tLeaderdNormal2dU2( idim, tCounter ) =
                                    ( tNtilde3Om52 * tLeaderNormalTilde( idim ) * tLeaderNormalAux( jd ) * tLeaderNormalAux( id )    //
                                            - tNtildeOm32 * ( RotMat( idim, id ) * tLeaderNormalAux( jd )                            //
                                                              + RotMat( idim, jd ) * tLeaderNormalAux( id )                          //
                                                              + tLeaderNormalTilde( idim ) * tIdentity( id, jd ) ) )                 //
                                    * tLeaderDifNUIGNodes( in ) * tLeaderDifNUIGNodes( jn );
                            tCounter++;
                        }
                    }
                }
            }
            tCounter = 0;
        }

        // get displacements at IG Nodes of follower
        Matrix< DDRMat > tFollowerVIGNodes( tSpaceDim, tNumIGNodes );
        Matrix< DDRMat > tFollowerNVIGNodes( tSpaceDim, tNumDofs );

        // since follower IP is not initialized use the one of the leader
        const Matrix< DDRMat >& tFollowerParamSpaceTime        = tLeaderParamSpaceTime;
        Matrix< DDRMat >        tFollowerIGNodesParamSpaceTime = tFollowerParamSpaceTime;

        for ( uint in = 0; in < tNumIGNodes; in++ )
        {
            tFollowerIGNodesParamSpaceTime( { 0, tSpaceDim - 1 }, { 0, 0 } ) =
                    trans( tFollowerIGNodesParam( { in, in }, { 0, tSpaceDim - 1 } ) );
            aFollowerFieldInterpolator->set_space_time( tFollowerIGNodesParamSpaceTime );

            tFollowerVIGNodes.get_column( in ) =
                    aFollowerFieldInterpolator->val().matrix_data();

            const Matrix< DDRMat > tFollowerNVIGNodesTilde = aFollowerFieldInterpolator->N().matrix_data();

            tFollowerNVIGNodes( { in, in }, { 0, tNumDofs - 1 } ) = tFollowerNVIGNodesTilde( { 0, 0 }, { 0, tNumDofs - 1 } );
        }
        aFollowerFieldInterpolator->set_space_time( tFollowerParamSpaceTime );

        // compute difference of shape functions at IG nodes
        const Matrix< DDRMat > tFollowerDifNUIGNodes =
                tFollowerNVIGNodes( { tSpaceDim - 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } )    //
                - tFollowerNVIGNodes( { 0, 0 }, { 0, tNumDofs - 1 } );

        // compute ray intersection with follower element
        Matrix< DDRMat > tSolRay( tSpaceDim, 1, 0.0 );        // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tSolRayPrev( tSpaceDim, 1, 1.0 );    // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tJacRayInv;

        Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );

        real       tRayResRefNorm     = 0.0;
        const real tRayResRefNormTol  = 1e-10;
        const real tStagnationNormTol = 1e-12;
        const uint tMaxNewtonSteps    = 100;
        bool       tNewtonConverged   = false;

        for ( uint inew = 0; inew < tMaxNewtonSteps; inew++ )
        {
            // clip eta to be between -1 and 1
            tSolRay( 1 ) = std::min( 1.0, std::max( -1.0, tSolRay( 1 ) ) );    // FIXME: only works in 2D

            // set local parameter point for follower
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = tSolRay( { 1, tSpaceDim - 1 }, { 0, 0 } );

            this->mFollowerFIManager->set_space_time_from_local_IG_point( tFollowerParamPoint );

            // get geometry and displacements for follower quadrature point
            const Matrix< DDRMat >& tFollowerYgp       = tFollowerIGGI->valx();
            const Matrix< DDRMat >& tFollowerNYgp      = tFollowerIGGI->NXi();
            const Matrix< DDRMat >& tFollowerdNYgpdeta = tFollowerIGGI->dNdXi();

            const Matrix< DDRMat > tFollowerVgp = tFollowerVIGNodes * trans( tFollowerNYgp );

            Matrix< DDRMat > tFollowerdYgpdeta = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
            Matrix< DDRMat > tFollowerdVgpdeta = tFollowerVIGNodes * trans( tFollowerdNYgpdeta );

            // compute the residual
            Matrix< DDRMat > tResRay =
                    trans( tLeaderXgp - tFollowerYgp ) + tLeaderUgp - tFollowerVgp + tSolRay( 0 ) * tLeaderNormal;

            // store reference norm of the residual
            real tResRayNorm = norm( tResRay );
            if ( inew == 0 )
            {
                tRayResRefNorm = tResRayNorm;
            }

            // compute jacobian and its inverse
            Matrix< DDRMat > tJacRay( tSpaceDim, tSpaceDim );
            tJacRay( { 0, tSpaceDim - 1 }, { 0, 0 } )             = tLeaderNormal.matrix_data();
            tJacRay( { 0, tSpaceDim - 1 }, { 1, tSpaceDim - 1 } ) = -tFollowerdVgpdeta - tFollowerdYgpdeta;

            // check determination of jacobian
            if ( std::abs( det( tJacRay ) ) < MORIS_REAL_EPS )
            {
                break;
            }

            // compute inverse of jacobian
            tJacRayInv = inv( tJacRay );

            // check for convergence
            if ( tResRayNorm < tRayResRefNormTol * tRayResRefNorm )
            {
                // MORIS_LOG_INFO( "Remap - Newton converged in %d steps.\n", inew + 1 );
                tNewtonConverged = true;
                break;
            }

            // check for stagnation - stagnation is not considered convergence
            if ( norm( tSolRay - tSolRayPrev ) < tStagnationNormTol )
            {
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
                //                fprintf( stdout, "inew = %i, gap = %e, eta = %e, Resnorm = %e\n", inew, tSolRay( 0 ), tSolRay( 1 ), norm( tResRay ) );
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
                break;
            }
            tSolRayPrev = tSolRay;

            // update solution
            tSolRay -= tJacRayInv * tResRay;

            // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
            //            if ( std::abs( tSolRay( 1 ) ) > 1.0 )
            //            {
            //                fprintf( stdout, "inew = %i, gap = %e, eta = %e, Resnorm = %e\n", inew, tSolRay( 0 ), tSolRay( 1 ), norm( tResRay ) );
            //            }
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
        }

        bool tPlotRays = true;

        // check if Newton converged
        if ( !tNewtonConverged )
        {
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxx
            //            const Matrix< DDRMat > tLeaderDeformedIGNodes = trans( tLeaderIGNodes ) + tLeaderUIGNodes;
            //            print( tLeaderDeformedIGNodes, "leaderPoints" );
            //            const Matrix< DDRMat > tFollowerDeformedIGNodes = trans( tFollowerIGNodes ) + tFollowerVIGNodes;
            //            print( tFollowerDeformedIGNodes, "followerPoints" );
            //
            //            print( trans( tLeaderXgp ) + tLeaderUgp, "raycastingPoint" );
            //            print( tLeaderNormal, "normalVector" );
            //
            //            print( tLeaderIGGI->valx(), "udefqpPoint" );
            //            print( aLeaderFieldInterpolator->val(), "dispPoint" );
            //
            //            print( tJacRayInv, "tJacRayInv" );
            //            print( inv( tJacRayInv ), "tJacRay" );
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxx

            // set eta to outside bounds
            mGapData->mEta.set_size( tSpaceDim - 1, 1, -2.0 );

            // set gap data evaluation flag to true to avoid re-evaluations
            mGapData->mEval = false;

            // return a point with eta set to -2.0
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = -2.0;

            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if ( tPlotRays )
            {
                sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
                Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

                fprintf( stdout, "Niter = %d NotConverged %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTgtPoint( 0 ),
                        tTgtPoint( 1 ) );
            }
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            return tFollowerParamPoint;
        }

        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if ( tPlotRays )
        {
            sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
            Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
            Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

            fprintf( stdout, "Niter = %d Converged %e  %e  %e  %e\n",    //
                    tNiter,
                    tRayCastPoint( 0 ),
                    tRayCastPoint( 1 ),
                    tTgtPoint( 0 ),
                    tTgtPoint( 1 ) );
        }
        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        MORIS_ASSERT( tNewtonConverged, "IWG::remap_nonconformal_rays - Newton did not converge." );

        // check if valid intersection was found
        MORIS_ASSERT( tSolRay( 1 ) >= -1.0 && tSolRay( 1 ) <= 1.0,
                "IWG::remap_nonconformal_rays - No valid intersection found (eta(0) out side bounds)." );

        if ( tSpaceDim == 3 )
        {
            MORIS_ASSERT( tSolRay( 2 ) >= -1.0 && tSolRay( 2 ) <= 1.0,
                    "IWG::remap_nonconformal_rays - No valid intersection found (eta(1) out side bounds)." );
        }

        // store solution of ray intersection
        const real             tGap = tSolRay( 0 );
        const Matrix< DDRMat > tEta = tSolRay( { 1, tSpaceDim - 1 } );

        const Matrix< DDRMat > tGapVector = tGap * tLeaderNormal;

        // get global coordinates of  follower IP nodes
        // const Matrix< DDRMat >& tFollowerIPNodes = tFollowerIPGI->get_space_coeff();

        // get geometry and displacements for follower quadrature point
        // Note: the follower parameter point was already set within Newton loop
        // xxx const Matrix< DDRMat >& tFollowerYgp         = tFollowerIGGI->valx();
        const Matrix< DDRMat >& tFollowerNYgp        = tFollowerIGGI->NXi();
        const Matrix< DDRMat >& tFollowerdNYgpdeta   = tFollowerIGGI->dNdXi();
        const Matrix< DDRMat >& tFollowerdNYgp2deta2 = tFollowerIGGI->d2NdXi2();

        // xxxx const Matrix< DDRMat >& tFollowerVgp       = aFollowerFieldInterpolator->val();
        //        const Matrix< DDRMat >& tFollowerNVgp      = tFollowerIPGI->NXi();
        //        const Matrix< DDRMat >& tFollowerdNVgpdq   = tFollowerIPGI->dNdXi();
        //        const Matrix< DDRMat >& tFollowerdNVgp2dq2 = tFollowerIPGI->d2NdXi2();
        //
        //        const Matrix< DDRMat > tFollowerdYgpdq   = tFollowerdNVgpdq * tFollowerIPNodes;
        //        const Matrix< DDRMat > tFollowerdYgp2dq2 = tFollowerdNVgp2dq2 * tFollowerIPNodes;
        //
        //        const Matrix< DDRMat > tFollowerdYgpdeta = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
        //        // xxxx const Matrix< DDRMat > tFollowerQgp        = tFollowerNYgp * tFollowerIGNodesParam;
        //        const Matrix< DDRMat > tFollowerdQgpdeta   = tFollowerdNYgpdeta * tFollowerIGNodesParam;
        //        const Matrix< DDRMat > tFollowerdQgp2deta2 = tFollowerdNYgp2deta2 * tFollowerIGNodesParam;
        //        const Matrix< DDRMat > tFollowerdVgpdq     = tFollowerdNVgpdq * tFollowerUHat;
        //        const Matrix< DDRMat > tFollowerdVgp2dq2   = tFollowerdNVgp2dq2 * tFollowerUHat;

        // xxx const Matrix< DDRMat > tFollowerdVgpdeta  = trans( tFollowerdQgpdeta * tFollowerdVgpdq );
        //        const Matrix< DDRMat > tFollowerdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerdNVgpdq );
        //        const Matrix< DDRMat > tFollowerVgp = tFollowerVIGNodes * trans( tFollowerNYgp );

        Matrix< DDRMat > tFollowerdYgpdeta  = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
        Matrix< DDRMat > tFollowerdNVgpdeta = tFollowerVIGNodes * trans( tFollowerdNYgpdeta );

        const Matrix< DDRMat > tFollowerVgp = tFollowerVIGNodes * trans( tFollowerNYgp );

        // compute second derivative of follower velocity with respect to local coordinate on follower side
        const Matrix< DDRMat > tFollowerdVgp2deta2 = trans( tFollowerdNYgp2deta2 * tFollowerVIGNodes );

        // compute RHS for sensitivity equations to compute the derivatives of gap and eta wrt follower displacement dofs
        Matrix< DDRMat > tFollowerdResRaydv( tSpaceDim, tNumDofs );
        Matrix< DDRMat > tFollowerdResRay2detadv( tSpaceDim, tNumDofs );

        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tFollowerdResRaydv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    dot( tFollowerNYgp, tFollowerNVIGNodes.get_column( in ) ) * tIdentity;
            tFollowerdResRay2detadv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    dot( tFollowerdNYgpdeta, tFollowerNVIGNodes.get_column( in ) ) * tIdentity;
        }

        // compute the first derivatives of gap and eta wrt follower displacement dofs
        const Matrix< DDRMat > tdSolRaydv = tJacRayInv * tFollowerdResRaydv;

        const Matrix< DDRMat > tdGapdv = tdSolRaydv.get_row( 0 );
        const Matrix< DDRMat > tdEtadv = tdSolRaydv( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );

        // compute the second derivatives of gap and eta wrt follower displacement dofs
        Matrix< DDRMat > tdGap2dv2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2dv2( tNumDofs, tNumDofs );

        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                Matrix< DDRMat > tdSolRay2dv2 =
                        tJacRayInv * ( tdEtadv( 0, idof ) * tFollowerdResRay2detadv.get_column( jdof ) +    //
                                       tdEtadv( 0, jdof ) * tFollowerdResRay2detadv.get_column( idof ) +    //
                                       tFollowerdVgp2deta2 * tdEtadv( 0, idof ) * tdEtadv( 0, jdof ) );

                tdGap2dv2( idof, jdof ) = tdSolRay2dv2( 0 );
                tdEta2dv2( idof, jdof ) = tdSolRay2dv2( 1 );    // FIXME: this is only valid for 2D line elements
            }
        }

        // compute RHS for sensitivity equations to compute the derivatives of gap and eta wrt leader displacement dofs
        Matrix< DDRMat > tLeaderdResRaydu( tSpaceDim, tNumDofs );

        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tLeaderdResRaydu( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =    //
                    -dot( tLeaderNXgp, tLeaderNUIGNodes.get_column( in ) ) * tIdentity            //
                    - tGap * tLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } );
        }

        // compute the derivatives of gap and eta wrt leader displacement dofs
        const Matrix< DDRMat > tdSolRaydu = tJacRayInv * tLeaderdResRaydu;

        const Matrix< DDRMat > tdGapdu = tdSolRaydu.get_row( 0 );
        const Matrix< DDRMat > tdEtadu = tdSolRaydu( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );

        // compute the second derivatives of gap and eta wrt leader displacement dofs
        Matrix< DDRMat > tdGap2du2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2du2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdGap2duv( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2duv( tNumDofs, tNumDofs );

        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                Matrix< DDRMat > tdSolRay2du2 =
                        tJacRayInv * ( -tGap * tLeaderdNormal2dU2.get_column( idof * tNumDofs + jdof )    //
                                       - tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )            //
                                       - tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )            //
                                       + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadu( 0, jdof ) );

                tdGap2du2( idof, jdof ) = tdSolRay2du2( 0 );
                tdEta2du2( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements

                tdSolRay2du2 = tJacRayInv * ( -tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof )               //
                                              + tFollowerdResRay2detadv.get_column( jdof ) * tdEtadu( 0, idof )    //
                                              + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadv( 0, jdof ) );

                tdGap2duv( idof, jdof ) = tdSolRay2du2( 0 );
                tdEta2duv( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements
            }
        }

        // compute the derivatives of the gap vector
        Matrix< DDRMat > tdGapvecdu = tGap * tLeaderdNormaldU + tLeaderNormal * tdGapdu;
        Matrix< DDRMat > tdGapvecdv = tLeaderNormal * tdGapdv;

        Matrix< DDRMat > tdGapvec2du2( tSpaceDim, tNumDofs * tNumDofs );
        Matrix< DDRMat > tdGapvec2duv( tSpaceDim, tNumDofs * tNumDofs );
        Matrix< DDRMat > tdGapvec2dv2 = tLeaderNormal * trans( vectorize( tdGap2dv2 ) );

        tCounter = 0;
        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                tdGapvec2du2.get_column( tCounter ) = tGap * tLeaderdNormal2dU2.get_column( tCounter )         //
                                                    + tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )    //
                                                    + tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )    //
                                                    + tLeaderNormal * tdGap2du2( idof, jdof );

                tdGapvec2duv.get_column( tCounter ) = tLeaderNormal * tdGap2duv( idof, jdof )    //
                                                    + tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof );
                tCounter++;
            }
        }

        // save results in GapData resorted such that nodes first
        uint tIcounter = 0;
        uint tJcounter = 0;

        mGapData->mGap             = tGap;
        mGapData->mEta             = tEta;
        mGapData->mLeaderNormal    = tLeaderNormal;
        mGapData->mLeaderRefNormal = tLeaderRefNormal;
        mGapData->mGapVec          = tGapVector;

        mGapData->mdGapdu.set_size( 1, tNumDofs );
        mGapData->mdGapdv.set_size( 1, tNumDofs );
        mGapData->mdEtadu.set_size( 1, tNumDofs );
        mGapData->mdEtadv.set_size( 1, tNumDofs );

        mGapData->mLeaderdNormaldu.set_size( tSpaceDim, tNumDofs );
        mGapData->mdGapvecdu.set_size( tSpaceDim, tNumDofs );
        mGapData->mdGapvecdv.set_size( tSpaceDim, tNumDofs );

        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                mGapData->mdGapdu( tIcounter ) = tdGapdu( idim + in * tSpaceDim );
                mGapData->mdGapdv( tIcounter ) = tdGapdv( idim + in * tSpaceDim );
                mGapData->mdEtadu( tIcounter ) = tdEtadu( idim + in * tSpaceDim );
                mGapData->mdEtadv( tIcounter ) = tdEtadv( idim + in * tSpaceDim );

                mGapData->mLeaderdNormaldu.get_column( tIcounter ) = tLeaderdNormaldU.get_column( idim + in * tSpaceDim );
                mGapData->mdGapvecdu.get_column( tIcounter )       = tdGapvecdu.get_column( idim + in * tSpaceDim );
                mGapData->mdGapvecdv.get_column( tIcounter )       = tdGapvecdv.get_column( idim + in * tSpaceDim );

                tIcounter++;
            }
        }

        mGapData->mdGap2du2.set_size( tNumDofs, tNumDofs );
        mGapData->mdGap2dv2.set_size( tNumDofs, tNumDofs );
        mGapData->mdGap2duv.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2du2.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2dv2.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2duv.set_size( tNumDofs, tNumDofs );

        mGapData->mLeaderdNormal2du2.set_size( tSpaceDim, tNumDofs * tNumDofs );    // tLeaderdNormal2dU2
        mGapData->mdGapvec2du2.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2du2
        mGapData->mdGapvec2dv2.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2dv2
        mGapData->mdGapvec2duv.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2duv

        tIcounter = 0;
        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                tJcounter = 0;
                for ( uint jdim = 0; jdim < tSpaceDim; jdim++ )
                {
                    for ( uint jn = 0; jn < tNumNodes; jn++ )
                    {
                        mGapData->mdGap2du2( tIcounter, tJcounter ) = tdGap2du2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdGap2dv2( tIcounter, tJcounter ) = tdGap2dv2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdGap2duv( tIcounter, tJcounter ) = tdGap2duv( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2du2( tIcounter, tJcounter ) = tdEta2du2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2dv2( tIcounter, tJcounter ) = tdEta2dv2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2duv( tIcounter, tJcounter ) = tdEta2duv( idim + in * tSpaceDim, jdim + jn * tSpaceDim );

                        mGapData->mLeaderdNormal2du2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tLeaderdNormal2dU2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2du2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2du2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2dv2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2dv2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2duv.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2duv.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        tJcounter++;
                    }
                }
                tIcounter++;
            }
        }

        // gap data evaluation flag to false, i.e. does not need to be recomputed
        mGapData->mEval = false;

        // return parametric point on follower side
        return tFollowerParamPoint;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IWG::remap_nonconformal_rays_consistent_deformed_geometry(
            Field_Interpolator* aLeaderFieldInterpolator,
            Field_Interpolator* aFollowerFieldInterpolator )
    {
        // get the space dimension
        uint const tSpaceDim     = aLeaderFieldInterpolator->get_space_dim();
        uint const tNumTimeBases = aLeaderFieldInterpolator->get_number_of_time_bases();

        // create GapData pointer
        if ( mGapData == nullptr )
        {
            mGapData = std::make_unique< GapData >();
        }

        // check if the gap data need to be evaluated
        if ( !mGapData->mEval )
        {
            Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = mGapData->mEta.matrix_data();
            return tFollowerParamPoint;
        }

        // check that we are dealing with line elements
        MORIS_ERROR( tSpaceDim == 2,
                "IWG::remap_nonconformal_rays - only implemented for line elements." );
        MORIS_ERROR( tNumTimeBases == 1,
                "IWG::remap_nonconformal_rays - only implemented for constant time slabs." );

        // get IG geometry interpolator for leader and follower
        Geometry_Interpolator* tLeaderIGGI   = this->mLeaderFIManager->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI = this->mFollowerFIManager->get_IG_geometry_interpolator();

        // get coordinates (global and parametric) of IG nodes
        const Matrix< DDRMat >& tLeaderIGNodes   = tLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat >& tFollowerIGNodes = tFollowerIGGI->get_space_coeff();

        const Matrix< DDRMat >& tLeaderIGNodesParam   = tLeaderIGGI->get_space_param_coeff();
        const Matrix< DDRMat >& tFollowerIGNodesParam = tFollowerIGGI->get_space_param_coeff();

        // get IP geometry interpolator for leader and follower
        Geometry_Interpolator* tLeaderIPGI   = aLeaderFieldInterpolator->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIPGI = aFollowerFieldInterpolator->get_IG_geometry_interpolator();

        // get displacements for leader and follower
        const Matrix< DDRMat >& tLeaderUHat   = aLeaderFieldInterpolator->get_coeff();
        const Matrix< DDRMat >& tFollowerUHat = aFollowerFieldInterpolator->get_coeff();

        // determine number of IP nodes and dofs (assume same for leader and follower)
        const uint tNumNodes = tLeaderUHat.n_rows();
        const uint tNumDofs  = tLeaderUHat.numel();

        MORIS_ASSERT( tNumDofs == tNumNodes * tSpaceDim,
                "IWG::remap_nonconformal_rays - Number of dofs does not match number of nodes." );

        // get geometry and displacements for leader quadrature point
        const Matrix< DDRMat >& tLeaderXgp      = tLeaderIGGI->valx();
        const Matrix< DDRMat >& tLeaderNXgp     = tLeaderIGGI->NXi();
        const Matrix< DDRMat >& tLeaderdNXgpdxi = tLeaderIGGI->dNdXi();

        const Matrix< DDRMat >& tLeaderUgp     = aLeaderFieldInterpolator->val();
        const Matrix< DDRMat >& tLeaderNUgp    = tLeaderIPGI->NXi();
        const Matrix< DDRMat >& tLeaderdNUgpdr = tLeaderIPGI->dNdXi();

        // compute derivatives along IG element side
        const Matrix< DDRMat > tLeaderdXgpdxi = trans( tLeaderdNXgpdxi * tLeaderIGNodes );
        const Matrix< DDRMat > tLeaderRgp     = tLeaderNXgp * tLeaderIGNodesParam;
        const Matrix< DDRMat > tLeaderdRgpdxi = tLeaderdNXgpdxi * tLeaderIGNodesParam;
        const Matrix< DDRMat > tLeaderdUgpdr  = tLeaderdNUgpdr * tLeaderUHat;

        const Matrix< DDRMat > tLeaderdUgpdxi  = trans( tLeaderdRgpdxi * tLeaderdUgpdr );
        const Matrix< DDRMat > tLeaderdNUgpdxi = trans( tLeaderdRgpdxi * tLeaderdNUgpdr );

        // compute outwards pointing normal vector at the leader quadrature point
        const Matrix< DDRMat > RotMat             = { { 0, 1 }, { -1, 0 } };    // FIXME: this is only valid for 2D line elements
        const Matrix< DDRMat > tLeaderNormalTilde = RotMat * ( tLeaderdXgpdxi + tLeaderdUgpdxi );

        Matrix< DDRMat > tLeaderRefNormal = RotMat * ( tLeaderdXgpdxi );
        tLeaderRefNormal                  = 1.0 / norm( tLeaderRefNormal ) * tLeaderRefNormal;

        const real tNtildeSquared = dot( tLeaderNormalTilde, tLeaderNormalTilde );
        const real tNtildeOm12    = 1.0 / std::sqrt( tNtildeSquared );
        const real tNtildeOm32    = std::pow( tNtildeSquared, -1.5 );
        const real tNtilde3Om52   = 3.0 * std::pow( tNtildeSquared, -2.5 );

        const Matrix< DDRMat > tLeaderNormal = tNtildeOm12 * tLeaderNormalTilde;

        // compute first derivative of the normal vector with respect to the leader dofs
        const Matrix< DDRMat > tIdentity = eye( tSpaceDim, tSpaceDim );

        const Matrix< DDRMat > tPreMultiply =
                ( tNtildeOm12 * tIdentity - tNtildeOm32 * tLeaderNormalTilde * trans( tLeaderNormalTilde ) ) * RotMat;

        Matrix< DDRMat > tLeaderdNormaldU( tSpaceDim, tNumDofs );
        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    tLeaderdNUgpdxi( in ) * tPreMultiply * tIdentity;
        }

        // compute second derivative of the normal vector with respect to the leader dofs
        const Matrix< DDRMat > tLeaderNormalAux = trans( RotMat ) * tLeaderNormalTilde;

        Matrix< DDRMat > tLeaderdNormal2dU2( tSpaceDim, tNumDofs * tNumDofs );

        uint tCounter = 0;
        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                for ( uint id = 0; id < tSpaceDim; id++ )
                {
                    for ( uint jn = 0; jn < tNumNodes; jn++ )
                    {
                        for ( uint jd = 0; jd < tSpaceDim; jd++ )
                        {
                            tLeaderdNormal2dU2( idim, tCounter ) =
                                    ( tNtilde3Om52 * tLeaderNormalTilde( idim ) * tLeaderNormalAux( jd ) * tLeaderNormalAux( id )    //
                                            - tNtildeOm32 * ( RotMat( idim, id ) * tLeaderNormalAux( jd )                            //
                                                              + RotMat( idim, jd ) * tLeaderNormalAux( id )                          //
                                                              + tLeaderNormalTilde( idim ) * tIdentity( id, jd ) ) )                 //
                                    * tLeaderdNUgpdxi( in ) * tLeaderdNUgpdxi( jn );
                            tCounter++;
                        }
                    }
                }
            }
            tCounter = 0;
        }

        // compute ray intersection with follower element
        Matrix< DDRMat > tSolRay( tSpaceDim, 1, 0.0 );        // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tSolRayPrev( tSpaceDim, 1, 1.0 );    // FIXME: only implemented for constant time slabs
        Matrix< DDRMat > tJacRayInv;

        Matrix< DDRMat > tFollowerParamPoint( tSpaceDim, 1, 0.0 );

        real       tRayResRefNorm     = 0.0;
        const real tRayResRefNormTol  = 1e-10;
        const real tStagnationNormTol = 1e-12;
        const uint tMaxNewtonSteps    = 100;
        bool       tNewtonConverged   = false;

        for ( uint inew = 0; inew < tMaxNewtonSteps; inew++ )
        {
            // clip eta to be between -1 and 1
            tSolRay( 1 ) = std::min( 1.0, std::max( -1.0, tSolRay( 1 ) ) );    // FIXME: only works in 2D

            // set local parameter point for follower
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = tSolRay( { 1, tSpaceDim - 1 }, { 0, 0 } );

            this->mFollowerFIManager->set_space_time_from_local_IG_point( tFollowerParamPoint );

            // get geometry and displacements for follower quadrature point
            const Matrix< DDRMat >& tFollowerYgp = tFollowerIGGI->valx();
            // xxxx const Matrix< DDRMat >& tFollowerNYgp      = tFollowerIGGI->NXi();
            const Matrix< DDRMat >& tFollowerdNYgpdeta = tFollowerIGGI->dNdXi();
            const Matrix< DDRMat >& tFollowerVgp       = aFollowerFieldInterpolator->val();
            // xxxx const Matrix< DDRMat >& tFollowerNVgp      = tFollowerIPGI->NXi();
            const Matrix< DDRMat >& tFollowerdNVgpdq = tFollowerIPGI->dNdXi();

            Matrix< DDRMat > tFollowerdYgpdeta = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
            // xxxx Matrix< DDRMat > tFollowerQgp      = tFollowerNYgp * tFollowerIGNodesParam;
            Matrix< DDRMat > tFollowerdQgpdeta = tFollowerdNYgpdeta * tFollowerIGNodesParam;
            Matrix< DDRMat > tFollowerdVgpdq   = tFollowerdNVgpdq * tFollowerUHat;

            Matrix< DDRMat > tFollowerdVgpdeta = trans( tFollowerdQgpdeta * tFollowerdVgpdq );
            // xxxx Matrix< DDRMat > tFollowerdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerdNVgpdq );

            // compute the residual
            Matrix< DDRMat > tResRay =
                    trans( tLeaderXgp - tFollowerYgp ) + tLeaderUgp - tFollowerVgp + tSolRay( 0 ) * tLeaderNormal;

            // store reference norm of the residual
            real tResRayNorm = norm( tResRay );
            if ( inew == 0 )
            {
                tRayResRefNorm = tResRayNorm;
            }

            // compute jacobian and its inverse
            Matrix< DDRMat > tJacRay( tSpaceDim, tSpaceDim );
            tJacRay( { 0, tSpaceDim - 1 }, { 0, 0 } )             = tLeaderNormal.matrix_data();
            tJacRay( { 0, tSpaceDim - 1 }, { 1, tSpaceDim - 1 } ) = -tFollowerdVgpdeta - tFollowerdYgpdeta;

            // check determination of jacobian
            if ( std::abs( det( tJacRay ) ) < MORIS_REAL_EPS )
            {
                break;
            }

            // compute inverse of jacobian
            tJacRayInv = inv( tJacRay );

            // check for convergence
            if ( tResRayNorm < tRayResRefNormTol * tRayResRefNorm )
            {
                // MORIS_LOG_INFO( "Remap - Newton converged in %d steps.\n", inew + 1 );
                tNewtonConverged = true;
                break;
            }

            // check for stagnation - stagnation is not considered convergence
            if ( norm( tSolRay - tSolRayPrev ) < tStagnationNormTol )
            {
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                //                 fprintf( stdout, "inew = %i, gap = %e, eta = %e, Resnorm = %e\n", inew, tSolRay( 0 ), tSolRay( 1 ), norm( tResRay ) );
                // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                break;
            }
            tSolRayPrev = tSolRay;

            // update solution
            tSolRay -= tJacRayInv * tResRay;
        }

        bool tPlotRays = true;

        // check if Newton converged
        if ( !tNewtonConverged )
        {
            // set eta to outside bounds
            mGapData->mEta.set_size( tSpaceDim - 1, 1, -2.0 );

            // set gap data evaluation flag to true to avoid re-evaluations
            mGapData->mEval = false;

            // return a point with eta set to -2.0
            tFollowerParamPoint( { 0, tSpaceDim - 2 }, { 0, 0 } ) = -2.0;

            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
            if ( tPlotRays )
            {
                sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
                Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
                Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

                fprintf( stdout, "Niter = %d NotConverged %e  %e  %e  %e\n",    //
                        tNiter,
                        tRayCastPoint( 0 ),
                        tRayCastPoint( 1 ),
                        tTgtPoint( 0 ),
                        tTgtPoint( 1 ) );
            }
            // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

            return tFollowerParamPoint;
        }

        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        if ( tPlotRays )
        {
            sint             tNiter        = (sint)gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve", true );
            Matrix< DDRMat > tRayCastPoint = trans( tLeaderXgp ) + tLeaderUgp;
            Matrix< DDRMat > tTgtPoint     = tRayCastPoint + tSolRay( 0 ) * tLeaderNormal;

            fprintf( stdout, "Niter = %d Converged %e  %e  %e  %e\n",    //
                    tNiter,
                    tRayCastPoint( 0 ),
                    tRayCastPoint( 1 ),
                    tTgtPoint( 0 ),
                    tTgtPoint( 1 ) );
        }
        // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

        MORIS_ASSERT( tNewtonConverged, "IWG::remap_nonconformal_rays - Newton did not converge." );

        // check if valid intersection was found
        MORIS_ASSERT( tSolRay( 1 ) >= -1.0 && tSolRay( 1 ) <= 1.0,
                "IWG::remap_nonconformal_rays - No valid intersection found (eta(0) out side bounds)." );

        if ( tSpaceDim == 3 )
        {
            MORIS_ASSERT( tSolRay( 2 ) >= -1.0 && tSolRay( 2 ) <= 1.0,
                    "IWG::remap_nonconformal_rays - No valid intersection found (eta(1) out side bounds)." );
        }

        // store solution of ray intersection
        const real             tGap = tSolRay( 0 );
        const Matrix< DDRMat > tEta = tSolRay( { 1, tSpaceDim - 1 } );

        const Matrix< DDRMat > tGapVector = tGap * tLeaderNormal;

        // get global coordinates of  follower IP nodes
        const Matrix< DDRMat >& tFollowerIPNodes = tFollowerIPGI->get_space_coeff();

        // get geometry and displacements for follower quadrature point
        // Note: the follower parameter point was already set within Newton loop
        // xxx const Matrix< DDRMat >& tFollowerYgp         = tFollowerIGGI->valx();
        // xxx const Matrix< DDRMat >& tFollowerNYgp        = tFollowerIGGI->NXi();
        const Matrix< DDRMat >& tFollowerdNYgpdeta   = tFollowerIGGI->dNdXi();
        const Matrix< DDRMat >& tFollowerdNYgp2deta2 = tFollowerIGGI->d2NdXi2();

        // xxxx const Matrix< DDRMat >& tFollowerVgp       = aFollowerFieldInterpolator->val();
        const Matrix< DDRMat >& tFollowerNVgp      = tFollowerIPGI->NXi();
        const Matrix< DDRMat >& tFollowerdNVgpdq   = tFollowerIPGI->dNdXi();
        const Matrix< DDRMat >& tFollowerdNVgp2dq2 = tFollowerIPGI->d2NdXi2();

        const Matrix< DDRMat > tFollowerdYgpdq   = tFollowerdNVgpdq * tFollowerIPNodes;
        const Matrix< DDRMat > tFollowerdYgp2dq2 = tFollowerdNVgp2dq2 * tFollowerIPNodes;

        const Matrix< DDRMat > tFollowerdYgpdeta = trans( tFollowerdNYgpdeta * tFollowerIGNodes );
        // xxxx const Matrix< DDRMat > tFollowerQgp        = tFollowerNYgp * tFollowerIGNodesParam;
        const Matrix< DDRMat > tFollowerdQgpdeta   = tFollowerdNYgpdeta * tFollowerIGNodesParam;
        const Matrix< DDRMat > tFollowerdQgp2deta2 = tFollowerdNYgp2deta2 * tFollowerIGNodesParam;
        const Matrix< DDRMat > tFollowerdVgpdq     = tFollowerdNVgpdq * tFollowerUHat;
        const Matrix< DDRMat > tFollowerdVgp2dq2   = tFollowerdNVgp2dq2 * tFollowerUHat;

        // xxx const Matrix< DDRMat > tFollowerdVgpdeta  = trans( tFollowerdQgpdeta * tFollowerdVgpdq );
        const Matrix< DDRMat > tFollowerdNVgpdeta = trans( tFollowerdQgpdeta * tFollowerdNVgpdq );

        // Define matrix to convert 2nd derivative matrix to vector
        // FIXME: this is only valid for 2D line elements
        const Matrix< DDRMat > tExtractMat = { { 1, 0, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 0, 1, 0 } };

        // compute second derivative of follower velocity with respect to local coordinate on follower side
        const Matrix< DDRMat > tFollowerdVgp2deta2 =
                trans( tExtractMat * ( tFollowerdYgp2dq2 + tFollowerdVgp2dq2 ) )         //
                        * vectorize( trans( tFollowerdQgpdeta ) * tFollowerdQgpdeta )    //
                + trans( tFollowerdQgp2deta2 * ( tFollowerdYgpdq + tFollowerdVgpdq ) );

        // compute RHS for sensitivity equations to compute the derivatives of gap and eta wrt follower displacement dofs
        Matrix< DDRMat > tFollowerdResRaydv( tSpaceDim, tNumDofs );
        Matrix< DDRMat > tFollowerdResRay2detadv( tSpaceDim, tNumDofs );

        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tFollowerdResRaydv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    tFollowerNVgp( in ) * tIdentity;
            tFollowerdResRay2detadv( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =
                    tFollowerdNVgpdeta( in ) * tIdentity;
        }

        // compute the first derivatives of gap and eta wrt follower displacement dofs
        const Matrix< DDRMat > tdSolRaydv = tJacRayInv * tFollowerdResRaydv;

        const Matrix< DDRMat > tdGapdv = tdSolRaydv.get_row( 0 );
        const Matrix< DDRMat > tdEtadv = tdSolRaydv( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );

        // compute the second derivatives of gap and eta wrt follower displacement dofs
        Matrix< DDRMat > tdGap2dv2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2dv2( tNumDofs, tNumDofs );

        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                Matrix< DDRMat > tdSolRay2dv2 =
                        tJacRayInv * ( tdEtadv( 0, idof ) * tFollowerdResRay2detadv.get_column( jdof ) +    //
                                       tdEtadv( 0, jdof ) * tFollowerdResRay2detadv.get_column( idof ) +    //
                                       tFollowerdVgp2deta2 * tdEtadv( 0, idof ) * tdEtadv( 0, jdof ) );

                tdGap2dv2( idof, jdof ) = tdSolRay2dv2( 0 );
                tdEta2dv2( idof, jdof ) = tdSolRay2dv2( 1 );    // FIXME: this is only valid for 2D line elements
            }
        }

        // compute RHS for sensitivity equations to compute the derivatives of gap and eta wrt leader displacement dofs
        Matrix< DDRMat > tLeaderdResRaydu( tSpaceDim, tNumDofs );

        for ( uint in = 0; in < tNumNodes; in++ )
        {
            tLeaderdResRaydu( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } ) =    //
                    -tLeaderNUgp( in ) * tIdentity                                                //
                    - tGap * tLeaderdNormaldU( { 0, tSpaceDim - 1 }, { in * tSpaceDim, in * tSpaceDim + 1 } );
        }

        // compute the derivatives of gap and eta wrt leader displacement dofs
        const Matrix< DDRMat > tdSolRaydu = tJacRayInv * tLeaderdResRaydu;

        const Matrix< DDRMat > tdGapdu = tdSolRaydu.get_row( 0 );
        const Matrix< DDRMat > tdEtadu = tdSolRaydu( { 1, tSpaceDim - 1 }, { 0, tNumDofs - 1 } );

        // compute the second derivatives of gap and eta wrt leader displacement dofs
        Matrix< DDRMat > tdGap2du2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2du2( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdGap2duv( tNumDofs, tNumDofs );
        Matrix< DDRMat > tdEta2duv( tNumDofs, tNumDofs );

        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                Matrix< DDRMat > tdSolRay2du2 =
                        tJacRayInv * ( -tGap * tLeaderdNormal2dU2.get_column( idof * tNumDofs + jdof )    //
                                       - tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )            //
                                       - tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )            //
                                       + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadu( 0, jdof ) );

                tdGap2du2( idof, jdof ) = tdSolRay2du2( 0 );
                tdEta2du2( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements

                tdSolRay2du2 = tJacRayInv * ( -tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof )               //
                                              + tFollowerdResRay2detadv.get_column( jdof ) * tdEtadu( 0, idof )    //
                                              + tFollowerdVgp2deta2 * tdEtadu( 0, idof ) * tdEtadv( 0, jdof ) );

                tdGap2duv( idof, jdof ) = tdSolRay2du2( 0 );
                tdEta2duv( idof, jdof ) = tdSolRay2du2( 1 );    // FIXME: this is only valid for 2D line elements
            }
        }

        // compute the derivatives of the gap vector
        Matrix< DDRMat > tdGapvecdu = tGap * tLeaderdNormaldU + tLeaderNormal * tdGapdu;
        Matrix< DDRMat > tdGapvecdv = tLeaderNormal * tdGapdv;

        Matrix< DDRMat > tdGapvec2du2( tSpaceDim, tNumDofs * tNumDofs );
        Matrix< DDRMat > tdGapvec2duv( tSpaceDim, tNumDofs * tNumDofs );
        Matrix< DDRMat > tdGapvec2dv2 = tLeaderNormal * trans( vectorize( tdGap2dv2 ) );

        tCounter = 0;
        for ( uint idof = 0; idof < tNumDofs; idof++ )
        {
            for ( uint jdof = 0; jdof < tNumDofs; jdof++ )
            {
                tdGapvec2du2.get_column( tCounter ) = tGap * tLeaderdNormal2dU2.get_column( tCounter )         //
                                                    + tLeaderdNormaldU.get_column( idof ) * tdGapdu( jdof )    //
                                                    + tLeaderdNormaldU.get_column( jdof ) * tdGapdu( idof )    //
                                                    + tLeaderNormal * tdGap2du2( idof, jdof );

                tdGapvec2duv.get_column( tCounter ) = tLeaderNormal * tdGap2duv( idof, jdof )    //
                                                    + tLeaderdNormaldU.get_column( idof ) * tdGapdv( jdof );
                tCounter++;
            }
        }

        // save results in GapData resorted such that nodes first
        uint tIcounter = 0;
        uint tJcounter = 0;

        mGapData->mGap             = tGap;
        mGapData->mEta             = tEta;
        mGapData->mLeaderNormal    = tLeaderNormal;
        mGapData->mLeaderRefNormal = tLeaderRefNormal;
        mGapData->mGapVec          = tGapVector;

        mGapData->mdGapdu.set_size( 1, tNumDofs );
        mGapData->mdGapdv.set_size( 1, tNumDofs );
        mGapData->mdEtadu.set_size( 1, tNumDofs );
        mGapData->mdEtadv.set_size( 1, tNumDofs );

        mGapData->mLeaderdNormaldu.set_size( tSpaceDim, tNumDofs );
        mGapData->mdGapvecdu.set_size( tSpaceDim, tNumDofs );
        mGapData->mdGapvecdv.set_size( tSpaceDim, tNumDofs );

        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                mGapData->mdGapdu( tIcounter ) = tdGapdu( idim + in * tSpaceDim );
                mGapData->mdGapdv( tIcounter ) = tdGapdv( idim + in * tSpaceDim );
                mGapData->mdEtadu( tIcounter ) = tdEtadu( idim + in * tSpaceDim );
                mGapData->mdEtadv( tIcounter ) = tdEtadv( idim + in * tSpaceDim );

                mGapData->mLeaderdNormaldu.get_column( tIcounter ) = tLeaderdNormaldU.get_column( idim + in * tSpaceDim );
                mGapData->mdGapvecdu.get_column( tIcounter )       = tdGapvecdu.get_column( idim + in * tSpaceDim );
                mGapData->mdGapvecdv.get_column( tIcounter )       = tdGapvecdv.get_column( idim + in * tSpaceDim );

                tIcounter++;
            }
        }

        mGapData->mdGap2du2.set_size( tNumDofs, tNumDofs );
        mGapData->mdGap2dv2.set_size( tNumDofs, tNumDofs );
        mGapData->mdGap2duv.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2du2.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2dv2.set_size( tNumDofs, tNumDofs );
        mGapData->mdEta2duv.set_size( tNumDofs, tNumDofs );

        mGapData->mLeaderdNormal2du2.set_size( tSpaceDim, tNumDofs * tNumDofs );    // tLeaderdNormal2dU2
        mGapData->mdGapvec2du2.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2du2
        mGapData->mdGapvec2dv2.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2dv2
        mGapData->mdGapvec2duv.set_size( tSpaceDim, tNumDofs * tNumDofs );          // tdGapvec2duv

        tIcounter = 0;
        for ( uint idim = 0; idim < tSpaceDim; idim++ )
        {
            for ( uint in = 0; in < tNumNodes; in++ )
            {
                tJcounter = 0;
                for ( uint jdim = 0; jdim < tSpaceDim; jdim++ )
                {
                    for ( uint jn = 0; jn < tNumNodes; jn++ )
                    {
                        mGapData->mdGap2du2( tIcounter, tJcounter ) = tdGap2du2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdGap2dv2( tIcounter, tJcounter ) = tdGap2dv2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdGap2duv( tIcounter, tJcounter ) = tdGap2duv( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2du2( tIcounter, tJcounter ) = tdEta2du2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2dv2( tIcounter, tJcounter ) = tdEta2dv2( idim + in * tSpaceDim, jdim + jn * tSpaceDim );
                        mGapData->mdEta2duv( tIcounter, tJcounter ) = tdEta2duv( idim + in * tSpaceDim, jdim + jn * tSpaceDim );

                        mGapData->mLeaderdNormal2du2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tLeaderdNormal2dU2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2du2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2du2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2dv2.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2dv2.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        mGapData->mdGapvec2duv.get_column( tIcounter * tNumDofs + tJcounter ) =
                                tdGapvec2duv.get_column( ( idim + in * tSpaceDim ) * tNumDofs + jdim + jn * tSpaceDim );
                        tJcounter++;
                    }
                }
                tIcounter++;
            }
        }

        // gap data evaluation flag to false, i.e. does not need to be recomputed
        mGapData->mEval = false;

        // return parametric point on follower side
        return tFollowerParamPoint;
    }

    //------ ------------------------------------------------------------------------

    bool IWG::check_jacobian(
            real              aPerturbation,
            real              aEpsilon,
            real              aWStar,
            Matrix< DDRMat >& aJacobian,
            Matrix< DDRMat >& aJacobianFD,
            bool              aErrorPrint,
            bool              aUseAbsolutePerturbations )
    {
        // get residual dof type index in set, start and end indices for residual dof type
        uint tLeaderDofIndex    = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartRow = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResEndRow   = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        uint tFollowerDofIndex;
        uint tFollowerResStartRow;
        uint tFollowerResEndRow;
        uint tFollowerNumRows = 0;

        if ( mFollowerGlobalDofTypes.size() > 0 )
        {
            tFollowerDofIndex    = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
            tFollowerResStartRow = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 0 );
            tFollowerResEndRow   = mSet->get_res_dof_assembly_map()( tFollowerDofIndex )( 0, 1 );
            tFollowerNumRows     = tFollowerResEndRow - tFollowerResStartRow + 1;
        }

        // get number of leader and follower rows
        uint tLeaderNumRows = tLeaderResEndRow - tLeaderResStartRow + 1;

        // get number of cols for jacobian
        uint tNumCols = mSet->get_jacobian().n_cols();

        // set size for analytical and FD jacobians
        aJacobian.set_size( tLeaderNumRows + tFollowerNumRows, tNumCols, 0.0 );
        aJacobianFD.set_size( tLeaderNumRows + tFollowerNumRows, tNumCols, 0.0 );

        // compute jacobian with IWG
        this->compute_jacobian( aWStar );

        // get the computed jacobian
        aJacobian( { 0, tLeaderNumRows - 1 }, { 0, tNumCols - 1 } ) =
                mSet->get_jacobian()( { tLeaderResStartRow, tLeaderResEndRow }, { 0, tNumCols - 1 } );
        if ( tFollowerNumRows > 0 )
        {
            aJacobian( { tLeaderNumRows, tLeaderNumRows + tFollowerNumRows - 1 }, { 0, tNumCols - 1 } ) =
                    mSet->get_jacobian()( { tFollowerResStartRow, tFollowerResEndRow }, { 0, tNumCols - 1 } );
        }

        // reset the jacobian
        mSet->get_jacobian().fill( 0.0 );

        // compute jacobian by FD
        this->compute_jacobian_FD( aWStar, aPerturbation, fem::FDScheme_Type::POINT_5, aUseAbsolutePerturbations );

        // get the computed jacobian
        aJacobianFD( { 0, tLeaderNumRows - 1 }, { 0, tNumCols - 1 } ) =
                mSet->get_jacobian()( { tLeaderResStartRow, tLeaderResEndRow }, { 0, tNumCols - 1 } );
        if ( tFollowerNumRows > 0 )
        {
            aJacobianFD( { tLeaderNumRows, tLeaderNumRows + tFollowerNumRows - 1 }, { 0, tNumCols - 1 } ) =
                    mSet->get_jacobian()( { tFollowerResStartRow, tFollowerResEndRow }, { 0, tNumCols - 1 } );
        }

        // check that matrices to compare have same size
        MORIS_ERROR(
                ( aJacobian.n_rows() == aJacobianFD.n_rows() ) && ( aJacobian.n_cols() == aJacobianFD.n_cols() ),
                "IWG::check_jacobian - matrices to check do not share same dimensions." );

        // define a boolean for check
        bool tCheckJacobian = true;

        // define a real for absolute difference
        real tAbsolute = 0.0;

        // define a real for relative difference
        real tRelative = 0.0;

        for ( uint iiJac = 0; iiJac < aJacobian.n_rows(); iiJac++ )
        {
            for ( uint jjJac = 0; jjJac < aJacobian.n_cols(); jjJac++ )
            {
                // get absolute difference
                tAbsolute = std::abs( aJacobian( iiJac, jjJac ) - aJacobianFD( iiJac, jjJac ) );

                // get relative difference
                tRelative = std::abs( ( aJacobianFD( iiJac, jjJac ) - aJacobian( iiJac, jjJac ) ) / aJacobianFD( iiJac, jjJac ) );

                // update check value
                tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) );

                // debug print
                if ( ( ( tAbsolute < aEpsilon ) || ( tRelative < aEpsilon ) ) == false )
                {
                    if ( aErrorPrint )
                    {
                        std::cout << "iiJac " << iiJac << " - jjJac " << jjJac << "\n";
                        std::cout << "aJacobian( iiJac, jjJac )   " << std::setprecision( 12 ) << aJacobian( iiJac, jjJac ) << "\n";
                        std::cout << "aJacobianFD( iiJac, jjJac ) " << std::setprecision( 12 ) << aJacobianFD( iiJac, jjJac ) << "\n";
                        std::cout << "Absolute difference " << tAbsolute << "\n";
                        std::cout << "Relative difference " << tRelative << "\n"
                                  << std::flush;
                    }
                }
            }
        }

        // return bool
        return tCheckJacobian;
    }

    //------------------------------------------------------------------------------

    // FIXME: This function needs to go, functionality will be integrated into the usual check jacobian function
    bool IWG::check_jacobian_multi_residual(
            real              aPerturbation,
            real              aEpsilon,
            real              aWStar,
            Matrix< DDRMat >& aJacobian,
            Matrix< DDRMat >& aJacobianFD,
            bool              aErrorPrint,
            bool              aMaxErrorPrint,
            moris::real       aFDtolerance,
            bool              aUseAbsolutePerturbations )
    {
        // check if FD tolerance is defined
        if ( aFDtolerance < 0.0 )
        {
            aFDtolerance = aEpsilon;
        }

        // get residual dof type index in set, start and end indices for residual dof type
        uint tNumResidualDofs     = mResidualDofType.size();
        uint tLeaderFirstDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderLastDofIndex  = mSet->get_dof_index_for_type( mResidualDofType( tNumResidualDofs - 1 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResStartRow   = mSet->get_res_dof_assembly_map()( tLeaderFirstDofIndex )( 0, 0 );
        uint tLeaderResEndRow     = mSet->get_res_dof_assembly_map()( tLeaderLastDofIndex )( 0, 1 );

        // get number of leader and follower rows
        uint tNumRows = tLeaderResEndRow - tLeaderResStartRow + 1;

        // get number of cols for jacobian
        uint tNumCols = mSet->get_jacobian().n_cols();

        // set size for analytical and FD jacobians
        aJacobian.set_size( tNumRows, tNumCols, 0.0 );
        aJacobianFD.set_size( tNumRows, tNumCols, 0.0 );

        // compute jacobian with IWG
        this->compute_jacobian( aWStar );

        // get the computed jacobian
        aJacobian = mSet->get_jacobian();

        // reset the jacobian
        mSet->get_jacobian().fill( 0.0 );

        // compute jacobian by FD
        // this->compute_jacobian_FD( aWStar, aPerturbation );
        {
            // storage residual value
            Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

            // get the FD scheme info
            Vector< Vector< real > > tFDScheme;
            fd_scheme( fem::FDScheme_Type::POINT_5, tFDScheme );
            uint tNumFDPoints = tFDScheme( 0 ).size();

            // get leader number of dof types
            uint tLeaderNumDofTypes = mRequestedLeaderGlobalDofTypes.size();

            // reset and evaluate the residual plus
            mSet->get_residual()( 0 ).fill( 0.0 );
            this->compute_residual( aWStar );
            Matrix< DDRMat > tResidual = mSet->get_residual()( 0 );

            // loop over the IWG dof types
            for ( uint iFI = 0; iFI < tLeaderNumDofTypes; iFI++ )
            {
                // init dof counter
                uint tDofCounter = 0;

                // get the dof type
                Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iFI );

                // get the index for the dof type
                sint tLeaderDepDofIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderFirstDofIndex )( tLeaderDepDofIndex, 0 );

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
                        if ( aUseAbsolutePerturbations )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // check that perturbation is not zero
                        if ( std::abs( tDeltaH ) < 1e-12 )
                        {
                            tDeltaH = aPerturbation;
                        }

                        // set starting point for FD
                        uint tStartPoint = 0;

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

                            // reset the residual
                            mSet->get_residual()( 0 ).fill( 0.0 );

                            // compute the residual
                            this->compute_residual( aWStar );

                            // assemble the Jacobian
                            mSet->get_jacobian()(
                                    { tLeaderResStartRow, tLeaderResEndRow },
                                    { tLeaderDepStartIndex + tDofCounter, tLeaderDepStartIndex + tDofCounter } ) +=
                                    tFDScheme( 1 )( iPoint ) * mSet->get_residual()( 0 ) / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                        // update dof counter
                        tDofCounter++;
                    }
                }
                // reset the coefficients values
                tFI->set_coeff( tCoeff );
            }

            // reset the value of the residual
            mSet->get_residual()( 0 ) = tResidualStore;
        }

        // get the computed jacobian
        aJacobianFD = mSet->get_jacobian();

        // check that matrices to compare have same size
        MORIS_ERROR(
                ( aJacobian.n_rows() == aJacobianFD.n_rows() ) && ( aJacobian.n_cols() == aJacobianFD.n_cols() ),
                "IWG::check_jacobian - matrices to check do not share same dimensions." );

        // define a boolean for check
        bool tCheckJacobian = true;

        // define a real for absolute difference
        real tAbsolute    = 0.0;
        real tMaxAbsolute = 0.0;

        // define a real for relative difference
        real tRelative    = 0.0;
        real tMaxRelative = 0.0;

        for ( uint iiJac = 0; iiJac < aJacobian.n_rows(); iiJac++ )
        {
            for ( uint jjJac = 0; jjJac < aJacobian.n_cols(); jjJac++ )
            {
                // get absolute difference
                tAbsolute    = std::abs( aJacobian( iiJac, jjJac ) - aJacobianFD( iiJac, jjJac ) );
                tMaxAbsolute = std::max( tAbsolute, tMaxAbsolute );

                // get relative difference
                tRelative = std::abs( ( aJacobianFD( iiJac, jjJac ) - aJacobian( iiJac, jjJac ) ) / aJacobianFD( iiJac, jjJac ) );
                if ( tAbsolute > aFDtolerance )
                {
                    tMaxRelative = std::max( tRelative, tMaxRelative );
                }

                // update check value
                tCheckJacobian = tCheckJacobian && ( ( tAbsolute < aFDtolerance ) || ( tRelative < aEpsilon ) );

                // debug print
                if ( ( ( tAbsolute < aFDtolerance ) || ( tRelative < aEpsilon ) ) == false )
                {
                    if ( aErrorPrint )
                    {
                        std::cout << "iiJac " << iiJac << " - jjJac " << jjJac << "\n"
                                  << std::flush;
                        std::cout << "aJacobian( iiJac, jjJac )   " << std::setprecision( 12 ) << aJacobian( iiJac, jjJac ) << "\n"
                                  << std::flush;
                        std::cout << "aJacobianFD( iiJac, jjJac ) " << std::setprecision( 12 ) << aJacobianFD( iiJac, jjJac ) << "\n"
                                  << std::flush;
                        std::cout << "Absolute difference " << tAbsolute << "\n"
                                  << std::flush;
                        std::cout << "Relative difference " << tRelative << "\n"
                                  << std::flush;
                    }
                }
            }
        }

        // print maximum difference
        if ( aMaxErrorPrint )
        {
            std::cout << "Maximum absolute difference " << tMaxAbsolute << "\n"
                      << std::flush;
            std::cout << "Maximum relative difference " << tMaxRelative << "\n"
                      << std::flush;
        }

        // return bool
        return tCheckJacobian;
    }

    //------------------------------------------------------------------------------

    real IWG::build_perturbation_size(
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

    real IWG::build_perturbation_size_relative(
            const real& aPerturbation,
            const real& aCoefficientToPerturb,
            const real& aMaxPerturbation,
            const real& aTolerance )
    {
        // check that maximum perturbation size is larger than tolerance
        MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                "IWG::build_perturbation_size_relative - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
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

    real IWG::build_perturbation_size_absolute(
            const real& aPerturbation,
            const real& aCoefficientToPerturb,
            const real& aMaxPerturbation,
            const real& aTolerance )
    {
        // check that maximum perturbation size is larger than tolerance
        MORIS_ASSERT( aMaxPerturbation >= aTolerance,
                "IWG::build_perturbation_size_absolute - maximum perturbation size is smaller than tolerance: max = %e  tol = %e\n",
                aMaxPerturbation,
                aTolerance );

        // determine actual tolerance (only useful when above assert inactive)
        const real tActualTol = std::min( aMaxPerturbation, aTolerance );

        // check that absolute value of perturbation is not smaller than tolerance
        // and not larger than maximum value
        return std::max( std::min( std::abs( aPerturbation ), aMaxPerturbation ), tActualTol );
    }

    //------------------------------------------------------------------------------

    real IWG::check_ig_coordinates_inside_ip_element(
            const real&          aPerturbation,
            const real&          aCoefficientToPerturb,
            const uint&          aSpatialDirection,
            fem::FDScheme_Type&  aUsedFDScheme,
            mtk::Leader_Follower aIsLeader )
    {
        // FIXME: only works for rectangular IP elements
        // FIXME: only works for forward, backward, central, not for higher as 5-point FD

        // get the IP element geometry interpolator
        Geometry_Interpolator const * tIPGI =
                mSet->get_field_interpolator_manager( aIsLeader )->get_IP_geometry_interpolator();

        // IP element max/min
        // get maximum values of coordinates of IP nodes
        real const tMaxIP = max( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );

        // get minimum values of coordinates of IP nodes
        real const tMinIP = min( tIPGI->get_space_coeff().get_column( aSpatialDirection ) );

        // get maximum possible perturbation
        real const tMaxPerturb = ( tMaxIP - tMinIP ) / 3.0;

        // compute the perturbation value
        real const tDeltaH = build_perturbation_size(
                aPerturbation,
                aCoefficientToPerturb,
                tMaxPerturb,
                mToleranceFD );

        // check that IG node coordinate is consistent with minimum and maximum IP coordinates
        MORIS_ASSERT(
                tMaxIP >= aCoefficientToPerturb - tDeltaH && tMinIP <= aCoefficientToPerturb + tDeltaH,
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

                // check for correctness of perturbation size for forward FD
                MORIS_ASSERT( tDeltaH < tMaxIP - aCoefficientToPerturb,
                        "ERROR: forward perturbation size exceeds limits of interpolation element:\n"
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
                // check for correctness of perturbation size for central FD
                MORIS_ASSERT(
                        tDeltaH < tMaxIP - aCoefficientToPerturb && tDeltaH < aCoefficientToPerturb - tMinIP,
                        "ERROR: central perturbation size exceed limits of interpolation element:\n"
                        "dim: %d  minIP: %e  maxIP: %e  cordIG: %e  maxPert: %e  delta: %e  precPert: %e\n.",
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

    void IWG::select_dRdp_FD_geometry_bulk(
            moris::real                   aWStar,
            moris::real                   aPerturbation,
            fem::FDScheme_Type            aFDSchemeType,
            Matrix< DDSMat >&             aGeoLocalAssembly,
            Vector< Matrix< IndexMat > >& aVertexIndices )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the GI for the IG element considered
        Geometry_Interpolator* tIGGI =
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
        Geometry_Interpolator* tIPGI =
                mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

        // get the residual dof type index in the set
        uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
        uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual =
                mSet->get_residual()( 0 )(
                        { tResDofAssemblyStart, tResDofAssemblyStop },
                        { 0, 0 } );

        // get number of leader GI bases and space dimensions
        uint tDerNumBases      = tIGGI->get_number_of_space_bases();
        uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

        // coefficients for dv type wrt which derivative is computed
        Matrix< DDRMat > tCoeff      = tIGGI->get_space_coeff();          // get nodal coordinates of integration element
        Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();    // get IP natural coordinate of integration nodes

        Matrix< DDRMat > tEvaluationPoint;    // get IG natural coordinates of quadrature points
        tIGGI->get_space_time( tEvaluationPoint );

        real tGPWeight = aWStar / tIGGI->det_J();    // reconstruct integration weight

        // init perturbation
        real tDeltaH = 0.0;

        // init FD scheme
        Vector< Vector< real > > tFDScheme;

        // loop over the spatial directions
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
        {
            // loop over the IG nodes
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // get the geometry pdv assembly index
                sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                // if pdv is active
                if ( tPdvAssemblyIndex != -1 )
                {
                    // provide adapted perturbation and FD scheme considering ip element boundaries
                    fem::FDScheme_Type tUsedFDSchemeType = aFDSchemeType;

                    // compute step size and change FD scheme if needed
                    tDeltaH = this->check_ig_coordinates_inside_ip_element(
                            aPerturbation,
                            tCoeff( iCoeffRow, iCoeffCol ),
                            iCoeffCol,
                            tUsedFDSchemeType );

                    // finalize FD scheme
                    fd_scheme( tUsedFDSchemeType, tFDScheme );
                    uint tNumFDPoints = tFDScheme( 0 ).size();

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed residual contribution to dRdp
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over point of FD scheme
                    for ( uint iPoint = tStartPoint; iPoint < tNumFDPoints; iPoint++ )
                    {
                        // reset the perturbed coefficients, i.e. the nodal coordinate of integration element
                        Matrix< DDRMat > tCoeffPert = tCoeff;

                        // perturb the coefficient, i.e. the nodal coordinate of integration element
                        tCoeffPert( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;

                        // setting the perturbed coefficients, i.e. the nodal coordinate of integration element
                        tIGGI->set_space_coeff( tCoeffPert );

                        // update natural coordinates of IG nodes in IP element
                        Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                        Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                        tIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                        Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                        tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                        tIGGI->set_space_param_coeff( tParamCoeffPert );

                        // set evaluation point for interpolators (FIs and GIs)
                        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // reset and evaluate the residual plus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        real tWStarPert = tGPWeight * tIGGI->det_J();
                        this->compute_residual( tWStarPert );

                        // evaluate dRdpGeo
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                }
            }
        }

        // reset the coefficients values
        tIGGI->set_space_coeff( tCoeff );
        tIGGI->set_space_param_coeff( tParamCoeff );
        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // add contribution of cluster measure to dRdp
        if ( mActiveCMEAFlag )
        {
            // add their contribution to dQIdp
            this->add_cluster_measure_dRdp_FD_geometry(
                    aWStar,
                    aPerturbation,
                    aFDSchemeType );
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::select_dRdp_FD_geometry_sideset(
            moris::real                   aWStar,
            moris::real                   aPerturbation,
            fem::FDScheme_Type            aFDSchemeType,
            Matrix< DDSMat >&             aGeoLocalAssembly,
            Vector< Matrix< IndexMat > >& aVertexIndices )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the GI for the IG element considered
        Geometry_Interpolator* tIGGI =
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
        Geometry_Interpolator* tIPGI =
                mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();

        // get the residual dof type index in the set
        uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
        uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual = mSet->get_residual()( 0 )(
                { tResDofAssemblyStart, tResDofAssemblyStop },
                { 0, 0 } );

        // store unperturbed xyz
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

        // init perturbation
        real tDeltaH = 0.0;

        // init FD scheme
        Vector< Vector< real > > tFDScheme;

        // get number of leader GI bases and space dimensions
        uint tDerNumBases      = tIGGI->get_number_of_space_bases();
        uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

        // loop over the spatial directions
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
        {
            // loop over the IG nodes
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // get the geometry pdv assembly index
                sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                // if pdv is active
                if ( tPdvAssemblyIndex != -1 )
                {
                    // provide adapted perturbation and FD scheme considering ip element boundaries
                    fem::FDScheme_Type tUsedFDSchemeType = aFDSchemeType;

                    // compute step size and change FD scheme if needed
                    tDeltaH = this->check_ig_coordinates_inside_ip_element(
                            aPerturbation,
                            tCoeff( iCoeffRow, iCoeffCol ),
                            iCoeffCol,
                            tUsedFDSchemeType );

                    // finalize FD scheme
                    fd_scheme( tUsedFDSchemeType, tFDScheme );
                    uint tNumFDPoints = tFDScheme( 0 ).size();

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed residual contribution to dRdp
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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
                        tIPGI->update_parametric_coordinates( tXCoords, tXiCoords );
                        Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                        tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();
                        tIGGI->set_space_param_coeff( tParamCoeffPert );

                        // set evaluation point for interpolators (FIs and GIs)
                        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

                        // reset the normal
                        Matrix< DDRMat > tNormalPert;
                        tIGGI->get_normal( tNormalPert );
                        this->set_normal( tNormalPert );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // reset and evaluate the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        real tWStarPert = tGPWeight * tIGGI->det_J();
                        this->compute_residual( tWStarPert );

                        // evaluate dRdpGeo
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                }
            }
        }

        // reset xyz values
        tIGGI->set_space_coeff( tCoeff );

        // reset local coordinates values
        tIGGI->set_space_param_coeff( tParamCoeff );

        // reset evaluation point
        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );

        // reset normal
        this->set_normal( tNormal );

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // add contribution of cluster measure to dRdp
        if ( mActiveCMEAFlag )
        {
            // add their contribution to dQIdp
            this->add_cluster_measure_dRdp_FD_geometry(
                    aWStar,
                    aPerturbation,
                    aFDSchemeType );
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::select_dRdp_FD_geometry_time_sideset(
            moris::real                   aWStar,
            moris::real                   aPerturbation,
            fem::FDScheme_Type            aFDSchemeType,
            Matrix< DDSMat >&             aGeoLocalAssembly,
            Vector< Matrix< IndexMat > >& aVertexIndices )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the GI for the IG element considered
        Geometry_Interpolator* tIGGI =
                mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator();
        Geometry_Interpolator* tIPGI =
                mSet->get_field_interpolator_manager()->get_IP_geometry_interpolator();
        Geometry_Interpolator* tIGGIPrevious =
                mSet->get_field_interpolator_manager_previous_time()->get_IG_geometry_interpolator();

        // get the residual dof type index in the set
        uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
        uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual = mSet->get_residual()( 0 )(
                { tResDofAssemblyStart, tResDofAssemblyStop },
                { 0, 0 } );

        // store unperturbed xyz
        Matrix< DDRMat > tCoeff = tIGGI->get_space_coeff();

        // store unperturbed local coordinates
        Matrix< DDRMat > tParamCoeff = tIGGI->get_space_param_coeff();

        // store unperturbed evaluation point
        Matrix< DDRMat > tEvaluationPoint;
        tIGGI->get_space_time( tEvaluationPoint );

        // store unperturbed evaluation point weight
        real tGPWeight = aWStar / tIGGI->det_J();

        // init perturbation
        real tDeltaH = 0.0;

        // init FD scheme
        Vector< Vector< real > > tFDScheme;

        // get number of leader GI bases and space dimensions
        uint tDerNumBases      = tIGGI->get_number_of_space_bases();
        uint tDerNumDimensions = tIPGI->get_number_of_space_dimensions();

        // loop over the spatial directions
        for ( uint iCoeffCol = 0; iCoeffCol < tDerNumDimensions; iCoeffCol++ )
        {
            // loop over the IG nodes
            for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
            {
                // get the geometry pdv assembly index
                sint tPdvAssemblyIndex = aGeoLocalAssembly( iCoeffRow, iCoeffCol );

                // if pdv is active
                if ( tPdvAssemblyIndex != -1 )
                {
                    // provide adapted perturbation and FD scheme considering ip element boundaries
                    fem::FDScheme_Type tUsedFDSchemeType = aFDSchemeType;

                    // compute step size and change FD scheme if needed
                    tDeltaH = this->check_ig_coordinates_inside_ip_element(
                            aPerturbation,
                            tCoeff( iCoeffRow, iCoeffCol ),
                            iCoeffCol,
                            tUsedFDSchemeType );

                    // finalize FD scheme
                    fd_scheme( tUsedFDSchemeType, tFDScheme );
                    uint tNumFDPoints = tFDScheme( 0 ).size();

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward add unperturbed contribution
                    if ( ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed residual contribution to dRdp
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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
                        tIGGIPrevious->set_space_coeff( tCoeffPert );

                        // update local coordinates
                        Matrix< DDRMat > tXCoords  = tCoeffPert.get_row( iCoeffRow );
                        Matrix< DDRMat > tXiCoords = tParamCoeff.get_row( iCoeffRow );

                        tIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                        Matrix< DDRMat > tParamCoeffPert     = tParamCoeff;
                        tParamCoeffPert.get_row( iCoeffRow ) = tXiCoords.matrix_data();

                        tIGGI->set_space_param_coeff( tParamCoeffPert );

                        tIGGIPrevious->set_space_param_coeff( tParamCoeffPert );

                        // set evaluation point for interpolators (FIs and GIs)
                        mSet->get_field_interpolator_manager()->    //
                                set_space_time_from_local_IG_point( tEvaluationPoint );

                        mSet->get_field_interpolator_manager_previous_time()->    //
                                set_space_time_from_local_IG_point( tEvaluationPoint );

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // reset and evaluate the residual plus
                        mSet->get_residual()( 0 ).fill( 0.0 );
                        real tWStarPert = tGPWeight * tIGGI->det_J();
                        this->compute_residual( tWStarPert );

                        // evaluate dRdpGeo
                        mSet->get_drdpgeo()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                }
            }
        }

        // reset xyz values
        tIGGI->set_space_coeff( tCoeff );
        tIGGIPrevious->set_space_coeff( tCoeff );

        // reset local coordinates values
        tIGGI->set_space_param_coeff( tParamCoeff );
        tIGGIPrevious->set_space_param_coeff( tParamCoeff );

        // reset evaluation point
        mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tEvaluationPoint );
        mSet->get_field_interpolator_manager_previous_time()->set_space_time_from_local_IG_point( tEvaluationPoint );

        // reset the coefficients values
        tIGGI->set_space_coeff( tCoeff );
        tIGGIPrevious->set_space_coeff( tCoeff );

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // add contribution of cluster measure to dRdp
        if ( mActiveCMEAFlag )
        {
            // add their contribution to dQIdp
            this->add_cluster_measure_dRdp_FD_geometry(
                    aWStar,
                    aPerturbation,
                    aFDSchemeType );
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::select_dRdp_FD_geometry_double(
            moris::real                   aWStar,
            moris::real                   aPerturbation,
            fem::FDScheme_Type            aFDSchemeType,
            Matrix< DDSMat >&             aGeoLocalAssembly,
            Vector< Matrix< IndexMat > >& aVertexIndices )
    {
        // check if the element type is conformal or nonconformal
        bool const tIsNonconformal = ( mSet->get_element_type() == Element_Type::NONCONFORMAL_SIDESET );

        // get the leader GI for the IG and IP element considered
        Geometry_Interpolator* tLeaderIGGI = mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IG_geometry_interpolator();
        Geometry_Interpolator* tLeaderIPGI = mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IP_geometry_interpolator();

        // get the follower GI for the IG and IP element considered
        Geometry_Interpolator* tFollowerIGGI = mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIPGI = mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_IP_geometry_interpolator();

        // if nonconform check that mapping for nominal geometry is successful
        if ( tIsNonconformal && mResidualDofType( 0 )( 0 ) == MSI::Dof_Type::UX )
        {
            const Matrix< DDRMat > tRemappedFollowerCoords = this->remap_nonconformal_rays(
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) ),
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) ) );

            // check whether the remapping is successful
            if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
            {
                return;    // exit if remapping was not successful
            }
        }

        // unpack vertex indices
        Matrix< IndexMat >& aLeaderVertexIndices   = aVertexIndices( 0 );
        Matrix< IndexMat >& aFollowerVertexIndices = aVertexIndices( 1 );

        // build list of unique vertex indices
        Matrix< IndexMat > tUniqueVertexIndices;
        unique( join_horiz( aLeaderVertexIndices, aFollowerVertexIndices ), tUniqueVertexIndices );

        // initialize additional follower vertex flag
        // bool tHasUniqueFollowerVertices = false;

        // build table for unique vertices
        // 1. column: vertex index
        // 2. column: leader side=0; follower side=1; both sides = 2
        // 3. column: local index of vertex in leader or follower vertex index vector
        Matrix< DDSMat > tUniqueVertexTable( tUniqueVertexIndices.numel(), 4, -1 );

        // storage for unique geometry local assembly
        Matrix< DDSMat > tUniqueGeoLocalAssembly( tUniqueVertexIndices.numel(), aGeoLocalAssembly.n_cols() );

        for ( uint iUnique = 0; iUnique < tUniqueVertexIndices.numel(); ++iUnique )
        {
            // get the vertex index
            moris_index tVertexIndex = tUniqueVertexIndices( iUnique );

            bool tFoundUniqueNode = false;

            // check if the vertex is in the leader side
            for ( uint iLeader = 0; iLeader < aLeaderVertexIndices.numel(); iLeader++ )
            {
                if ( tVertexIndex == aLeaderVertexIndices( iLeader ) )
                {
                    tUniqueVertexTable( iUnique, 0 ) = tVertexIndex;
                    tUniqueVertexTable( iUnique, 1 ) = 0;
                    tUniqueVertexTable( iUnique, 2 ) = iLeader;

                    tUniqueGeoLocalAssembly.get_row( iUnique ) = aGeoLocalAssembly.get_row( iLeader );

                    tFoundUniqueNode = true;
                    break;
                }
            }

            // check if the vertex is in the follower side
            for ( uint iFollower = 0; iFollower < aFollowerVertexIndices.numel(); iFollower++ )
            {
                if ( tVertexIndex == aFollowerVertexIndices( iFollower ) )
                {
                    tUniqueVertexTable( iUnique, 0 ) = tVertexIndex;
                    tUniqueVertexTable( iUnique, 1 ) += 2;
                    tUniqueVertexTable( iUnique, 3 ) = iFollower;

                    // check if the vertex is exclusively in the follower side
                    if ( tUniqueVertexTable( iUnique, 1 ) == 1 )
                    {
                        // tHasUniqueFollowerVertices = true;

                        tUniqueGeoLocalAssembly.get_row( iUnique ) = aGeoLocalAssembly.get_row( iFollower + aLeaderVertexIndices.numel() );
                    }

                    tFoundUniqueNode = true;
                    break;
                }
            }
            MORIS_ERROR( tFoundUniqueNode,
                    "IWG::select_dRdp_FD_geometry_double - vertex index %d not found in leader or follower vertex indices.",
                    tVertexIndex );
        }

        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get requested geometry pdv types
        Vector< gen::PDV_Type > tRequestedGeoPdvType;
        mSet->get_ig_unique_dv_types_for_set( tRequestedGeoPdvType );

        // get the leader residual dof type index in the set
        uint const tLeaderResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint const tLeaderResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tLeaderResDofIndex )( 0, 0 );
        uint const tLeaderResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tLeaderResDofIndex )( 0, 1 );

        // get the follower residual dof type index in the set (if exists)
        sint const tFollowerResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );

        uint tFollowerResDofAssemblyStart = MORIS_UINT_MAX;
        uint tFollowerResDofAssemblyStop  = MORIS_UINT_MAX;

        if ( tFollowerResDofIndex != -1 )
        {
            tFollowerResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 0 );
            tFollowerResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 1 );
        }

        // get GP weight
        real const tGPWeight = aWStar / tLeaderIGGI->det_J();

        // save leader physcial and parametric coordinates
        const Matrix< DDRMat > tLeaderCoeff      = tLeaderIGGI->get_space_coeff();
        const Matrix< DDRMat > tLeaderParamCoeff = tLeaderIGGI->get_space_param_coeff();

        Matrix< DDRMat > tLeaderEvaluationPoint;
        tLeaderIGGI->get_space_time( tLeaderEvaluationPoint );

        Matrix< DDRMat > tLeaderNormal;
        tLeaderIGGI->get_normal( tLeaderNormal );

        // save leader physcial and parametric coordinates
        const Matrix< DDRMat > tFollowerCoeff      = tFollowerIGGI->get_space_coeff();
        const Matrix< DDRMat > tFollowerParamCoeff = tFollowerIGGI->get_space_param_coeff();

        Matrix< DDRMat > tFollowerEvaluationPoint;
        tFollowerIGGI->get_space_time( tFollowerEvaluationPoint );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );

        this->compute_residual( aWStar );

        Matrix< DDRMat > const tLeaderResidual =
                mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } );

        Matrix< DDRMat > tFollowerResidual;
        if ( tFollowerResDofIndex != -1 )
        {
            tFollowerResidual =
                    mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } );
        }

        // get number of leader GI bases and space dimensions
        uint const tNumDimensions = tLeaderIPGI->get_number_of_space_dimensions();

        // initialize FD scheme
        Vector< Vector< real > > tFDScheme;

        // loop over unique IG nodes
        for ( uint iUniqueNode = 0; iUniqueNode < tUniqueVertexTable.n_rows(); iUniqueNode++ )
        {
            // find the node on the leader and/or follower side
            sint iLeaderNode   = tUniqueVertexTable( iUniqueNode, 1 ) == 0 or tUniqueVertexTable( iUniqueNode, 1 ) == 2
                                       ? tUniqueVertexTable( iUniqueNode, 2 )
                                       : -1;
            sint iFollowerNode = tUniqueVertexTable( iUniqueNode, 1 ) == 1 or tUniqueVertexTable( iUniqueNode, 1 ) == 2
                                       ? tUniqueVertexTable( iUniqueNode, 3 )
                                       : -1;

            // check that at least one of the nodes is valid
            MORIS_ASSERT( iLeaderNode != -1 or iFollowerNode != -1,
                    "IWG::select_dRdp_FD_geometry_double - unique node %d does not have a valid leader or follower node.",
                    iUniqueNode );

            // loop over the spatial directions
            for ( uint iSpatialDir = 0; iSpatialDir < tNumDimensions; iSpatialDir++ )
            {
                // reset valid perturbation flag
                bool tUsePerturbation = true;

                // get the geometry pdv assembly index
                sint const tPdvAssemblyIndex = tUniqueGeoLocalAssembly( iUniqueNode, iSpatialDir );

                if ( tPdvAssemblyIndex != -1 )
                {
                    // determine perturbation FD scheme and compute step size considering ip element boundaries
                    fem::FDScheme_Type tLeaderFDSchemeType   = aFDSchemeType;
                    fem::FDScheme_Type tFollowerFDSchemeType = aFDSchemeType;

                    real tLeaderDeltaH   = MORIS_REAL_MAX;
                    real tFollowerDeltaH = MORIS_REAL_MAX;

                    if ( iLeaderNode > -1 )
                    {
                        tLeaderDeltaH = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tLeaderCoeff( iLeaderNode, iSpatialDir ),
                                iSpatialDir,
                                tLeaderFDSchemeType,
                                mtk::Leader_Follower::LEADER );
                    }
                    else
                    {
                        tLeaderFDSchemeType = fem::FDScheme_Type::END_FD_SCHEME;
                    }

                    if ( iFollowerNode > -1 )
                    {
                        tFollowerDeltaH = this->check_ig_coordinates_inside_ip_element(
                                aPerturbation,
                                tFollowerCoeff( iFollowerNode, iSpatialDir ),
                                iSpatialDir,
                                tFollowerFDSchemeType,
                                mtk::Leader_Follower::FOLLOWER );
                    }
                    else
                    {
                        tFollowerFDSchemeType = fem::FDScheme_Type::END_FD_SCHEME;
                    }

                    // finalize FD scheme
                    real tDeltaH = std::min( tLeaderDeltaH, tFollowerDeltaH );

                    // check if conflict in FD schemes
                    if ( ( tLeaderFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD &&    //
                                 tFollowerFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD )
                            or                                                                  //
                            ( tLeaderFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD &&    //
                                    tFollowerFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // invalid perturbation; skip FD for this node and direction
                        MORIS_LOG_INFO(
                                "IWG::select_dRdp_FD_geometry_double - unique node has conflicting FD schemes" );
                        continue;
                    }

                    // use the simpler of both schemes which corresponds to one with the smaller enum value
                    fem::FDScheme_Type tUsedFDSchemeType = std::min( tLeaderFDSchemeType, tFollowerFDSchemeType );

                    fd_scheme( tUsedFDSchemeType, tFDScheme );
                    uint const tNumFDPoints = tFDScheme( 0 ).size();

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // create local vector for derivatives for residual
                    Matrix< DDRMat > tLeaderDrDpgeo( tLeaderResDofAssemblyStop - tLeaderResDofAssemblyStart + 1, 1, 0.0 );

                    Matrix< DDRMat > tFollowerDrDpgeo;
                    if ( tFollowerResDofIndex != -1 )
                    {
                        tFollowerDrDpgeo.set_size( tFollowerResDofAssemblyStop - tFollowerResDofAssemblyStart + 1, 1, 0.0 );
                    }

                    // if backward or forward add unperturbed contribution
                    if ( ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( tUsedFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        tLeaderDrDpgeo +=
                                tFDScheme( 1 )( 0 ) * tLeaderResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // add unperturbed follower residual contribution to dRdp
                        if ( tFollowerResDofIndex != -1 )
                        {
                            tFollowerDrDpgeo +=
                                    tFDScheme( 1 )( 0 ) * tFollowerResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }

                        // skip first point in FD
                        tStartPoint = 1;
                    }

                    // loop over point of FD scheme
                    for ( uint iFDPoint = tStartPoint; iFDPoint < tNumFDPoints; iFDPoint++ )
                    {
                        // check that perturbation is valid
                        if ( !tUsePerturbation )
                        {
                            continue;
                        }

                        // perturb the coefficients of the leader or follower node
                        if ( iLeaderNode > -1 )
                        {
                            // perturb physcial coordinate of leader node and update in IGGI
                            Matrix< DDRMat > tLeaderCoeffPert = tLeaderCoeff;
                            tLeaderCoeffPert( iLeaderNode, iSpatialDir ) += tFDScheme( 0 )( iFDPoint ) * tDeltaH;
                            tLeaderIGGI->set_space_coeff( tLeaderCoeffPert );

                            // get parameteric coordinates of perturbed leader node in the IPGI
                            Matrix< DDRMat > tXCoords  = tLeaderCoeffPert.get_row( iLeaderNode );
                            Matrix< DDRMat > tXiCoords = tLeaderParamCoeff.get_row( iLeaderNode );
                            tLeaderIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                            // update parametric coordinates of leader node in IGGI
                            Matrix< DDRMat > tLeaderParamCoeffPert       = tLeaderParamCoeff;
                            tLeaderParamCoeffPert.get_row( iLeaderNode ) = tXiCoords.matrix_data();
                            tLeaderIGGI->set_space_param_coeff( tLeaderParamCoeffPert );

                            // update physical and parametric coordinates of quadrature point
                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                                    ->set_space_time_from_local_IG_point( tLeaderEvaluationPoint );

                            // update normal of the leader side
                            Matrix< DDRMat > tNormalPert;
                            tLeaderIGGI->get_normal( tNormalPert );
                            this->set_normal( tNormalPert );
                        }

                        if ( iFollowerNode > -1 )
                        {
                            // perturb physcial coordinate of follower node and update in IGGI
                            Matrix< DDRMat > tFollowerCoeffPert = tFollowerCoeff;
                            tFollowerCoeffPert( iFollowerNode, iSpatialDir ) += tFDScheme( 0 )( iFDPoint ) * tDeltaH;
                            tFollowerIGGI->set_space_coeff( tFollowerCoeffPert );

                            // get parameteric coordinates of perturbed follower node in the IPGI
                            Matrix< DDRMat > tXCoords  = tFollowerCoeffPert.get_row( iFollowerNode );
                            Matrix< DDRMat > tXiCoords = tFollowerParamCoeff.get_row( iFollowerNode );
                            tFollowerIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

                            // update parametric coordinates of follower node in IGGI
                            Matrix< DDRMat > tFollowerParamCoeffPert         = tFollowerParamCoeff;
                            tFollowerParamCoeffPert.get_row( iFollowerNode ) = tXiCoords.matrix_data();
                            tFollowerIGGI->set_space_param_coeff( tFollowerParamCoeffPert );

                            // update physical and parametric coordinates of quadrature point
                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                                    ->set_space_time_from_local_IG_point( tFollowerEvaluationPoint );
                        }

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // if nonconform check that mapping for perturbed geometry is successful
                        if ( tIsNonconformal && mResidualDofType( 0 )( 0 ) == MSI::Dof_Type::UX )
                        {
                            const Matrix< DDRMat > tRemappedFollowerCoords = this->remap_nonconformal_rays(
                                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) ),
                                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) ) );

                            // check whether the remapping is successful
                            if ( std::abs( tRemappedFollowerCoords( 0 ) ) > 1 )
                            {
                                tUsePerturbation = false;
                                continue;
                            }
                        }

                        // compute residual for perturbed geometry
                        if ( tUsePerturbation )
                        {
                            // reset and evaluate the residual plus
                            mSet->get_residual()( 0 ).fill( 0.0 );
                            real tWStarPert = tGPWeight * tLeaderIGGI->det_J();

                            // compute residual for perturbed geometry
                            this->compute_residual( tWStarPert );

                            // evaluate dLeaderRdpGeo
                            tLeaderDrDpgeo +=
                                    tFDScheme( 1 )( iFDPoint ) *                                                                          //
                                    mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );

                            // evaluate dFollowerRdpGeo (not needed in the nonconformal case)
                            if ( tFollowerResDofIndex != -1 )
                            {
                                tFollowerDrDpgeo +=
                                        tFDScheme( 1 )( iFDPoint ) *                                                                              //
                                        mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } ) /    //
                                        ( tFDScheme( 2 )( 0 ) * tDeltaH );
                            }
                        }

                        // reset the perturbed coefficients
                        if ( iLeaderNode > -1 )
                        {
                            // reset the coefficients in leader IGGI
                            tLeaderIGGI->set_space_coeff( tLeaderCoeff );
                            tLeaderIGGI->set_space_param_coeff( tLeaderParamCoeff );

                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->    //
                                    set_space_time_from_local_IG_point( tLeaderEvaluationPoint );

                            this->set_normal( tLeaderNormal );
                        }

                        if ( iFollowerNode > -1 )
                        {
                            // reset the coefficients in follower IGGI
                            tFollowerIGGI->set_space_coeff( tFollowerCoeff );
                            tFollowerIGGI->set_space_param_coeff( tFollowerParamCoeff );

                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                                    set_space_time_from_local_IG_point( tFollowerEvaluationPoint );
                        }
                    }
                    // if perturbation was valid, update resdual derivatives
                    if ( tUsePerturbation )
                    {
                        mSet->get_drdpgeo()(
                                { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                                { tPdvAssemblyIndex, tPdvAssemblyIndex } ) += tLeaderDrDpgeo.matrix_data();

                        // add unperturbed follower residual contribution to dRdp
                        if ( tFollowerResDofIndex != -1 )
                        {
                            mSet->get_drdpgeo()(
                                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                    { tPdvAssemblyIndex, tPdvAssemblyIndex } ) += tFollowerDrDpgeo.matrix_data();
                        }
                    }
                }
            }
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // add contribution of cluster measure to dRdp
        if ( mActiveCMEAFlag )
        {
            // add their contribution to dQIdp
            this->add_cluster_measure_dRdp_FD_geometry_double(
                    aWStar,
                    aPerturbation,
                    aFDSchemeType );
        }

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry_double - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::add_cluster_measure_dRdp_FD_geometry(
            moris::real        aWStar,
            moris::real        aPerturbation,
            fem::FDScheme_Type aFDSchemeType )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the residual dof type index in the set
        uint tResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
        uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual =
                mSet->get_residual()( 0 )(
                        { tResDofAssemblyStart, tResDofAssemblyStop },
                        { 0, 0 } );

        // initialize perturbation of cluster measure
        real tDeltaCM = 0.0;

        // init FD scheme
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // loop over the cluster measures
        for ( auto const& [ _, tClusterMeasure ] : mCluster->get_cluster_measures() )
        {
            // evaluate the perturbation of cluster measure
            tDeltaCM = this->build_perturbation_size(
                    aPerturbation,
                    tClusterMeasure->val()( 0 ),
                    std::max( tClusterMeasure->val()( 0 ), mToleranceFD ),
                    mToleranceFD );

            // number of pdv to assemble
            uint tNumPdvToAssemble = tClusterMeasure->dMEAdPDV().numel() - 1;

            // set starting point for FD
            uint tStartPoint = 0;

            // if backward or forward add unperturbed contribution
            if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                    ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
            {
                // add unperturbed residual contribution to dRdp
                mSet->get_drdpgeo()(
                        { tResDofAssemblyStart, tResDofAssemblyStop },
                        { 0, tNumPdvToAssemble } ) +=
                        tFDScheme( 1 )( 0 ) * tResidual *    //
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

                // reset and evaluate the residual plus
                mSet->get_residual()( 0 ).fill( 0.0 );
                this->compute_residual( aWStar );

                // evaluate dRdpGeo
                mSet->get_drdpgeo()(
                        { tResDofAssemblyStart, tResDofAssemblyStop },
                        { 0, tNumPdvToAssemble } ) +=
                        tFDScheme( 1 )( iPoint ) *                                                                //
                        mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) *    //
                        tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                // reset cluster measures
                mCluster->reset_cluster_measure();
            }
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::add_cluster_measure_dRdp_FD_geometry_double(
            moris::real        aWStar,
            moris::real        aPerturbation,
            fem::FDScheme_Type aFDSchemeType )
    {
        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the residual dof type index in the set
        uint tLeaderDofIndex            = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
        uint tLeaderResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

        // get the follower residual dof type index in the set
        uint tFollowerResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        uint tFollowerResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 0 );
        uint tFollowerResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tLeaderResidual =
                mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } );
        Matrix< DDRMat > tFollowerResidual =
                mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } );

        // initialize perturbation on cluster measure
        real tDeltaCM = 0.0;

        // initialize FD scheme
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // loop over the cluster measures
        for ( auto const& [ _, tClusterMeasure ] : mCluster->get_cluster_measures() )
        {
            // evaluate the perturbation of cluster measure
            tDeltaCM = this->build_perturbation_size(
                    aPerturbation,
                    tClusterMeasure->val()( 0 ),
                    std::max( tClusterMeasure->val()( 0 ), mToleranceFD ),
                    mToleranceFD );

            // get end pdv index
            uint tEndPdvIndex = tClusterMeasure->dMEAdPDV().numel() - 1;

            // set starting point for FD
            uint tStartPoint = 0;

            // if backward or forward add unperturbed contribution
            if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                    ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
            {
                // add unperturbed leader residual contribution to dRdp
                mSet->get_drdpgeo()(
                        { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                        { 0, tEndPdvIndex } ) +=
                        tFDScheme( 1 )( 0 ) * tLeaderResidual * tClusterMeasure->dMEAdPDV() /    //
                        ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                // add unperturbed follower residual contribution to dRdp
                mSet->get_drdpgeo()(
                        { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                        { 0, tEndPdvIndex } ) +=
                        tFDScheme( 1 )( 0 ) * tFollowerResidual * tClusterMeasure->dMEAdPDV() /    //
                        ( tFDScheme( 2 )( 0 ) * tDeltaCM );

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

                // reset and evaluate the residual plus
                mSet->get_residual()( 0 ).fill( 0.0 );
                this->compute_residual( aWStar );

                // evaluate dLeaderRdpGeo
                mSet->get_drdpgeo()(
                        { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                        { 0, tEndPdvIndex } ) +=
                        tFDScheme( 1 )( iPoint ) *                                                                          //
                        mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } )    //
                        * tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                // evaluate dFollowerRdpGeo
                mSet->get_drdpgeo()(
                        { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                        { 0, tEndPdvIndex } ) +=
                        tFDScheme( 1 )( iPoint ) *                                                                                //
                        mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } ) *    //
                        tClusterMeasure->dMEAdPDV() / ( tFDScheme( 2 )( 0 ) * tDeltaCM );

                // reset cluster measures
                mCluster->reset_cluster_measure();
            }
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpgeo() ),
                "IWG::compute_dRdp_FD_geometry - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::select_dRdp_FD_material(
            moris::real        aWStar,
            moris::real        aPerturbation,
            fem::FDScheme_Type aFDSchemeType )
    {
        // get the requested ip pdv types
        Vector< Vector< gen::PDV_Type > > tRequestedPdvTypes;
        mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

        // get number of requested dv types
        uint tNumDvType = tRequestedPdvTypes.size();

        // skip remainder if no dv types active on element
        if ( tNumDvType == 0 )
        {
            return;
        }

        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // get the residual dof type index in the set
        uint tResDofIndex = mSet->get_dof_index_for_type(
                mResidualDofType( 0 )( 0 ),
                mtk::Leader_Follower::LEADER );
        uint tResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
        uint tResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );
        this->compute_residual( aWStar );
        Matrix< DDRMat > tResidual =
                mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } );

        // init perturbation
        real tDeltaH = 0.0;

        // loop over the dv types associated with a FI
        for ( uint iFI = 0; iFI < tNumDvType; iFI++ )
        {
            // get dv index
            sint tDvDepIndex = mSet->get_dv_index_for_type(
                    tRequestedPdvTypes( iFI )( 0 ),
                    mtk::Leader_Follower::LEADER );

            // get the FI for the dv type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

            // get number of leader FI bases and fields
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // init coeff counter
            uint tCoeffCounter = 0;

            // loop over the coefficient column
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over the coefficient row
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // get mat pdv index
                    uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward FD, add unperturbed residual contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        mSet->get_drdpmat()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

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

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // evaluate dRdpMat
                        mSet->get_drdpmat()(
                                { tResDofAssemblyStart, tResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                //
                                mSet->get_residual()( 0 )( { tResDofAssemblyStart, tResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );
                    }
                    // update coefficient counter
                    tCoeffCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpmat() ),
                "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

    void IWG::select_dRdp_FD_material_double(
            moris::real        aWStar,
            moris::real        aPerturbation,
            fem::FDScheme_Type aFDSchemeType )
    {
        // get the requested ip pdv types
        Vector< Vector< gen::PDV_Type > > tRequestedPdvTypes;
        mSet->get_ip_dv_types_for_set( tRequestedPdvTypes );

        // get number of requested dv types
        uint tNumDvType = tRequestedPdvTypes.size();

        // skip remainder if no dv types active on element
        if ( tNumDvType == 0 )
        {
            return;
        }

        // storage residual value
        Matrix< DDRMat > tResidualStore = mSet->get_residual()( 0 );

        // get the FD scheme info
        Vector< Vector< real > > tFDScheme;
        fd_scheme( aFDSchemeType, tFDScheme );
        uint tNumFDPoints = tFDScheme( 0 ).size();

        // get the leader residual dof type index in the set
        uint tLeaderResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
        uint tLeaderResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tLeaderResDofIndex )( 0, 0 );
        uint tLeaderResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tLeaderResDofIndex )( 0, 1 );

        // get the follower residual dof type index in the set
        sint tFollowerResDofIndex         = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::FOLLOWER );
        uint tFollowerResDofAssemblyStart = MORIS_UINT_MAX;
        uint tFollowerResDofAssemblyStop  = MORIS_UINT_MAX;

        // check that follower side has residual dof type
        if ( tFollowerResDofIndex != -1 )
        {
            tFollowerResDofAssemblyStart = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 0 );
            tFollowerResDofAssemblyStop  = mSet->get_res_dof_assembly_map()( tFollowerResDofIndex )( 0, 1 );
        }

        // reset, evaluate and store the residual for unperturbed case
        mSet->get_residual()( 0 ).fill( 0.0 );

        this->compute_residual( aWStar );

        Matrix< DDRMat > tLeaderResidual =
                mSet->get_residual()( 0 )(
                        { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                        { 0, 0 } );

        Matrix< DDRMat > tFollowerResidual;
        if ( tFollowerResDofIndex != -1 )
        {
            tFollowerResidual = mSet->get_residual()( 0 )(
                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                    { 0, 0 } );
        }

        // init perturbation
        real tDeltaH = 0.0;

        // loop over the leader dv types associated with a FI
        for ( uint iFI = 0; iFI < tNumDvType; iFI++ )
        {
            // get dv index
            sint tDvDepIndex = mSet->get_dv_index_for_type(
                    tRequestedPdvTypes( iFI )( 0 ),
                    mtk::Leader_Follower::LEADER );

            // get the FI for the dv type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

            // get number of leader FI bases and fields
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // init coeff counter
            uint tCoeffCounter = 0;

            // loop over the coefficient column
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over the coefficient row
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // get mat pdv index
                    uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward FD, add unperturbed residual contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        mSet->get_drdpmat()(
                                { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tLeaderResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // add unperturbed leader residual contribution to dRdp
                        if ( tFollowerResDofIndex != -1 )
                        {
                            mSet->get_drdpmat()(
                                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tFollowerResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }

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

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // assemble dRLeaderdpMat
                        mSet->get_drdpmat()(
                                { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                            //
                                mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // assemble dRFollowerdpMat
                        if ( tFollowerResDofIndex != -1 )
                        {
                            mSet->get_drdpmat()(
                                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) *                                                                                //
                                    mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                    // update coefficient counter
                    tCoeffCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        // get the requested ip pdv types for the follower side
        mSet->get_ip_dv_types_for_set( tRequestedPdvTypes, mtk::Leader_Follower::FOLLOWER );

        // loop over the follower dv types associated with a FI
        for ( uint iFI = 0; iFI < tRequestedPdvTypes.size(); iFI++ )
        {
            // get dv index
            sint tDvDepIndex = mSet->get_dv_index_for_type(
                    tRequestedPdvTypes( iFI )( 0 ),
                    mtk::Leader_Follower::FOLLOWER );

            // get the FI for the dv type
            Field_Interpolator* tFI =
                    mFollowerFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) );

            // get number of leader FI bases and fields
            uint tDerNumBases  = tFI->get_number_of_space_time_bases();
            uint tDerNumFields = tFI->get_number_of_fields();

            // coefficients for dof type wrt which derivative is computed
            Matrix< DDRMat > tCoeff = tFI->get_coeff();

            // reset properties, CM and SP for IWG
            this->reset_eval_flags();

            // reset and evaluate the residual
            mSet->get_residual()( 0 ).fill( 0.0 );

            this->compute_residual( aWStar );

            Matrix< DDRMat > tLeaderResidual =
                    mSet->get_residual()( 0 )(
                            { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                            { 0, 0 } );

            Matrix< DDRMat > tFollowerResidual;
            if ( tFollowerResDofIndex != -1 )
            {
                tFollowerResidual =
                        mSet->get_residual()( 0 )(
                                { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                { 0, 0 } );
            }

            // init coeff counter
            uint tCoeffCounter = 0;

            // loop over the coefficient column
            for ( uint iCoeffCol = 0; iCoeffCol < tDerNumFields; iCoeffCol++ )
            {
                // loop over the coefficient row
                for ( uint iCoeffRow = 0; iCoeffRow < tDerNumBases; iCoeffRow++ )
                {
                    // compute the perturbation absolute value
                    tDeltaH = aPerturbation * tCoeff( iCoeffRow, iCoeffCol );

                    // check that perturbation is not zero
                    if ( std::abs( tDeltaH ) < 1e-12 )
                    {
                        tDeltaH = aPerturbation;
                    }

                    // get mat pdv index
                    uint tPdvIndex = mSet->get_mat_pdv_assembly_map()( tDvDepIndex )( 0, 0 ) + tCoeffCounter;

                    // set starting point for FD
                    uint tStartPoint = 0;

                    // if backward or forward FD, add unperturbed residual contribution
                    if ( ( aFDSchemeType == fem::FDScheme_Type::POINT_1_BACKWARD ) ||    //
                            ( aFDSchemeType == fem::FDScheme_Type::POINT_1_FORWARD ) )
                    {
                        // add unperturbed leader residual contribution to dRdp
                        mSet->get_drdpmat()(
                                { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( 0 ) * tLeaderResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // add unperturbed leader residual contribution to dRdp
                        if ( tFollowerResDofIndex != -1 )
                        {
                            mSet->get_drdpmat()(
                                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( 0 ) * tFollowerResidual / ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }

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

                        // reset the residual
                        mSet->get_residual()( 0 ).fill( 0.0 );

                        // compute the residual
                        this->compute_residual( aWStar );

                        // assemble dRLeaderdpMat
                        mSet->get_drdpmat()(
                                { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop },
                                { tPdvIndex, tPdvIndex } ) +=
                                tFDScheme( 1 )( iPoint ) *                                                                            //
                                mSet->get_residual()( 0 )( { tLeaderResDofAssemblyStart, tLeaderResDofAssemblyStop }, { 0, 0 } ) /    //
                                ( tFDScheme( 2 )( 0 ) * tDeltaH );

                        // assemble dRFollowerdpMat
                        if ( tFollowerResDofIndex != -1 )
                        {
                            mSet->get_drdpmat()(
                                    { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop },
                                    { tPdvIndex, tPdvIndex } ) +=
                                    tFDScheme( 1 )( iPoint ) *                                                                                //
                                    mSet->get_residual()( 0 )( { tFollowerResDofAssemblyStart, tFollowerResDofAssemblyStop }, { 0, 0 } ) /    //
                                    ( tFDScheme( 2 )( 0 ) * tDeltaH );
                        }
                    }
                    // update coefficient counter
                    tCoeffCounter++;
                }
            }
            // reset the coefficients values
            tFI->set_coeff( tCoeff );
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpmat() ),
                "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
