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
#include "cl_FEM_Contact_Gap_Data.hpp"
#include "cl_FEM_Contact_Geometry.hpp"
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

#include <iomanip>
#include <utility>

namespace moris::fem
{
    //------------------------------------------------------------------------------

    // GapData is deliberately incomplete in IWG.hpp. Construct and destroy
    // its owning pointer here, keeping contact data out of generic IWG headers.
    IWG::IWG()  = default;
    IWG::~IWG() = default;

    //------------------------------------------------------------------------------

    void
    IWG::compute_lumped_matrix(
            Field_Interpolator* aFIRow,
            Field_Interpolator* aFIColumn,
            Lumped_Matrix&      aMatrix )
    {
        MORIS_ERROR( aFIRow != nullptr && aFIColumn != nullptr,
                "IWG::compute_lumped_matrix - Field interpolators must be provided." );

        const uint tNumSpace  = aFIRow->get_number_of_space_bases();
        const uint tNumTime   = aFIRow->get_number_of_time_bases();
        const uint tNumFields = aFIRow->get_number_of_fields();

        MORIS_ERROR( aFIColumn->get_number_of_space_bases() == tNumSpace
                             && aFIColumn->get_number_of_time_bases() == tNumTime
                             && aFIColumn->get_number_of_fields() == tNumFields,
                "IWG::compute_lumped_matrix - Field interpolator dimensions must match." );

        const Matrix< DDRMat >& tRowN    = aFIRow->N();
        const Matrix< DDRMat >& tColumnN = aFIColumn->N();
        aMatrix.set_size( tNumSpace, tNumTime, tNumFields );

        for ( uint iField = 0; iField < tNumFields; ++iField )
        {
            const uint tOffset = iField * tNumTime * tNumSpace;

            for ( uint jTime = 0; jTime < tNumTime; ++jTime )
            {
                // Sum only this field/time block. Multiplying by a row shape
                // value gives the corresponding row sum of N_row^T N_column.
                real tColumnSum = 0.0;

                for ( uint iSpace = 0; iSpace < tNumSpace; ++iSpace )
                {
                    tColumnSum += tColumnN( iField, tOffset + jTime * tNumSpace + iSpace );
                }

                for ( uint iTime = 0; iTime < tNumTime; ++iTime )
                {
                    for ( uint iSpace = 0; iSpace < tNumSpace; ++iSpace )
                    {
                        aMatrix.diagonal( iField, iTime, jTime, iSpace ) =
                                tRowN( iField, tOffset + iTime * tNumSpace + iSpace ) * tColumnSum;
                    }
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IWG::compute_lumped_matrices(
            Field_Interpolator*     aFICurrent,
            Field_Interpolator*     aFIPrevious,
            Lumped_Matrix&          aLumpCurrent,
            Lumped_Matrix&          aLumpPrevious,
            Matrix< DDRMat >&       aInitialConditionCoeff,
            const Matrix< DDRMat >& aInitialConditions,
            bool                    aInitialTime )
    {
        compute_lumped_matrix( aFICurrent, aFICurrent, aLumpCurrent );

        // Constant time interpolation uses the current operator on both sides.
        if ( aFICurrent->get_number_of_time_bases() == 1 )
        {
            aLumpPrevious = aLumpCurrent;
        }
        else
        {
            compute_lumped_matrix( aFICurrent, aFIPrevious, aLumpPrevious );
        }

        aInitialConditionCoeff.set_size( 0, 0 );

        if ( aInitialTime )
        {
            const uint tNumFields = aFICurrent->get_number_of_fields();
            MORIS_ERROR( aInitialConditions.numel() == tNumFields,
                    "IWG::compute_lumped_matrices - One initial value per field is required." );

            const uint tNumBases = aFICurrent->get_number_of_space_bases() * aFICurrent->get_number_of_time_bases();
            aInitialConditionCoeff.set_size( tNumBases, tNumFields );

            for ( uint iField = 0; iField < tNumFields; ++iField )
            {
                for ( uint iBasis = 0; iBasis < tNumBases; ++iBasis )
                {
                    aInitialConditionCoeff( iBasis, iField ) = aInitialConditions( iField );
                }
            }
        }
    }

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
            case fem::Element_Type::TIME_FINAL_SIDESET:
            case fem::Element_Type::TIME_BOUNDARY:
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
            case mtk::Leader_Follower::FOLLOWER:
            {
                mFollowerPreviousFIManager = aFieldInterpolatorManager;
                break;
            }

            default:
            {
                MORIS_ERROR( false, "IWG::set_field_interpolator_manager - can only be leader or follower" );
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

                        // Invalidate cached gap geometry after each perturbation so the
                        // contact residual remaps using the perturbed displacement field.
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

                        // Invalidate cached gap geometry after each perturbation so the
                        // contact residual remaps using the perturbed displacement field.
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
    Matrix< DDRMat > IWG::remap_nonconformal_rays_undeformed_geometry(
            Field_Interpolator_Manager* aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager* aFollowerFieldInterpolatorManager ) const
    {
        return Contact_Geometry::remap_nonconformal_rays_undeformed_geometry(
                aLeaderFieldInterpolatorManager,
                aFollowerFieldInterpolatorManager );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IWG::remap_nonconformal_rays(
            const bool                     aUseDeformedGeometryForGap,
            const Vector< MSI::Dof_Type >& aDisplDofTypes,
            Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
            Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
            std::unique_ptr< GapData >&    aGapData ) const
    {
        return Contact_Geometry::remap_nonconformal_rays(
                aUseDeformedGeometryForGap,
                aDisplDofTypes,
                aLeaderFieldInterpolatorManager,
                aFollowerFieldInterpolatorManager,
                aLeaderPreviousFieldInterpolatorManager,
                aFollowerPreviousFieldInterpolatorManager,
                aGapData );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IWG::remap_nonconformal_rays_deformed_geometry(
            const Vector< MSI::Dof_Type >& aDisplDofTypes,
            Field_Interpolator_Manager*    aLeaderFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerFieldInterpolatorManager,
            Field_Interpolator_Manager*    aLeaderPreviousFieldInterpolatorManager,
            Field_Interpolator_Manager*    aFollowerPreviousFieldInterpolatorManager,
            std::unique_ptr< GapData >&    aGapData ) const
    {
        return Contact_Geometry::remap_nonconformal_rays_deformed_geometry(
                aDisplDofTypes,
                aLeaderFieldInterpolatorManager,
                aFollowerFieldInterpolatorManager,
                aLeaderPreviousFieldInterpolatorManager,
                aFollowerPreviousFieldInterpolatorManager,
                aGapData );
    }

    //------------------------------------------------------------------------------

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

        Geometry_Interpolator* tLeaderPreviousIGGI   = nullptr;
        Geometry_Interpolator* tFollowerPreviousIGGI = nullptr;

        if ( mSet->get_time_continuity() )
        {
            tLeaderPreviousIGGI   = mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )->get_IG_geometry_interpolator();
            tFollowerPreviousIGGI = mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER )->get_IG_geometry_interpolator();
        }

        // if nonconform check that mapping for nominal geometry is successful
        const Vector< MSI::Dof_Type > tDisplDofTypes = mResidualDofType( 0 );
        if ( tIsNonconformal && tDisplDofTypes( 0 ) == MSI::Dof_Type::UX )
        {
            const Matrix< DDRMat > tRemappedFollowerCoords = Contact_Geometry::remap_nonconformal_rays(
                    mUseDeformedGeometryForGap,
                    tDisplDofTypes,
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER ),
                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER ),
                    mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER ),
                    mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER ),
                    mGapData,
                    false );

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

        // lock contact conditions
        if ( tIsNonconformal && mResidualDofType( 0 )( 0 ) == MSI::Dof_Type::UX )
        {
            mGapData->mIsLocked = true;
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

                            if ( mSet->get_time_continuity() )
                            {
                                tLeaderPreviousIGGI->set_space_coeff( tLeaderCoeffPert );
                                tLeaderPreviousIGGI->set_space_param_coeff( tLeaderParamCoeffPert );

                                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )
                                        ->set_space_time_from_local_IG_point( tLeaderEvaluationPoint );
                            }

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

                            if ( mSet->get_time_continuity() )
                            {
                                tFollowerPreviousIGGI->set_space_coeff( tFollowerCoeffPert );
                                tFollowerPreviousIGGI->set_space_param_coeff( tFollowerParamCoeffPert );

                                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER )
                                        ->set_space_time_from_local_IG_point( tFollowerEvaluationPoint );
                            }
                        }

                        // reset properties, CM and SP for IWG
                        this->reset_eval_flags();

                        // if nonconform check that mapping for perturbed geometry is successful
                        if ( tIsNonconformal && mResidualDofType( 0 )( 0 ) == MSI::Dof_Type::UX )
                        {
                            const Matrix< DDRMat > tRemappedFollowerCoords = Contact_Geometry::remap_nonconformal_rays(
                                    mUseDeformedGeometryForGap,
                                    tDisplDofTypes,
                                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER ),
                                    mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER ),
                                    mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER ),
                                    mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER ),
                                    mGapData,
                                    false );

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

                            if ( mSet->get_time_continuity() )
                            {
                                tLeaderPreviousIGGI->set_space_coeff( tLeaderCoeff );
                                tLeaderPreviousIGGI->set_space_param_coeff( tLeaderParamCoeff );

                                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::LEADER )->    //
                                        set_space_time_from_local_IG_point( tLeaderEvaluationPoint );
                            }

                            // reset normal
                            this->set_normal( tLeaderNormal );
                        }

                        if ( iFollowerNode > -1 )
                        {
                            // reset the coefficients in follower IGGI
                            tFollowerIGGI->set_space_coeff( tFollowerCoeff );

                            tFollowerIGGI->set_space_param_coeff( tFollowerParamCoeff );

                            mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->    //
                                    set_space_time_from_local_IG_point( tFollowerEvaluationPoint );

                            if ( mSet->get_time_continuity() )
                            {
                                tFollowerPreviousIGGI->set_space_coeff( tFollowerCoeff );
                                tFollowerPreviousIGGI->set_space_param_coeff( tFollowerParamCoeff );

                                mSet->get_field_interpolator_manager_previous_time( mtk::Leader_Follower::FOLLOWER )->    //
                                        set_space_time_from_local_IG_point( tFollowerEvaluationPoint );
                            }
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

        // unlock contact conditions in gap data
        if ( tIsNonconformal && mResidualDofType( 0 )( 0 ) == MSI::Dof_Type::UX )
        {
            mGapData->mIsLocked = false;
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

        Vector< Vector< gen::PDV_Type > > tFollowerPdvTypes;
        mSet->get_ip_dv_types_for_set( tFollowerPdvTypes, mtk::Leader_Follower::FOLLOWER );

        // A material ADV can occur only on the follower of a contact pair.
        if ( tNumDvType == 0 && tFollowerPdvTypes.size() == 0 )
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

            // The same material design coefficient controls constitutive
            // response on both slabs. Its partial includes previous traction,
            // while displacement/history coefficients remain fixed.
            Field_Interpolator*    tPreviousFI    = get_time_continuity()
                                                          ? mLeaderPreviousFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) )
                                                          : nullptr;
            const Matrix< DDRMat > tPreviousCoeff = tPreviousFI != nullptr
                                                          ? tPreviousFI->get_coeff()
                                                          : Matrix< DDRMat >();

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
                        if ( tPreviousFI != nullptr )
                        {
                            Matrix< DDRMat > tPreviousPerturbed = tPreviousCoeff;
                            tPreviousPerturbed( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;
                            tPreviousFI->set_coeff( tPreviousPerturbed );
                        }

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
            if ( tPreviousFI != nullptr )
            {
                tPreviousFI->set_coeff( tPreviousCoeff );
            }
            this->reset_eval_flags();
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

            // The same material design coefficient controls constitutive
            // response on both slabs. Its partial includes previous traction,
            // while displacement/history coefficients remain fixed.
            Field_Interpolator*    tPreviousFI    = get_time_continuity()
                                                          ? mFollowerPreviousFIManager->get_field_interpolators_for_type( tRequestedPdvTypes( iFI )( 0 ) )
                                                          : nullptr;
            const Matrix< DDRMat > tPreviousCoeff = tPreviousFI != nullptr
                                                          ? tPreviousFI->get_coeff()
                                                          : Matrix< DDRMat >();

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
                        if ( tPreviousFI != nullptr )
                        {
                            Matrix< DDRMat > tPreviousPerturbed = tPreviousCoeff;
                            tPreviousPerturbed( iCoeffRow, iCoeffCol ) += tFDScheme( 0 )( iPoint ) * tDeltaH;
                            tPreviousFI->set_coeff( tPreviousPerturbed );
                        }

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
            if ( tPreviousFI != nullptr )
            {
                tPreviousFI->set_coeff( tPreviousCoeff );
            }
            this->reset_eval_flags();
        }

        // reset the value of the residual
        mSet->get_residual()( 0 ) = tResidualStore;

        // check for nan, infinity
        MORIS_ASSERT( isfinite( mSet->get_drdpmat() ),
                "IWG::compute_dRdp_FD_material - dRdp contains NAN or INF, exiting!" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
