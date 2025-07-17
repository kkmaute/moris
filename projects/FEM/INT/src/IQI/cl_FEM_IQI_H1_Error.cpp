/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_H1_Error.cpp
 *
 */

#include "cl_FEM_IQI_H1_Error.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_vectorize.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_H1_Error::IQI_H1_Error()
    {
        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::H1_ERROR;

        // set the property pointer cell size
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "L2_Reference" ]  = static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE );
        mPropertyMap[ "H1S_Reference" ] = static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE );
    }

    //------------------------------------------------------------------------------

    void IQI_H1_Error::initialize()
    {
        if ( !mIsInitialized )
        {
            // size of parameter list
            uint tParamSize = mParameters.size();

            // check for proper size of constant function parameters
            MORIS_ERROR( tParamSize == 2 || tParamSize == 3,
                    "IQI_H1_Error::initialize - either 2 or 3 constant parameters need to be set." );

            // get weights for L2 and H1 semi-norm contributions
            mL2Weight  = mParameters( 0 )( 0 );
            mH1SWeight = mParameters( 1 )( 0 );

            // extract parameter whether to skip computing dQIdu
            if ( tParamSize > 2 )
            {
                mSkipComputeDQIDU = mParameters( 2 )( 0 ) > 0;
            }

            // check mQuantityDofType is defined
            MORIS_ERROR( mQuantityDofType.size() > 0,
                    "IQI_H1_Error::initialize - dof_quantity parameter needs to be defined." );

            // set initialize flag to true
            mIsInitialized = true;
        }
    }

    //------------------------------------------------------------------------------

    void IQI_H1_Error::compute_QI( Matrix< DDRMat > &aQI )
    {
        // initialize if needed
        this->initialize();

        // get field interpolator
        Field_Interpolator *tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // initialize QI
        aQI.fill( 0.0 );

        // L2 contributions
        if ( mL2Weight > 0.0 )
        {
            MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::L2_REFERENCE_VALUE ),
                    "IQI_H1_Error::compute_QI - no weight for L2 contribution provided.\n" );

            const std::shared_ptr< Property > &tPropL2Value =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

            // compute difference between dof value and reference value
            auto tL2error = tFI->val() - tPropL2Value->val();

            // compute L2 error
            aQI += mL2Weight * trans( tL2error ) * tL2error;
        }

        // H1 semi-norm contribution
        if ( mH1SWeight > 0.0 )
        {
            MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::H1S_REFERENCE_VALUE ),
                    "IQI_H1_Error::compute_QI - no weight for H1 semi-norm contribution provided.\n" );

            const std::shared_ptr< Property > &tPropH1SValue =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

            // compute difference between dof spatial gradient and reference value and flatten it
            Matrix< DDRMat > tH1Serror = vectorize( tFI->gradx( 1 ) - tPropH1SValue->val() );

            // compute H1 semi-norm error
            aQI += mH1SWeight * trans( tH1Serror ) * tH1Serror;
        }
    }

    //------------------------------------------------------------------------------

    void IQI_H1_Error::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        Matrix< DDRMat > tQI( 1, 1 );

        this->compute_QI( tQI );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
    }

    //------------------------------------------------------------------------------

    void IQI_H1_Error::compute_dQIdu( real aWStar )
    {
        // check whether to compute derivatives
        if ( mSkipComputeDQIDU )
        {
            return;
        }

        // initialize if needed
        this->initialize();

        // get field interpolator
        Field_Interpolator *tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get sub-vector
            auto tRes = mSet->get_residual()( tQIIndex )(
                    { tLeaderDepStartIndex, tLeaderDepStopIndex } );

            // L2 contributions
            if ( mL2Weight > 0.0 )
            {
                MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::L2_REFERENCE_VALUE ),
                        "IQI_H1_Error::compute_QI - no weight for L2 contribution provided.\n" );

                const std::shared_ptr< Property > &tPropL2Value =
                        mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

                // compute difference between dof value and reference value
                auto tL2error = tFI->val() - tPropL2Value->val();

                // derivative with respect to quantity dof type
                if ( tDofType( 0 ) == mQuantityDofType( 0 ) )
                {
                    // compute dQIdDof
                    tRes += 2.0 * aWStar * mL2Weight * tFI->N_trans() * tL2error;
                }

                // if L2 reference value depends on dof type
                if ( tPropL2Value->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdof
                    tRes -= 2.0 * aWStar * mL2Weight * trans( tPropL2Value->dPropdDOF( tDofType ) ) * tL2error;
                }
            }

            // H1 semi-norm contribution
            if ( mH1SWeight > 0.0 )
            {
                MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::H1S_REFERENCE_VALUE ),
                        "IQI_H1_Error::compute_QI - no weight for H1 semi-norm contribution provided.\n" );

                const std::shared_ptr< Property > &tPropH1SValue =
                        mLeaderProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

                // compute difference between dof spatial gradient and reference value
                auto tH1Serror = tFI->gradx( 1 ) - tPropH1SValue->val();

                // derivative with respect to quantity dof type
                if ( tDofType( 0 ) == mQuantityDofType( 0 ) )
                {
                    // compute dQIdDof
                    tRes += 2.0 * aWStar * mH1SWeight * trans( tFI->dnNdxn( 1 ) ) * tH1Serror;
                }

                // if H1S reference value depends on dof type
                if ( tPropH1SValue->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdDof
                    tRes -= 2.0 * aWStar * mH1SWeight * trans( tPropH1SValue->dPropdDOF( tDofType ) ) * tH1Serror;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void IQI_H1_Error::compute_dQIdu(
            Vector< MSI::Dof_Type > &aDofType,
            Matrix< DDRMat >        &adQIdu )
    {
        // check whether to compute derivatives
        if ( mSkipComputeDQIDU )
        {
            return;
        }

        // initialize if needed
        this->initialize();

        // get field interpolator
        Field_Interpolator *tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // L2 contributions
        if ( mL2Weight > 0.0 )
        {
            MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::L2_REFERENCE_VALUE ),
                    "IQI_H1_Error::compute_QI - no weight for L2 contribution provided.\n" );

            const std::shared_ptr< Property > &tPropL2Value =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_REFERENCE_VALUE ) );

            // compute difference between dof value and reference value
            auto tL2error = tFI->val() - tPropL2Value->val();

            // derivative with respect to quantity dof type
            if ( aDofType( 0 ) == mQuantityDofType( 0 ) )
            {
                // compute dQIdDof
                adQIdu += 2.0 * mL2Weight * tFI->N_trans() * tL2error;
            }

            // if L2 reference value depends on dof type
            if ( tPropL2Value->check_dof_dependency( aDofType ) )
            {
                // compute dQIdof
                adQIdu -= 2.0 * mL2Weight * tPropL2Value->dPropdDOF( aDofType ) * tL2error;
            }
        }

        // H1 semi-norm contribution
        if ( mH1SWeight > 0.0 )
        {
            MORIS_ASSERT( mLeaderProp( (uint)IQI_Property_Type::H1S_REFERENCE_VALUE ),
                    "IQI_H1_Error::compute_QI - no weight for H1 semi-norm contribution provided.\n" );

            const std::shared_ptr< Property > &tPropH1SValue =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::H1S_REFERENCE_VALUE ) );

            // compute difference between dof spatial gradient and reference value
            auto tH1Serror = tFI->gradx( 1 ) - tPropH1SValue->val();

            // derivative with respect to quantity dof type
            if ( aDofType( 0 ) == mQuantityDofType( 0 ) )
            {
                // compute dQIdDof
                adQIdu += 2.0 * mH1SWeight * trans( tFI->dnNdxn( 1 ) ) * tH1Serror;
            }

            // if H1S reference value depends on dof type
            if ( tPropH1SValue->check_dof_dependency( aDofType ) )
            {
                // compute dQIdDof
                adQIdu -= 2.0 * mH1SWeight * tPropH1SValue->dPropdDOF( aDofType ) * tH1Serror;
            }
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
