/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_ALM_Dof.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_ALM_Dof.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IQI_ALM_Dof::IQI_ALM_Dof()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "LagrangeMultiplier" ] = static_cast< uint >( IQI_Property_Type::LAGRANGE_MULTIPLIER );
            mPropertyMap[ "PenaltyFactor" ]      = static_cast< uint >( IQI_Property_Type::PENALTY_FACTOR );

            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::ALM_DOF;
        }

        //------------------------------------------------------------------------------

        void
        IQI_ALM_Dof::initialize()
        {
            if ( !mIsInitialized )
            {
                // size of parameter list
                uint tParamSize = mParameters.size();

                // check for proper size of constant function parameters
                MORIS_ERROR( tParamSize <= 2,
                        "IQI_ALM_Dof::initialize - either 1 or 2 constant parameters need to be set." );

                mRefValue = mParameters( 0 )( 0 );

                // shift parameter
                mShift = 1.0;

                if ( tParamSize > 1 )
                {
                    mShift = mParameters( 1 )( 0 );
                }

                // check mQuantityDofType is defined
                MORIS_ERROR( mQuantityDofType.size() > 0,
                        "IQI_ALM_Dof::initialize - dof_quantity parameter needs to be defined." );

                // check if dof index was set (for the case of vector field)
                if ( mQuantityDofType.size() > 1 )
                {
                    MORIS_ERROR( mIQITypeIndex != -1, "IQI_ALM_Dof::compute_QI - mIQITypeIndex not set." );
                }
                else
                {
                    mIQITypeIndex = 0;
                }

                // check that Lagrange multiplier and penalty factor are defined
                const std::shared_ptr< Property >& tPropLagrangeMultiplier =
                        mMasterProp( static_cast< uint >( IQI_Property_Type::LAGRANGE_MULTIPLIER ) );

                const std::shared_ptr< Property >& tPropPenaltyFactor =
                        mMasterProp( static_cast< uint >( IQI_Property_Type::PENALTY_FACTOR ) );

                MORIS_ERROR( tPropLagrangeMultiplier && tPropPenaltyFactor,
                        "IQI_ALM_Dof::initialize - properties for Lagrange multiplier and penalty factor need to be defined." );

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_ALM_Dof::compute_QI( Matrix< DDRMat >& aQI )
        {
            // initialize if needed
            this->initialize();

            // get Lagrange multiplier and penalty factor
            const std::shared_ptr< Property >& tPropLagrangeMultiplier =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::LAGRANGE_MULTIPLIER ) );

            const std::shared_ptr< Property >& tPropPenaltyFactor =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PENALTY_FACTOR ) );

            // get field interpolator for a given dof type
            Field_Interpolator* tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // evaluate Lagrange multiplier and penalty factor
            real tLagrangeMultplier = tPropLagrangeMultiplier->val()( 0 );
            real tPenaltyFactor     = tPropPenaltyFactor->val()( 0 );

            // evaluate constraint
            real tConstraint = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;

            // evaluate gplus
            real tGplus = std::max( tConstraint, -( tLagrangeMultplier / tPenaltyFactor ) );

            // evaluate integrand of IQI
            aQI = { { tLagrangeMultplier * tGplus + 0.5 * tPenaltyFactor * std::pow( tGplus, 2 ) } };

            MORIS_ASSERT( isfinite( aQI ),
                    "IQI_ALM_Dof::compute_QI - QI is nan, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI_ALM_Dof::compute_QI( real aWStar )
        {
            // define matrix to store temporarily integrand of IQI
            Matrix< DDRMat > tQI;

            // evaluate integrand
            this->compute_QI( tQI );

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // compute IQI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );

            MORIS_ASSERT( isfinite( mSet->get_QI()( tQIIndex ) ),
                    "IQI_ALM_Dof::compute_QI - QI is nan, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IQI_ALM_Dof::compute_dQIdu( real aWStar )
        {
            // get Lagrange multiplier and penalty factor
            const std::shared_ptr< Property >& tPropLagrangeMultiplier =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::LAGRANGE_MULTIPLIER ) );

            const std::shared_ptr< Property >& tPropPenaltyFactor =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PENALTY_FACTOR ) );

            // get field interpolator for max dof type
            Field_Interpolator* tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

                // if derivative dof type is max dof type
                if ( tDofType( 0 ) == mQuantityDofType( 0 ) )
                {
                    // evaluate Lagrange multiplier and penalty factor
                    real tLagrangeMultplier = tPropLagrangeMultiplier->val()( 0 );
                    real tPenaltyFactor     = tPropPenaltyFactor->val()( 0 );

                    // evaluate constraint
                    real tContraint = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;

                    // evaluate gplus
                    if ( tContraint >= -( tLagrangeMultplier / tPenaltyFactor ) )
                    {
                        // build selection matrix for vector dofs
                        uint tNumVecFieldComps = tFIMaxDof->val().numel();

                        Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                        tSelect( mIQITypeIndex, 0 ) = 1.0;

                        // evaluate derivative of integrand wrt constraint
                        real tdQI = tLagrangeMultplier + tPenaltyFactor * tContraint;

                        // evaluate derivative of IQI wrt dof
                        mSet->get_residual()( tQIIndex )(
                                { tMasterDepStartIndex, tMasterDepStopIndex }, { 0, 0 } ) +=
                                aWStar * ( tdQI * tFIMaxDof->N_trans() * tSelect / mRefValue );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void
        IQI_ALM_Dof::compute_dQIdu(
                moris::Cell< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&             adQIdu )
        {
            // get Lagrange multiplier and penalty factor
            const std::shared_ptr< Property >& tPropLagrangeMultiplier =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::LAGRANGE_MULTIPLIER ) );

            const std::shared_ptr< Property >& tPropPenaltyFactor =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PENALTY_FACTOR ) );

            // get field interpolator for max dof type
            Field_Interpolator* tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // if derivative dof type is max dof type
            if ( aDofType( 0 ) == mQuantityDofType( 0 ) )
            {
                // evaluate Lagrange multiplier and penalty factor
                real tLagrangeMultplier = tPropLagrangeMultiplier->val()( 0 );
                real tPenaltyFactor     = tPropPenaltyFactor->val()( 0 );

                // evaluate constraint
                real tContraint = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;

                // evaluate gplus
                if ( tContraint >= -( tLagrangeMultplier / tPenaltyFactor ) )
                {
                    // build selection matrix for vector dofs
                    uint tNumVecFieldComps = tFIMaxDof->val().numel();

                    Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                    tSelect( mIQITypeIndex, 0 ) = 1.0;

                    // evaluate derivative of integrand wrt constraint
                    real tdQI = tLagrangeMultplier + tPenaltyFactor * tContraint;

                    // evaluate derivative of IQI wrt dof
                    adQIdu = tdQI * tFIMaxDof->N_trans() * tSelect / mRefValue;
                }
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
