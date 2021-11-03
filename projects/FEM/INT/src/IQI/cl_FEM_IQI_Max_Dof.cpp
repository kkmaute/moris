/*
 * cl_FEM_IQI_Max_Dof.cpp
 *
 *  Created on: Jul 10, 2020
 *      Author: wunsch
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Max_Dof.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void IQI_Max_Dof::initialize()
        {
            if ( ! mIsInitialized )
            {
                // size of parameter list
                uint tParamSize = mParameters.size();

                // check for proper size of constant function parameters
                MORIS_ERROR( tParamSize >= 2 && tParamSize <= 4,
                        "IQI_Max_Dof::initialize - either 2, 3 or 4 constant parameters need to be set." );

                if( tParamSize > 2 )
                {
                    mShift = mParameters( 2 )( 0 );
                }

                mRefValue = mParameters( 0 )( 0 );
                mExponent = mParameters( 1 )( 0 );

                // shift parameter
                mShift = 1.0;

                if( tParamSize > 2 )
                {
                    mShift = mParameters( 2 )( 0 );
                }

                // sign parameter
                mSign = 0;

                if( tParamSize > 3 )
                {
                    mSign = (sint) mParameters( 3 )( 0 );
                }

                // check mQuantityDofType is defined
                MORIS_ERROR( mQuantityDofType.size() > 0,
                        "IQI_Max_Dof::initialize - dof_quantity parameter needs to be defined." );

                // check if dof index was set (for the case of vector field)
                if( mQuantityDofType.size() > 1 )
                {
                    MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
                }
                else
                {
                    mIQITypeIndex = 0;
                }

                // set initialize flag to true
                mIsInitialized = true;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_QI( Matrix< DDRMat > & aQI )
        {
            // initialize if needed
            this->initialize();

            // get field interpolator for a given dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // evaluate the QI
            real tArgument = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;

            switch ( mSign )
            {
                case 1:
                    tArgument = tArgument > 0 ? tArgument : 0.0;
                    break;
                case -1:
                    tArgument = tArgument < 0 ? tArgument : 0.0;
                    break;
            }

            const std::shared_ptr< Property > & tPropWeight =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::WEIGHT ) );

            // check if the weight property has been stipulated
            if ( tPropWeight != nullptr )
            {
                aQI = {{ tPropWeight->val() * std::pow( tArgument, mExponent ) }};
            }
            else
            {
                aQI = {{ std::pow( tArgument, mExponent ) }};
            }

            MORIS_ASSERT( isfinite( aQI ),
                    "IQI_Max_Dof::compute_QI - QI is nan, exiting!");
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_QI( real aWStar )
        {
            // initialize if needed
            this->initialize();

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator for a given dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // evaluate the QI
            real tArgument = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;

            switch ( mSign )
            {
                case 1:
                    tArgument = tArgument > 0 ? tArgument : 0.0;
                    break;
                case -1:
                    tArgument = tArgument < 0 ? tArgument : 0.0;
                    break;
            }

            const std::shared_ptr< Property > & tPropWeight =
            mMasterProp( static_cast< uint >( IQI_Property_Type::WEIGHT ) );

            // check if the weight property has been stipulated
            if ( tPropWeight != nullptr )
            {
                mSet->get_QI()( tQIIndex ) += tPropWeight->val() * aWStar * (
                    std::pow( tArgument, mExponent ) );
            }
            else
            {
                mSet->get_QI()( tQIIndex ) += aWStar * (
                    std::pow( tArgument, mExponent ) );
            }

            MORIS_ASSERT( isfinite( mSet->get_QI()( tQIIndex ) ),
                    "IQI_Max_Dof::compute_QI - QI is nan, exiting!");
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_dQIdu( real aWStar )
        {
            // get field interpolator for max dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

                // if derivative dof type is max dof type
                if( tDofType( 0 ) == mQuantityDofType( 0 ) )
                {
                    // build selection matrix
                    uint tNumVecFieldComps = tFIMaxDof->val().numel();

                    Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                    tSelect( mIQITypeIndex, 0 ) = 1.0;

                    // compute dQIdDof
                    real tArgument = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue  - mShift;

                    switch ( mSign )
                    {
                        case 1:
                            tArgument = tArgument > 0 ? tArgument : 0.0;
                            break;
                        case -1:
                            tArgument = tArgument < 0 ? tArgument : 0.0;
                            break;
                    }

                    real tdQI = std::pow( tArgument, mExponent - 1.0 );

                    const std::shared_ptr< Property > & tPropWeight =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::WEIGHT ) );

                    // check if the weight property has been stipulated
                    if ( tPropWeight != nullptr )
                    {
                        mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },{ 0, 0 } ) += 
                            tPropWeight->val()(0) * aWStar * ( mExponent * tdQI * tFIMaxDof->N_trans() * tSelect / mRefValue );

                        // check dof dependency on the property
                        if ( tPropWeight->check_dof_dependency( tDofType ) )
                        {
                            mSet->get_residual()( tQIIndex )(
                                { tMasterDepStartIndex, tMasterDepStopIndex },{ 0, 0 } ) += 
                                tPropWeight->dPropdDOF( tDofType ) * aWStar * (std::pow( tArgument, mExponent ) );
                        }

                    }
                    else
                    {
                        mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex }, { 0, 0 } ) += 
                            aWStar * ( mExponent * tdQI * tFIMaxDof->N_trans() * tSelect / mRefValue );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get field interpolator for max dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // if derivative dof type is max dof type
            if( aDofType( 0 ) == mQuantityDofType( 0 ) )
            {
                // build selection matrix
                uint tNumVecFieldComps = tFIMaxDof->val().numel();

                Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                tSelect( mIQITypeIndex, 0 ) = 1.0;

                // compute dQIdDof
                real tArgument = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue  - mShift;

                switch ( mSign )
                {
                    case 1:
                        tArgument = tArgument > 0 ? tArgument : 0.0;
                        break;
                    case -1:
                        tArgument = tArgument < 0 ? tArgument : 0.0;
                        break;
                }

                real tdQI = std::pow( tArgument, mExponent - 1.0 );

                const std::shared_ptr< Property > & tPropWeight =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::WEIGHT ) );

                // check if the weight property has been stipulated
                if ( tPropWeight != nullptr )
                {
                    adQIdu = tPropWeight->val()(0) * mExponent * tdQI * trans( tFIMaxDof->N() ) * tSelect / mRefValue;

                    // check dof dependency on the property
                    if (tPropWeight->check_dof_dependency( aDofType ))
                    {
                        adQIdu +=  tPropWeight->dPropdDOF(aDofType) * (std::pow(tArgument, mExponent));
                    }
                }
                else
                {
                    adQIdu = mExponent * tdQI * trans( tFIMaxDof->N() ) * tSelect / mRefValue;
                }
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
