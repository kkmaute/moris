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

        IQI_Max_Dof::IQI_Max_Dof(){}

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator for a given dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // unpack constant parameters
            MORIS_ERROR( mParameters.size() == 2 || mParameters.size() == 3,
                    "IQI_Max_Dof::compute_QI - constant parameters were not set." );
            real tRefValue = mParameters( 0 )( 0 );
            real tExponent = mParameters( 1 )( 0 );
            real tShift = 1.0;
            if( mParameters.size() > 2 )
            {
                tShift = mParameters( 2 )( 0 );
            }

            // check if dof index was set (for the case of vector field)
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI = {{ std::pow( ( tFIMaxDof->val()( mIQITypeIndex ) / tRefValue ) - tShift, tExponent ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator for a given dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // unpack constant parameters
            MORIS_ERROR( mParameters.size() == 2 || mParameters.size() == 3,
                    "IQI_Max_Dof::compute_QI - constant parameters were not set." );
            real tRefValue = mParameters( 0 )( 0 );
            real tExponent = mParameters( 1 )( 0 );
            real tShift = 1.0;
            if( mParameters.size() > 2 )
            {
                tShift = mParameters( 2 )( 0 );
            }

            // check if dof index was set (for the case of vector field)
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    std::pow( ( tFIMaxDof->val()( mIQITypeIndex ) / tRefValue ) - tShift, tExponent ) );
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_dQIdu( real aWStar )
        {
            // get field interpolator for max dof type
            Field_Interpolator * tFIMaxDof =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // unpack constant parameters
            MORIS_ERROR( mParameters.size() == 2 || mParameters.size() == 3,
                    "IQI_Max_Dof::compute_QI - constant parameters were not set." );
            real tRefValue = mParameters( 0 )( 0 );
            real tExponent = mParameters( 1 )( 0 );
            real tShift = 1.0;
            if( mParameters.size() > 2 )
            {
                tShift = mParameters( 2 )( 0 );
            }

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
                    // check if dof index was set (for the case of vector field)
                    if( mQuantityDofType.size() > 1 )
                    {
                        MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
                    }
                    else
                    {
                        mIQITypeIndex = 0;
                    }

                    // build selection matrix
                    uint tNumVecFieldComps = tFIMaxDof->val().numel();
                    Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );
                    tSelect( mIQITypeIndex, 0 ) = 1.0;

                    // compute dQIdDof
                    real tdQI = std::pow( ( tFIMaxDof->val()( mIQITypeIndex ) / tRefValue ) - tShift, tExponent - 1.0 );
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * ( tExponent * tdQI * tFIMaxDof->N_trans() * tSelect / tRefValue );
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

            // unpack constant parameters
            MORIS_ERROR( mParameters.size() == 2 || mParameters.size() == 3,
                    "IQI_Max_Dof::compute_QI - constant parameters were not set." );
            real tRefValue = mParameters( 0 )( 0 );
            real tExponent = mParameters( 1 )( 0 );
            real tShift = 1.0;
            if( mParameters.size() > 2 )
            {
                tShift = mParameters( 2 )( 0 );
            }

            // if derivative dof type is max dof type
            if( aDofType( 0 ) == mQuantityDofType( 0 ) )
            {
                // check if dof index was set (for the case of vector field)
                if( mQuantityDofType.size() > 1 )
                {
                    MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
                }
                else
                {
                    mIQITypeIndex = 0;
                }

                // build selection matrix
                uint tNumVecFieldComps = tFIMaxDof->val().numel();
                Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );
                tSelect( mIQITypeIndex, 0 ) = 1.0;

                // compute dQIdDof
                real tdQI = std::pow( ( tFIMaxDof->val()( mIQITypeIndex ) / tRefValue ) - tShift, tExponent - 1.0 );
                adQIdu = tExponent * tdQI * trans( tFIMaxDof->N() ) * tSelect / tRefValue;
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
