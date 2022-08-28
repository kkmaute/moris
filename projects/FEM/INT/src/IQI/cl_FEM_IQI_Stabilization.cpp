/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Stabilization.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Stabilization.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Stabilization::IQI_Stabilization()
        {
            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IQI_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "Stabilization" ] = static_cast< uint >( IQI_Stabilization_Type::STABILIZATION );
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

            // check if index was set
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Stabilization::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI = {{ tSP->val()( mIQITypeIndex ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

            // check if index was set
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Stabilization::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += {{ aWStar * tSP->val()( mIQITypeIndex ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

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

                // Dof dependency
                if ( tSP != nullptr && tSP->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * trans( tSP->dSPdMasterDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Stabilization::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IQI_Stabilization_Type::STABILIZATION ) );

            // Dof dependency
            if ( tSP != nullptr && tSP->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = trans( tSP->dSPdMasterDOF( aDofType ) );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

