/*
 * cl_FEM_IQI_Dof.cpp
 *
 *  Created on: Jan 23, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Dof.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Dof::IQI_Dof(){}

        //------------------------------------------------------------------------------

        void IQI_Dof::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // check if dof index was set (for the case of vector field)
            if( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tFI->val()( mIQITypeIndex ) );
        }

        //------------------------------------------------------------------------------

        void IQI_Dof::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

            // check if dof index was set (for the case of vector field)
            if( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 )
            {
                aQI = { tFI->val()( mIQITypeIndex ) };
            }
            else
            {
                aQI = tFI->val();
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
