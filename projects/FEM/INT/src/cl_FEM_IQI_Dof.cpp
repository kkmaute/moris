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
        void IQI_Dof::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator for a given dof type
            Field_Interpolator * tFI
            = mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // check if dof indx was set (for the case of vectorial field)
            if( mMasterDofTypes( 0 ).size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI.matrix_data() = { tFI->val()( mIQITypeIndex ) };
        }

//------------------------------------------------------------------------------
        void IQI_Dof::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR( false, "IQI_Dof::compute_dQIdDof - not implemented." );
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



