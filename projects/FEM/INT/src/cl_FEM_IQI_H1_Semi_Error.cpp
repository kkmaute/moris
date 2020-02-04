/*
 * cl_FEM_IQI_H1_Semi_Error.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_H1_Semi_Error.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IQI_H1_Semi_Error::IQI_H1_Semi_Error()
        {
            // set IQI type
            mIQIType = vis::Output_Type::H1_SEMI_ERROR;
        }

//------------------------------------------------------------------------------
        void IQI_H1_Semi_Error::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // evaluate the QI
            aQI = trans( tFI->gradx( 1 ) ) * tFI->gradx( 1 );
        }

//------------------------------------------------------------------------------
        void IQI_H1_Semi_Error::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR( false, "IQI_H1_Semi_Error::compute_dQIdDof - Not implemented." );
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



