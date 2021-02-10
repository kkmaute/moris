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
            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::H1_SEMI_ERROR;
        }
        
        //------------------------------------------------------------------------------

        void IQI_H1_Semi_Error::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // evaluate the QI
            aQI = trans( tFI->gradx( 1 ) ) * tFI->gradx( 1 );
        }
        
        //------------------------------------------------------------------------------

        void IQI_H1_Semi_Error::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( trans( tFI->gradx( 1 ) ) * tFI->gradx( 1 ) );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



