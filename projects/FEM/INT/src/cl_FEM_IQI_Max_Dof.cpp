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
            // FIXME: convert these values to parameters
            real tReferenceValue = 0.0;
            real tExponent = 10.0;

            // get field interpolator for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // check if dof index was set (for the case of vector field)
            if( mMasterDofTypes( 0 ).size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            aQI = {{ std::pow( tFI->val()( mIQITypeIndex ) - tReferenceValue, tExponent ) }};
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_QI( moris::real aWStar )
        {
            // FIXME: convert these values to parameters
            real tReferenceValue = 0.0;
            real tExponent = 10.0;

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get field interpolator for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // check if dof index was set (for the case of vector field)
            if( mMasterDofTypes( 0 ).size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1, "IQI_Max_Dof::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += { aWStar *  std::pow( tFI->val()( mIQITypeIndex ) - tReferenceValue, tExponent ) };
        }

        //------------------------------------------------------------------------------

        void IQI_Max_Dof::compute_dQIdu( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR( false, "IQI_Max_Dof::compute_dQIdu( Matrix< DDRMat > & adQIdDof ) --- Not implemented, yet. \n" );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */
