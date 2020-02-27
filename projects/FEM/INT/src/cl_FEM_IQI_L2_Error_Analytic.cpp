/*
 * cl_FEM_IQI_L2_Error_Analytic.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IQI_L2_Error_Analytic::IQI_L2_Error_Analytic()
        {
            // set IQI type
            mIQIType = vis::Output_Type::L2_ERROR_ANALYTIC;

            // set FEM IQI type
            mFEMIQIType = fem::IQI_Type::L2_ERROR_ANALYTIC;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "L2Check" ] = IQI_Property_Type::L2_CHECK;
        }

//------------------------------------------------------------------------------
        void IQI_L2_Error_Analytic::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // get property index
            uint tL2CheckIndex = static_cast< uint >( IQI_Property_Type::L2_CHECK );

            // evaluate the QI
            aQI = ( trans( tFI->val() - mMasterProp( tL2CheckIndex )->val() ) * ( tFI->val() - mMasterProp( tL2CheckIndex )->val() ) );
        }

//------------------------------------------------------------------------------
        void IQI_L2_Error_Analytic::compute_dQIdu( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR( false, "IQI_L2_Error_Analytic::compute_dQIdu - Not implemented." );
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



