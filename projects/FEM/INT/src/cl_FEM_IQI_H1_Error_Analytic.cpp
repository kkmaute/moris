/*
 * cl_FEM_IQI_H1_Error_Analytic.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IQI_H1_Error_Analytic::IQI_H1_Error_Analytic()
        {
            // set IQI type
            mIQIType = vis::Output_Type::H1_ERROR_ANALYTIC;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "H1Check" ] = IQI_Property_Type::H1_CHECK;
        }

//------------------------------------------------------------------------------
        void IQI_H1_Error_Analytic::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get field interpolator
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mMasterDofTypes( 0 )( 0 ) );

            // get property index
            uint tH1CheckIndex = static_cast< uint >( IQI_Property_Type::H1_CHECK );

            // evaluate the QI
            aQI = trans( tFI->gradx( 1 ) - mMasterProp( tH1CheckIndex )->val() ) * ( tFI->gradx( 1 ) - mMasterProp( tH1CheckIndex )->val() ) ;
        }

//------------------------------------------------------------------------------
        void IQI_H1_Error_Analytic::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR( false, "IQI_H1_Error_Analytic::compute_dQIdDof - Not implemented." );
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



