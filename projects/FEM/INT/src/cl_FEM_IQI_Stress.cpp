/*
 * cl_FEM_IQI_Stress.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Stress.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IQI_Stress::IQI_Stress()
        {
            // set IQI type
            mIQIType = vis::Output_Type::STRAIN_ENERGY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;
        }
//------------------------------------------------------------------------------
        void IQI_Stress::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            aQI.resize(1,1);

            // evaluate the QI
            aQI(0,0) = mMasterCM( tElastLinIsoIndex )->flux()(mIQITypeIndex);
        }

//------------------------------------------------------------------------------
        void IQI_Stress::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            MORIS_ERROR(0,"Derivates of stress wrt dof not implemented");
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



