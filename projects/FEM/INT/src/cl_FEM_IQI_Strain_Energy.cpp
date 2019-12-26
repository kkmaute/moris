/*
 * cl_FEM_IQI_Strain_energy.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include "cl_FEM_IQI_Strain_Energy.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IQI_Strain_Energy::IQI_Strain_Energy()
        {
            // set IQI type
            mIQIType = vis::Output_Type::STRAIN_ENERGY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IQI_Constitutive_Type::ELAST_LIN_ISO;
        }
//------------------------------------------------------------------------------
        void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // evaluate the QI
            aQI = trans( mMasterCM( tElastLinIsoIndex )->flux() ) * mMasterCM( tElastLinIsoIndex )->strain();
        }

//        void IQI_strain_Energy::compute_dQIdDof()
//        {
//            // get indices for properties, CM and SP
//            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
//
//            // compute dQIdDof for indirect dof dependencies
//            uint tNumDofDep = mRequestedMasterGlobalDofTypes.size();
//            for( uint iDof = 0; iDof < tNumDofDep; iDof++ )
//            {
//                // get the treated dof type
//                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDof );
//
//                // if constitutive model has dependency on the dof type
//                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
//                {
//                    // compute dQIdDof
//                	// FIXME how do we assemble this...
//                    adQIdDof = trans( mMasterCM( tElastLinIsoIndex )->dStressdDOF( tDofType ) ) * mMasterCM( tDiffLinIsoIndex )->strain( )
//                             + trans( mMasterCM( tElastLinIsoIndex )->stress() ) * mMasterCM( tDiffLinIsoIndex )->dStraindDOF( tDofType );
//                }
//            }
//        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



