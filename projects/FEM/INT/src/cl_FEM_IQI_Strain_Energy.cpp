/*
 * cl_FEM_IQI_Strain_energy.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
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

//------------------------------------------------------------------------------
        void IQI_Strain_Energy::compute_dQIdDof( Matrix< DDRMat > & adQIdDof )
        {
            // FIXME should not be done here
            // get number of dof coefficients
            uint tNumCoeff = 0;

            // get the requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes = this->get_requested_dof_types();

            // loop over the requested dof types
            for( uint Ik = 0; Ik < tRequestedDofTypes.size(); Ik++ )
            {
                // get the set index for dof type
                sint tDofIndex = mSet->get_dof_index_for_type( tRequestedDofTypes( Ik ), mtk::Master_Slave::MASTER );

                // if the dof type was set
                if( tDofIndex != -1 )
                {
                    // update number of coefficients
                    tNumCoeff += mMasterFIManager->get_field_interpolators_for_type( tRequestedDofTypes( Ik ) )
                                                 ->get_number_of_space_time_coefficients();
                }
            }

            // set size for dQIdDof
            adQIdDof.set_size( tNumCoeff, 1, 0.0 );
            // END FIXME

            // get indices for properties, CM and SP
            uint tElastLinIsoIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );

            // compute dQIdDof for indirect dof dependencies
            for( uint iDof = 0; iDof < tRequestedDofTypes.size(); iDof++ )
            {
                // get treated dof type
                MSI::Dof_Type tDofType = tRequestedDofTypes( iDof );

                // get the set index for dof type
                sint tDofIndex = mSet->get_dof_index_for_type( tDofType, mtk::Master_Slave::MASTER );

                // get start and end indices for assembly
                uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
                uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

//                uint tDvIndex = mSet->get_requested_dv_index_for_type(     ---------------------    );
                uint tDvIndex = 0;

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( { tDofType } ) )
                {
                    // compute dQIdDof
                   mSet->get_residual()( tDvIndex )( { tStartRow, tEndRow }, { 0, 0 } )
                    += trans( mMasterCM( tElastLinIsoIndex )->dFluxdDOF( { tDofType } ) ) * mMasterCM( tElastLinIsoIndex )->strain( )
                     + trans( trans(mMasterCM( tElastLinIsoIndex )->flux()) * mMasterCM( tElastLinIsoIndex )->dStraindDOF( { tDofType } ) );
                }
            }
        }

//------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



