/*
 * cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.cpp
 *
 *  Created on: Sep 27, 2020
 *      Author: noel
 */
#include "cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Thermal_Energy_Diffusive_Flux::IQI_Thermal_Energy_Diffusive_Flux()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IQI_Constitutive_Type::DIFFUSION );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Diffusive_Flux::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // evaluate the QI
            // note: because diffusive flux is defined as "conductivity tensor * gradx(Temp)"
            //       the numerical diffusive flux needs to be multiplied by (-1) to obtain
            //       the physical flux
            aQI = -1.0 * dot(tCMDiffusion->flux(),mNormal);
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Diffusive_Flux::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // evaluate the QI
            // note: because diffusive flux is defined as "conductivity tensor * gradx(Temp)"
             //       the numerical diffusive flux needs to be multiplied by (-1) to obtain
             //       the physical flux
            mSet->get_QI()( tQIIndex ) += -1.0 * aWStar * dot(tCMDiffusion->flux(),mNormal);
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Diffusive_Flux::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

                // if density depends on dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    -1.0 * trans( tCMDiffusion->dFluxdDOF( tDofType ) ) * mNormal );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Thermal_Energy_Diffusive_Flux::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

            // if density depends on dof type
            if ( tCMDiffusion->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu += -1.0 * trans( tCMDiffusion->dFluxdDOF( aDofType ) ) * mNormal;
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

