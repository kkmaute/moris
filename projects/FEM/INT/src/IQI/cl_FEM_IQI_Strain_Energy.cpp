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
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST );
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // evaluate the QI
            aQI = 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain();
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain() );
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

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

                // if elasticity CM depends on dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * 0.5 * (
                                    trans( tCMElasticity->dFluxdDOF( tDofType ) )   * tCMElasticity->strain() +
                                    trans( tCMElasticity->dStraindDOF( tDofType ) ) * tCMElasticity->flux() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // if elasticity CM depends on dof type
            if ( tCMElasticity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = 0.5 * (
                        trans( tCMElasticity->dFluxdDOF( aDofType ) )   * tCMElasticity->strain() +
                        trans( tCMElasticity->dStraindDOF( aDofType ) ) * tCMElasticity->flux() );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

