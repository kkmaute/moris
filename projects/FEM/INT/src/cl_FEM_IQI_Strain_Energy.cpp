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

            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Elast" ] = IQI_Constitutive_Type::ELAST;
        }

        //------------------------------------------------------------------------------
        void IQI_Strain_Energy::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(),
                    "IQI_Strain_Energy::set_constitutive_model - Unknown aConstitutiveString." );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Strain_Energy::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
        void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get indices for properties, CM and SP
            uint tElastIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST );

            // evaluate the QI
            aQI = 0.5 * trans( mMasterCM( tElastIndex )->flux() ) * mMasterCM( tElastIndex )->strain();
        }

        //------------------------------------------------------------------------------
        void IQI_Strain_Energy::compute_QI( moris::real aWStar )
        {
            // get indices for properties, CM and SP
            uint tElastIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST );

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += 0.5 * aWStar * trans( mMasterCM( tElastIndex )->flux() ) * mMasterCM( tElastIndex )->strain();

        }

        //------------------------------------------------------------------------------
        void IQI_Strain_Energy::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_requested_dof_types();

            // get indices for properties, CM and SP
            uint tElastIndex = static_cast< uint >( IQI_Constitutive_Type::ELAST );

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

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( tElastIndex )->check_dof_dependency( { tDofType } ) )
                {
                    // compute dQIdDof
                   mSet->get_residual()( tQIIndex )( { tStartRow, tEndRow }, { 0, 0 } )
                    += 0.5 * aWStar * ( trans( mMasterCM( tElastIndex )->dFluxdDOF( { tDofType } ) ) * mMasterCM( tElastIndex )->strain( )
                    +  trans( trans( mMasterCM( tElastIndex )->flux() )             * mMasterCM( tElastIndex )->dStraindDOF( { tDofType } ) ) );

                }
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



