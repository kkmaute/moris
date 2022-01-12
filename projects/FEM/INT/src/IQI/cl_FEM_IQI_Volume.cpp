/*
 * cl_FEM_IQI_Volume.cpp
 *
 *  Created on: Nov 20, 2019
 *      Author: noel
 */

#include "cl_FEM_IQI_Volume.hpp"
#include "cl_FEM_Set.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Volume::IQI_Volume()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Density" ] = static_cast< uint >( IQI_Property_Type::DENSITY );
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                aQI = tPropDensity->val();
            }
            else
            {
                // set density to 1
                aQI = {{ 1.0 }};
            }
        }
        
        //------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                mSet->get_QI()( tQIIndex ) += aWStar * ( tPropDensity->val() );
            }
            else
            {
                // set density to 1
                mSet->get_QI()( tQIIndex ) += aWStar;
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get density property
            const std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

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

                // Dof dependency
                if ( tPropDensity != nullptr && tPropDensity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * ( tPropDensity->dPropdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get density property
            std::shared_ptr< Property > & tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // Dof dependency
            if ( tPropDensity != nullptr && tPropDensity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = tPropDensity->dPropdDOF( aDofType );
            }
        }

        //------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */
