/*
 * cl_FEM_IQI_Property.cpp
 *
 *  Created on: Feb 2, 2020
 *      Author: noel
 */
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Property.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Property::IQI_Property()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Property" ] = static_cast< uint >( IQI_Property_Type::PROPERTY );
        }
        
        //------------------------------------------------------------------------------

        void IQI_Property::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get property index
            uint tPropertyIndex = static_cast< uint >( IQI_Property_Type::PROPERTY );

            // evaluate the QI
            aQI = mMasterProp( tPropertyIndex )->val();
        }

        //------------------------------------------------------------------------------

        void IQI_Property::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the property
            std::shared_ptr< Property > & tProperty =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PROPERTY ) );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * ( tProperty->val() );
        }

        //------------------------------------------------------------------------------

        void IQI_Property::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the property
            std::shared_ptr< Property > & tProperty =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PROPERTY ) );

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

                // if property depends on dof type
                if ( tProperty->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex },
                            { 0, 0 } ) += aWStar * trans( tProperty->dPropdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Property::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get the property
            std::shared_ptr< Property > & tProperty =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::PROPERTY ) );

            // if property depends on dof type
            if ( tProperty->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = trans( tProperty->dPropdDOF( aDofType ) );
            }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



