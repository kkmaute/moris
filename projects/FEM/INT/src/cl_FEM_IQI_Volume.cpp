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

        void IQI_Volume::compute_QI( Matrix< DDRMat > & aQI )
        {
            // get density property
            std::shared_ptr< Property > tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // init density value
            Matrix< DDRMat > tPropVal;

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                tPropVal = tPropDensity->val();
            }
            else
            {
                // set density to 1
                tPropVal = {{ 1.0 }};
            }

            // evaluate the QI
            aQI = tPropVal;
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( moris::real aWStar )
        {
            // get density property
            std::shared_ptr< Property > tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // init density value
            Matrix< DDRMat > tPropVal;

            // if density property
            if ( tPropDensity != nullptr )
            {
                // evaluate the density
                tPropVal = tPropDensity->val();
            }
            else
            {
                // set density to 1
                tPropVal = {{ 1.0 }};
            }

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += aWStar * tPropVal;
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Volume::set_property - can only be master." );

            // check that aPropertyString makes sense
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(),
                    "IQI_Volume::set_property - Unknown aPropertyString." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        void IQI_Volume::compute_dQIdu( real aWStar )
        {
            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get density property
            std::shared_ptr< Property > tPropDensity =
                    mMasterProp( static_cast< uint >( IQI_Property_Type::DENSITY ) );

            // get the requested dof types
            moris::Cell < enum MSI::Dof_Type > tRequestedDofTypes =
                    this->get_requested_dof_types();

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
                if ( tPropDensity != nullptr && tPropDensity->check_dof_dependency( { tDofType } ) )
                {
                    // compute dQIdDof
                    mSet->get_residual()( tQIIndex )( { tStartRow, tEndRow }, { 0, 0 } ) +=
                            aWStar * tPropDensity->dPropdDOF( { tDofType } );
                }
            }
        }

        //------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


