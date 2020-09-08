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

        void IQI_Volume::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // get density property
            std::shared_ptr< Property > tPropDensity =
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


