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

        void IQI_Property::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IQI_Property::set_property - can only be master." );

            // FIXME check that property type makes sense?

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------

        IQI_Property::IQI_Property()
        {
            // set IQI type
            mIQIType = vis::Output_Type::PROPERTY;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Property" ] = IQI_Property_Type::PROPERTY;
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

        void IQI_Property::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            MORIS_ERROR( false, "IQI_Property::compute_dQIdu - Not implemented." );
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */



