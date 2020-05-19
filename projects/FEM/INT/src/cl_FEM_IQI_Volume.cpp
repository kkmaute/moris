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
            // set the size for the QI
            aQI.set_size( 1, 1 );

            // evaluate the QI
            aQI = {{ 1.0 }};
        }

//------------------------------------------------------------------------------

        void IQI_Volume::compute_QI( moris::real aWStar )
        {
            // get property index
            uint tPropertyIndex = static_cast< uint >( IQI_Property_Type::DENSITY );

            std::shared_ptr< Property > tPropDensity = mMasterProp( tPropertyIndex );

            Matrix< DDRMat > tPropVal;

            if ( tPropDensity != nullptr )
            {
                // evaluate the QI
                tPropVal = mMasterProp( tPropertyIndex )->val();
            }
            else
            {
                tPropVal = {{ 1.0 }};
            }

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            //print( mSet->get_QI(), "mSet->get_QI()");
            // evaluate the QI
            mSet->get_QI()( tQIIndex ).matrix_data() += aWStar * tPropVal;
        }

//------------------------------------------------------------------------------

        void IQI_Volume::set_property( std::shared_ptr< Property > aProperty,
                                       std::string                 aPropertyString,
                                       mtk::Master_Slave           aIsMaster )
        {
            // can only be master
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                         "IQI::set_property - can only be master." );

            // FIXME check that property type makes sense?

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

//------------------------------------------------------------------------------
    }/* end namespace fem */
}/* end namespace moris */


