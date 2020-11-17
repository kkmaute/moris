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


