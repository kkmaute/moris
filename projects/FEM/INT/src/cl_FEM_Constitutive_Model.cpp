/*
 * cl_FEM_Constitutive_Model.cpp
 *
 *  Created on: Dec 10, 2019
 *      Author: noel
 */

#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void Constitutive_Model::set_field_interpolator_manager( Field_Interpolator_Manager * aFieldInterpolatorManager )
        {
            // set the field interpolator manager for the constitutive model
            mFIManager = aFieldInterpolatorManager;

            // FIXME
            // get the list of dof types for the CM
            moris::Cell< moris::Cell< MSI::Dof_Type > > tCMDofTypes
            = this->get_global_dof_type_list();

            // get the number of dof type for the CM
            uint tNumDofTypes = tCMDofTypes.size();

            // set the size of the field interpolators list for the CM
            mDofFI.resize( tNumDofTypes, nullptr );

            // loop over the dof types
            for( uint iDof = 0; iDof < tNumDofTypes; iDof++ )
            {
                // fill the field interpolators list for the CM
                mDofFI( iDof ) = mFIManager->get_field_interpolators_for_type( tCMDofTypes( iDof )( 0 ) );
            }
            // END FIXME

            // loop over the underlying properties
            for( std::shared_ptr< Property > tProp : this->get_properties() )
            {
                if (tProp != nullptr )
                {
                    // set the field interpolator manager for the property
                    tProp->set_field_interpolator_manager( mFIManager );
                }
            }
        }

    }/* end_fem_namespace */
}/* end_moris_namespace */
