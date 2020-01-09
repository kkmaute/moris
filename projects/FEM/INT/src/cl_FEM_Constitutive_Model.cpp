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

            // FIXME
            // get the list of dv types for the CM
            moris::Cell< moris::Cell< MSI::Dv_Type > > tCMDvTypes
            = this->get_global_dv_type_list();

            // get the number of dv type for the CM
            uint tNumDvTypes = tCMDvTypes.size();

            // set the size of the field interpolators list for the CM
            mDvFI.resize( tNumDvTypes, nullptr );

            // loop over the dof types
            for( uint iDv = 0; iDv < tNumDvTypes; iDv++ )
            {
                // fill the field interpolator list for the CM
                mDvFI( iDv ) = mFIManager->get_field_interpolators_for_type( tCMDvTypes( iDv )( 0 ) );
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

//------------------------------------------------------------------------------
        void Constitutive_Model::get_non_unique_dof_and_dv_types( moris::Cell< MSI::Dof_Type > & aDofTypes,
                                                                  moris::Cell< MSI::Dv_Type >  & aDvTypes )
        {
            // init dof counter
            uint tDofCounter = 0;
            uint tDvCounter  = 0;

            // loop over direct dof dependencies
            for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // update counter
                tDofCounter += mDofTypes( iDof ).size();
            }

            // loop over direct dv dependencies
            for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // update counter
                tDvCounter += mDvTypes( iDv ).size();
            }

            // loop over properties
            for ( std::shared_ptr< Property > tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // update counter
                    tDofCounter += tActiveDofTypes.size();
                    tDvCounter  += tActiveDvTypes.size();
                }
            }

            // reserve memory for the non unique dof and dv types
            aDofTypes.reserve( tDofCounter );
            aDvTypes.reserve( tDvCounter );

            // loop over direct dof dependencies
            for ( uint iDof = 0; iDof < mDofTypes.size(); iDof++ )
            {
                // populate the dof type list
                aDofTypes.append( mDofTypes( iDof ) );
            }

            // loop over direct dv dependencies
            for ( uint iDv = 0; iDv < mDvTypes.size(); iDv++ )
            {
                // populate the dv type list
                aDvTypes.append( mDvTypes( iDv ) );
            }

            // loop over the properties
            for ( std::shared_ptr< Property > tProperty : mProperties )
            {
                if ( tProperty != nullptr )
                {
                    // get property dof and dv type list
                    moris::Cell< MSI::Dof_Type > tActiveDofTypes;
                    moris::Cell< MSI::Dv_Type >  tActiveDvTypes;
                    tProperty->get_non_unique_dof_and_dv_types( tActiveDofTypes,
                                                                tActiveDvTypes );

                    // populate the dof and dv type lists
                    aDofTypes.append( tActiveDofTypes );
                    aDvTypes.append( tActiveDvTypes );
                }
            }
        }

    }/* end_fem_namespace */
}/* end_moris_namespace */
