/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Contact_Pressure.cpp
 *
 */

#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Contact_Pressure.hpp"
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include <iostream>
#include <memory>

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Contact_Pressure::IQI_Contact_Pressure( bool aInReferenceConfiguration )
                : mInReferenceConfiguration( aInReferenceConfiguration )
        {
            mFEMIQIType = fem::IQI_Type::JUMP_TRACTION;

            // set size for the constitutive model pointer cell
            mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );
            mFollowerCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Contact_Pressure::compute_QI( Matrix< DDRMat >& aQI )
        {
            MORIS_ASSERT( mNormal.numel() > 0,
                    "IQI_Contact_Pressure::compute_QI() - "
                    "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                    "and averaged such that there is a well-defined normal." );

            // get the constitutive model for computing the fluxes/tractions
            const std::shared_ptr< Constitutive_Model >& tCMLeader       = mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
            Field_Interpolator*                          tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Geometry_Interpolator*                       tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

            Matrix< DDRMat > tNormal;
            Matrix< DDRMat > tTraction;

            if ( mInReferenceConfiguration )
            {
                // get the pressure in the reference configuration (small deformation)!
                tNormal   = tLeaderGeometry->get_normal();
                tTraction = tCMLeader->traction( tNormal );
            }
            else
            {
                // get the pressure in the current configuration (large deformation)!
                tNormal   = tLeaderGeometry->get_normal_current( tLeaderDofs );
                tTraction = tCMLeader->traction( tNormal, CM_Function_Type::PK1 );
            }

            const real tPressure = dot( tTraction, tNormal );

            // pick the component otherwise (0,1,2)
            aQI = { tPressure };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Contact_Pressure::compute_QI( real aWStar )
        {
            Matrix< DDRMat > tMat( 1, 1 );
            this->compute_QI( tMat );

            // add the contribution
            sint tQIIndex = mSet->get_QI_assembly_index( mName );
            mSet->get_QI()( tQIIndex ) += aWStar * tMat;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Contact_Pressure::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Contact_Pressure::compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
