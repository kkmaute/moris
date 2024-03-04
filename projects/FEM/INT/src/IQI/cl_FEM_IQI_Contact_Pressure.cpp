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

        IQI_Contact_Pressure::IQI_Contact_Pressure()
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

            Matrix< DDRMat >       tNormalLeader = tLeaderGeometry->get_normal_current( tLeaderDofs );
            const Matrix< DDRMat > tTraction     = tCMLeader->traction( tNormalLeader );    // ( 2 x 1 )
            const real             tPressure     = dot( tTraction, tNormalLeader );

            // TODO @ff: remove this
            // [[maybe_unused]] uint const tIteration = gLogger.get_iteration( "NonLinearAlgorithm", "Newton", "Solve" );
            // std::cout << "CP:"
            //           << tIteration << ","
            //           << tPressure << ","
            //           << tLeaderGeometry->valx() << ","
            //           << tLeaderGeometry->valx_current( tLeaderDofs ) << ","
            //           << tNormalLeader << "\n";

            // pick the component otherwise (0,1,2)
            aQI = { tPressure };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Contact_Pressure::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            MORIS_ASSERT( mNormal.numel() > 0,
                    "IQI_Contact_Pressure::compute_QI() - "
                    "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                    "and averaged such that there is a well-defined normal." );

            // get the constitutive model for computing the fluxes/tractions
            const std::shared_ptr< Constitutive_Model >& tCMLeader       = mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );
            Field_Interpolator*                          tLeaderDofs     = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Geometry_Interpolator*                       tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

            Matrix< DDRMat >       tNormalLeader = tLeaderGeometry->get_normal_current( tLeaderDofs );
            const Matrix< DDRMat > tTraction     = tCMLeader->traction( tNormalLeader );    // ( 2 x 1 )
            const real             tPressure     = dot( tTraction, tNormalLeader );

            // initialize the output matrix
            Matrix< DDRMat > tMat( { { tPressure } } );

            // add the contribution
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
