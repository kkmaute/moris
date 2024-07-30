/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Gap.cpp
 *
 */

#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Gap.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Gap::IQI_Gap()
        {
            mFEMIQIType = fem::IQI_Type::GAP;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Gap::compute_QI( Matrix< DDRMat >& aQI )
        {
            Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            Geometry_Interpolator* tGILeader   = mLeaderFIManager->get_IG_geometry_interpolator();
            Geometry_Interpolator* tGIFollower = mFollowerFIManager->get_IG_geometry_interpolator();

            Matrix< DDRMat > tCurrentNormal = tGILeader->get_normal_current( tFILeader );

            real tGap = dot(
                    tGIFollower->valx_current( tFIFollower ) - tGILeader->valx_current( tFILeader ),
                    tCurrentNormal );

            aQI = { { tGap } };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Gap::compute_QI( real aWStar )
        {
            Field_Interpolator* tFILeader   = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            Field_Interpolator* tFIFollower = mFollowerFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            Geometry_Interpolator* tGILeader   = mLeaderFIManager->get_IG_geometry_interpolator();
            Geometry_Interpolator* tGIFollower = mFollowerFIManager->get_IG_geometry_interpolator();

            Matrix< DDRMat > tCurrentNormal = tGILeader->get_normal_current( tFILeader );

            real tGap = dot(
                    tGIFollower->valx_current( tFIFollower ) - tGILeader->valx_current( tFILeader ),
                    tCurrentNormal );

            sint tQIIndex = mSet->get_QI_assembly_index( mName );
            mSet->get_QI()( tQIIndex ) += aWStar * tGap;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Gap::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Gap::compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
