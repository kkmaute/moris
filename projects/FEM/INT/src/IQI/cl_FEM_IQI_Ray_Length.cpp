/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Ray_Length.cpp
 *
 */

#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Ray_Length.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Ray_Length::IQI_Ray_Length()
        {
            mFEMIQIType = fem::IQI_Type::RAY_LENGTH;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Ray_Length::compute_QI( Matrix< DDRMat >& aQI )
        {
            Geometry_Interpolator* tGILeader   = mLeaderFIManager->get_IG_geometry_interpolator();
            Geometry_Interpolator* tGIFollower = mFollowerFIManager->get_IG_geometry_interpolator();

            real tGap = norm( tGIFollower->valx() - tGILeader->valx() );

            aQI = { { tGap } };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Ray_Length::compute_QI( real aWStar )
        {
            Geometry_Interpolator* tGILeader   = mLeaderFIManager->get_IG_geometry_interpolator();
            Geometry_Interpolator* tGIFollower = mFollowerFIManager->get_IG_geometry_interpolator();

            real tGap     = norm( tGIFollower->valx() - tGILeader->valx() );
            sint tQIIndex = mSet->get_QI_assembly_index( mName );
            mSet->get_QI()( tQIIndex ) += aWStar * tGap;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Ray_Length::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Ray_Length::compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
