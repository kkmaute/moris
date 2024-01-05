/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Debug.cpp
 *
 */

#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Debug.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Debug::IQI_Debug()
        {
            mFEMIQIType = fem::IQI_Type::DEBUG;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Debug::compute_QI( Matrix< DDRMat >& aQI )
        {
            moris_index const iGP        = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->get_current_ig_point_index();
            Vector< real >    tDistances = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->get_integration_point_distances();

            MORIS_ASSERT( !tDistances.empty(), "IQI_Debug::compute_QI( real aWStar ): Integration point distances are empty. Are you using it on a Nonconformal Side Set?" );

            std::cout << "QI (aQI): " << mName << " " << iGP << " "
                      << " " << tDistances( iGP ) << "\n";
            aQI = { { tDistances( iGP ) } };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Debug::compute_QI( real aWStar )
        {
            moris_index const iGP        = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->get_current_ig_point_index();
            Vector< real >    tDistances = mSet->get_field_interpolator_manager()->get_IG_geometry_interpolator()->get_integration_point_distances();

            MORIS_ASSERT( !tDistances.empty(), "IQI_Debug::compute_QI( real aWStar ): Integration point distances are empty. Are you using it on a Nonconformal Side Set?" );

            std::cout << "QI (aWStar): " << mName << " " << iGP << " " << aWStar << " " << tDistances( iGP ) << "\n";
        }

        //------------------------------------------------------------------------------

        void
        IQI_Debug::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Debug::compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
