/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Normal_Vector.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Normal_Vector.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Normal_Vector::IQI_Normal_Vector()
        {
            mFEMIQIType = fem::IQI_Type::NORMAL_VECTOR;
        }

        //------------------------------------------------------------------------------

        void
        IQI_Normal_Vector::compute_QI( Matrix< DDRMat >& aQI )
        {
            // get the constitutive model for computing the fluxes/tractions

            MORIS_ASSERT( mNormal.numel() > 0,
                    "IQI_Normal_Vector::compute_QI() - "
                    "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                    "and averaged such that there is a well-defined normal." );

            aQI = { { mNormal( mIQITypeIndex ) } };
        }

        //------------------------------------------------------------------------------

        void
        IQI_Normal_Vector::compute_QI( real aWStar )
        {
            MORIS_ERROR( false, "Normal Vector IQI not implemented for integrated evalutaion." );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Normal_Vector::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        void
        IQI_Normal_Vector::compute_dQIdu(
                moris::Cell< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&             adQIdu )
        {
            MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
        }

        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------

    } /* end namespace fem */
} /* end namespace moris */
