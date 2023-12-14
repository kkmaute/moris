/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Debug.cpp
 *
 */

#include "cl_FEM_IWG_Debug.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Debug::IWG_Debug() = default;

        //------------------------------------------------------------------------------

        void IWG_Debug::compute_residual( real aWStar )
        {
        }

        //------------------------------------------------------------------------------
        void IWG_Debug::compute_jacobian( real aWStar )
        {
        }

        //------------------------------------------------------------------------------

        void IWG_Debug::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Debug::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void IWG_Debug::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Debug::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
