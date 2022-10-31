/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MSI_Design_Variable_Interface.cpp
 *
 */

#include "cl_MSI_Design_Variable_Interface.hpp"
#include "cl_MSI_Equation_Model.hpp"
// SOL/src
#include "cl_SOL_Dist_Vector.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace MSI
    {

        //------------------------------------------------------------------------------

        void
        Design_Variable_Interface::set_requested_IQIs(
                const moris::Cell< std::string >& aRequestedIQIs )
        {
            MORIS_ASSERT( mModel != nullptr,
                    "Design_Variable_Interface::set_requested_IQIs - mModel has not been set." );

            mModel->set_requested_IQI_names( aRequestedIQIs );
        }

        //------------------------------------------------------------------------------

        sol::Dist_Vector*
        Design_Variable_Interface::get_dQIdp()
        {
            MORIS_ASSERT( mModel != nullptr,
                    "Design_Variable_Interface::get_dQIdp - mModel has not been set." );

            return mModel->get_dQIdp();
        }

        //------------------------------------------------------------------------------

    }    // namespace MSI
}    // namespace moris
