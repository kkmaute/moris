/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_SP_Stab_Penalty_Contact.cpp
 *
 */

#include "cl_FEM_SP_Stab_Penalty_Contact.hpp"   //FEM/INT/src
#include "cl_FEM_Cluster.hpp"              //FEM/INT/src

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        SP_Stab_Penalty_Contact::SP_Stab_Penalty_Contact()
        {
            // set size for the property pointer cells
            mLeaderProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
            mFollowerProp.resize( static_cast< uint >( SP_Property_Type::MAX_ENUM ), nullptr );
        }

        //------------------------------------------------------------------------------
        void SP_Stab_Penalty_Contact::eval_SP()
        {
            // compute stabilization parameter value
            mPPVal = mParameters( 0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

