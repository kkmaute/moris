/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_MIG_Parameters.hpp
 *
 */

#pragma once

#include "cl_Submodule_Parameter_Lists.hpp"

namespace moris::prm
{
    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline Parameter_List
    create_mig_parameter_list()
    {
        Parameter_List tParameterList( "MIG" );

        tParameterList.insert( "periodic_side_set_pair", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
