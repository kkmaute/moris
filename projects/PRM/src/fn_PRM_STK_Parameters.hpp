/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_STK_Parameters.hpp
 *
 */

#pragma once

#include "cl_Parameter_List.hpp"

namespace moris::prm
{
    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline Parameter_List
    create_stk_parameter_list()
    {
        Parameter_List tParameterList( "STK" );

        // decomposition and decomposition related parameters
        tParameterList.insert( "input_file", "" );
        tParameterList.insert( "periodic_workspace", false );
        tParameterList.insert( "periodic_side_set_pair", "" );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
