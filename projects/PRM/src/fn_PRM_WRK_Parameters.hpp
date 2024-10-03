/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_WRK_Parameters.hpp
 *
 */

#pragma once

#include "cl_Submodule_Parameter_Lists.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline Parameter_List
    create_wrk_parameter_list()
    {
        Parameter_List tParameterList( "Workflow" );

        tParameterList.insert( "adv_field", "" );
        tParameterList.insert( "dof_type", "" );
        tParameterList.insert( "reinitialization_frequency", 1 );
        tParameterList.insert( "output_mesh_file", "" );
        tParameterList.insert( "time_offset", 0.0 );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
