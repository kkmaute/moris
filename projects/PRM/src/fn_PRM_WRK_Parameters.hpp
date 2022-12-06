/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_WRK_Parameters.hpp
 *
 */

#ifndef SRC_fn_PRM_WRK_Parameters
#define SRC_fn_PRM_WRK_Parameters

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        inline ParameterList
        create_wrk_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert( "adv_field", "" );
            tParameterList.insert( "dof_type", "" );
            tParameterList.insert( "reinitialization_frequency", 1 );
            tParameterList.insert( "output_mesh_file", "" );
            tParameterList.insert( "time_offset", 0.0 );

            return tParameterList;
        }
        
        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif
