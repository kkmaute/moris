/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_STK_Parameters.hpp
 *
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_STK_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_STK_PARAMETERS_HPP_

#include "cl_Param_List.hpp"
namespace moris
{
    namespace prm
    {
        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        inline ParameterList
        create_stk_parameter_list()
        {
            ParameterList tParameterList;

            // decomposition and decomposition related parameters
            tParameterList.insert( "input_file", "" );
            tParameterList.insert( "periodic_workspace", false );
            tParameterList.insert( "periodic_side_set_pair", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif // PROJECTS_PRM_SRC_FN_PRM_STK_PARAMETERS_HPP_
