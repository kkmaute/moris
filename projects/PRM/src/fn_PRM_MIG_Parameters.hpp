/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_MIG_Parameters.hpp
 *
 */

#ifndef SRC_fn_PRM_MIG_Parameters
#define SRC_fn_PRM_MIG_Parameters

#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace prm
    {
        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        inline Parameter_List
        create_mig_parameter_list()
        {
            Parameter_List tParameterList;

            tParameterList.insert( "periodic_side_set_pair", "" );

            return tParameterList;
        }

        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif    // SRC_fn_PRM_MIG_Parameters
