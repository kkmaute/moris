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

#include <string>
#include <cstdio>

#include "assert.hpp"
//#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Param_List.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline ParameterList
    create_mig_parameter_list()
    {
        ParameterList tParameterList;

        tParameterList.insert( "periodic_side_set_pair", "" );

        return tParameterList;
    }
    //------------------------------------------------------------------------------

}// namespace moris::prm

#endif
