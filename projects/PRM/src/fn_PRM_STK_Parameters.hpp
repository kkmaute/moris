/*
 * fn_PRM_STK_Parameters.hpp
 *
 *  Created on: March 10, 2020
 *      Author: doble
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_STK_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_STK_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
//#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {

//------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline
    ParameterList create_stk_parameter_list()
    {
        ParameterList tParameterList;

        // decomposition and decomposition related parameters
        tParameterList.insert( "input_file", "" );
        tParameterList.insert( "periodic_workspace", false);
        tParameterList.insert( "periodic_side_set_pair", "");

        return tParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /*  */
