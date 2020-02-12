/*
 * cl_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#ifndef PROJECTS_FEM_INT_SRC_CL_FEM_PARAMETERS_HPP_
#define PROJECTS_FEM_INT_SRC_CL_FEM_PARAMETERS_HPP_

#include <string>
#include <cstdio>

#include "assert.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_XML_Parser.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        // creates a property parameter list with default inputs
        ParameterList create_property_parameter_list();

//------------------------------------------------------------------------------
        // creates a constitutive model parameter list with default inputs
        ParameterList create_constitutive_model_parameter_list();

//------------------------------------------------------------------------------
        // creates a stabilization parameter parameter list with default inputs
        ParameterList create_stabilization_parameter_parameter_list();

//------------------------------------------------------------------------------
        // creates an IWG parameter list with default inputs
        ParameterList create_IWG_parameter_list();

//------------------------------------------------------------------------------
        // creates an IQI parameter list with default inputs
        ParameterList create_IQI_parameter_list();

//------------------------------------------------------------------------------

    }/* end_namespace_fem */
}/* end_namespace_moris */

#endif /* PROJECTS_FEM_INT_SRC_CL_FEM_PARAMETERS_HPP_ */
