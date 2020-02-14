/*
 * cl_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: noel
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_

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
        /*
         * creates a property parameter list with default inputs
         * @param [ out ] ParameterList a property parameter list
         */
        ParameterList create_property_parameter_list();

//------------------------------------------------------------------------------
        /*
         * creates a constitutive model parameter list with default inputs
         * @param [ out ] ParameterList a CM parameter list
         */
        ParameterList create_constitutive_model_parameter_list();

//------------------------------------------------------------------------------
        /*
         * creates a stabilization parameter parameter list with default inputs
         * @param [ out ] ParameterList a SP parameter list
         */
        ParameterList create_stabilization_parameter_parameter_list();

//------------------------------------------------------------------------------
        /*
         * creates an IWG parameter list with default inputs
         * @param [ out ] ParameterList a IWG parameter list
         */
        ParameterList create_IWG_parameter_list();

//------------------------------------------------------------------------------
        /*
         * creates an IQI parameter list with default inputs
         * @param [ out ] ParameterList a IQI parameter list
         */
        ParameterList create_IQI_parameter_list();

//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_FEM_PARAMETERS_HPP_ */
