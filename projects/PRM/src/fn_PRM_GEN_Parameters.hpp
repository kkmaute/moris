/*
 * cl_PRM_GEN_Parameters.hpp
 *
 *  Created on: Feb 18, 2020
 *      Author: sonne
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_

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

    //------------------------------------------------------------------------------------------------------------------

    ParameterList create_gen_parameter_list()
    {
        ParameterList tParameterList;

        tParameterList.insert("spatial_dimensions", 3); // Number of spatial dimensions
        tParameterList.insert("threshold_value", 0.0);
        tParameterList.insert("perturbation_value", 1E-6);
        tParameterList.insert("HMR_refinements", 0); // Number of HMR refinements to be performed
        tParameterList.insert("initial_advs", ""); // Initial advs, string converted into vector
        tParameterList.insert("lower_bounds", ""); // Lower bounds on advs, string converted into vector
        tParameterList.insert("upper_bounds", ""); // Upper bounds on advs, string converted into vector
        tParameterList.insert("phase_table", ""); // Construct phase table directly
        tParameterList.insert("phase_table_structure", "exp_base_2"); // Phase table structure (if not using phase_table)
        tParameterList.insert("library", ""); // File for loading user-defined functions

        return tParameterList;
    }

    //------------------------------------------------------------------------------------------------------------------

    ParameterList create_geometry_parameter_list()
    {
        ParameterList tParameterList;

        tParameterList.insert("type", ""); // Type (name) of geometry
        tParameterList.insert("geometry_variable_indices", ""); // String of uints converted to a vector;
                                                                // geometry variables to fill
        tParameterList.insert("adv_indices", ""); // String of uints converted to a vector;
                                                  // advs used to fill in variables
        tParameterList.insert("constant_parameters", ""); // String of reals converted to a vector;
                                                          // remaining geometry parameters that are constant

        return tParameterList;
    }

    //------------------------------------------------------------------------------------------------------------------

    ParameterList create_user_defined_geometry_parameter_list()
    {
        ParameterList tParameterList = create_geometry_parameter_list();

        tParameterList.set("type", "user_defined"); // User-defined geometry
        tParameterList.insert("field_function_name", ""); // Function name for evaluating the geometry field
        tParameterList.insert("sensitivity_function_name", ""); // Function name for evaluating the sensitivity of the field

        return tParameterList;
    }

    } // end prm namespace
} // end moris namespace



#endif /* PROJECTS_PRM_SRC_CL_PRM_GEN_PARAMETERS_HPP_ */
