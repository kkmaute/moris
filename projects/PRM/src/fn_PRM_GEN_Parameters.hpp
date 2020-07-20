#ifndef MORIS_CL_PRM_GEN_PARAMETERS_HPP_
#define MORIS_CL_PRM_GEN_PARAMETERS_HPP_

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {

        /**
         * Creates a parameter list used for the construction of the geometry engine. One of these parameter lists is
         * always required as the only parameter list in the first cell (index 0).
         *
         * @return GEN parameter list
         */
        ParameterList create_gen_parameter_list()
        {
            ParameterList tParameterList;

            // Level-set parameters
            tParameterList.insert("isocontour_threshold", 0.0); // Level-set isocontour level
            tParameterList.insert("isocontour_error_factor", 0.0); // Error factor for determining interface nodes
            tParameterList.insert("level_set_file", ""); // Base file name (without extension) for saving level-set data

            // ADVs/IQIs
            tParameterList.insert("initial_advs", "");       // Initial ADVs, string converted into vector
            tParameterList.insert("advs_size", 0);           // Specify size and fill value for ADVs in addition to
            tParameterList.insert("initial_advs_fill", 0.0); // explicitly defined ADVs (above)
            tParameterList.insert("lower_bounds", "");       // Lower bounds on advs, string converted into vector
            tParameterList.insert("lower_bounds_fill", 0.0); // Fill value for lower bounds up to size of ADV vector
            tParameterList.insert("upper_bounds", "");       // Upper bounds on advs, string converted into vector
            tParameterList.insert("upper_bounds_fill", 0.0); // Fill value for upper bounds up to size of ADV vector
            tParameterList.insert("IQI_types", "");          // Requested IQI types for sensitivity analysis
            tParameterList.insert("PDV_types", "");          // Requested PDV types for sensitivity analysis

            // Phase table
            tParameterList.insert("phase_table", ""); // Construct phase table directly
            tParameterList.insert("phase_table_structure", "exp_base_2"); // Phase table structure (if not using phase_table)

            return tParameterList;
        }

        /**
         * Creates a parameter list that can create a geometry field. Any number of these can be added to the second
         * cell (index 1) of then parameter lists for GEN.
         *
         * @return Geometry parameter list
         */
        ParameterList create_geometry_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert("type", "");                      // Type (name) of geometry
            tParameterList.insert("geometry_variable_indices", ""); // Geometry variables to fill
            tParameterList.insert("adv_indices", "");               // ADVs used to fill in variables
            tParameterList.insert("constant_parameters", "");       // Remaining geometry parameters that are constant
            tParameterList.insert("number_of_refinements", 0);      // Number of refinement steps using HMR
            tParameterList.insert("refinement_function_index", -1); // User-defined refinement function (-1 = none)

            return tParameterList;
        }

        /**
         * Same as a geometry parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a geometry field and to evaluate sensitivities.
         *
         * @return User-defined geometry parameter list
         */
        ParameterList create_user_defined_geometry_parameter_list()
        {
            ParameterList tParameterList = create_geometry_parameter_list();

            tParameterList.set("type", "user_defined");             // User-defined geometry
            tParameterList.insert("field_function_name", "");       // Function name for evaluating the geometry field
            tParameterList.insert("sensitivity_function_name", ""); // Function name for evaluating the sensitivity of the field

            return tParameterList;
        }

        /**
         * Creates a parameter list that can create a property field. Any number of these can be added to the third cell
         * (index 2) of the parameter lists for GEN.
         *
         * @return GEN property parameter list
         */
        ParameterList create_gen_property_parameter_list()
        {
            ParameterList tParameterList;

            tParameterList.insert("type", "");                      // Type of property
            tParameterList.insert("name", "");                      // Name of property
            tParameterList.insert("property_variable_indices", ""); // Property variables to fill
            tParameterList.insert("adv_indices", "");               // ADVs used to fill in variables
            tParameterList.insert("constant_parameters", "");       // Remaining property parameters that are constant
            tParameterList.insert("dependencies", "");              // Names of other properties that this property depends on

            // Assignment to PDVs
            tParameterList.insert("pdv_type", "");
            tParameterList.insert("pdv_mesh_type", "interpolation");
            tParameterList.insert("pdv_mesh_set_names", "");
            tParameterList.insert("pdv_mesh_set_indices", "");

            return tParameterList;
        }

        /**
         * Same as a property parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a property field and to evaluate sensitivities.
         *
         * @return User-defined GEN property parameter list
         */
        ParameterList create_user_defined_property_parameter_list()
        {
            ParameterList tParameterList = create_gen_property_parameter_list();

            tParameterList.set("type", "user_defined"); // User-defined property
            tParameterList.insert("field_function_name", ""); // Function name for evaluating the property field
            tParameterList.insert("sensitivity_function_name", ""); // Function name for evaluating the sensitivity of the field

            return tParameterList;
        }

    }
}

#endif /* MORIS_CL_PRM_GEN_PARAMETERS_HPP_ */
