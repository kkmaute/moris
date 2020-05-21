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
        ParameterList create_gen_parameter_list();

        /**
         * Creates a parameter list that can create a geometry field. Any number of these can be added to the second
         * cell (index 1) of then parameter lists for GEN.
         *
         * @return Geometry parameter list
         */
        ParameterList create_geometry_parameter_list();

        /**
         * Same as a geometry parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a geometry field and to evaluate sensitivities.
         *
         * @return User-defined geometry parameter list
         */
        ParameterList create_user_defined_geometry_parameter_list();

        /**
         * Creates a parameter list that can create a property field. Any number of these can be added to the third cell
         * (index 2) of the parameter lists for GEN.
         *
         * @return GEN property parameter list
         */
        ParameterList create_gen_property_parameter_list();

        /**
         * Same as a property parameter list, but forces the user to specify two more parameters which name the
         * user-defined functions used to evaluate a property field and to evaluate sensitivities.
         *
         * @return User-defined GEN property parameter list
         */
        ParameterList create_user_defined_property_parameter_list();


    } // end prm namespace
} // end moris namespace

#endif /* MORIS_CL_PRM_GEN_PARAMETERS_HPP_ */
