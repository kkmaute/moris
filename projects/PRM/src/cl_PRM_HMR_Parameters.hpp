/*
 * cl_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_HMR_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_HMR_PARAMETERS_HPP_

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
    ParameterList create_hmr_parameter_list()
    {
        ParameterList tParameterList;

        tParameterList.insert( "number_of_elements_per_dimension", std::string( "2, 2" ) );
        tParameterList.insert( "domain_dimensions", std::string( "1, 1" ) );
        tParameterList.insert( "domain_offset", std::string( "0, 0 ") );
        tParameterList.insert( "domain_sidesets", std::string( "" ) );
        tParameterList.insert( "lagrange_output_meshes", std::string( "" ) );

        tParameterList.insert( "lagrange_input_meshes", std::string( "" ) );

        tParameterList.insert( "refinement_buffer", 0 );
        tParameterList.insert( "staircase_buffer", 0 );

        tParameterList.insert( "lagrange_orders", std::string( "1" ) );
        tParameterList.insert( "lagrange_pattern", std::string( "0" ) );

        tParameterList.insert( "bspline_orders", std::string( "1" ) );
        tParameterList.insert( "bspline_pattern", std::string( "0" ) );

        tParameterList.insert( "union_pattern", 6 );
        tParameterList.insert( "working_pattern", 7 );

        tParameterList.insert( "lagrange_to_bspline", std::string( "0" ) );

        tParameterList.insert( "severity_level", 1 );
        tParameterList.insert( "truncate_bsplines", 1 );

        tParameterList.insert( "use_multigrid", 0 );
        tParameterList.insert( "use_refinement_interrelation", 0 );
        tParameterList.insert( "renumber_lagrange_nodes", 0 );
        tParameterList.insert( "use_number_aura", 0 );

        tParameterList.insert( "initial_refinement", 0 );
        tParameterList.insert( "additional_lagrange_refinement", 0 );

        tParameterList.insert( "max_refinement_level", -1 );

        tParameterList.insert( "adaptive_refinement_level", 0 );

        return tParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_ */
