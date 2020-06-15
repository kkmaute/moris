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

        //! number of elements per direction in overall mesh, without aura
        //! 2D or 3D is determined by length of this vector
        tParameterList.insert( "number_of_elements_per_dimension", std::string( "2, 2" ) );
        //! width, height and depth of domain (without aura)
        tParameterList.insert( "domain_dimensions", std::string( "1, 1" ) );
        //! offset from the origin
        tParameterList.insert( "domain_offset", std::string( "0, 0 ") );
        //! sidesets which should be built
        tParameterList.insert( "domain_sidesets", std::string( "" ) );

        //! Lagrange Meshes that are used for the output meshes
        tParameterList.insert( "lagrange_output_meshes", std::string( "" ) );
        //! Lagrange Meshes that are used as input meshes
        tParameterList.insert( "lagrange_input_meshes", std::string( "" ) );

        //! size of refinement buffer
        tParameterList.insert( "refinement_buffer", 0 );
        //! size of staircase buffer
        tParameterList.insert( "staircase_buffer", 0 );

        //! Lagrange orders
        tParameterList.insert( "lagrange_orders", std::string( "1" ) );
        //! Lagrange pattern
        tParameterList.insert( "lagrange_pattern", std::string( "0" ) );

        //! B-Spline orders
        tParameterList.insert( "bspline_orders", std::string( "1" ) );
        //! B-Spline orders
        tParameterList.insert( "bspline_pattern", std::string( "0" ) );

        tParameterList.insert( "union_pattern", 6 );
        tParameterList.insert( "working_pattern", 7 );

        //! defines which B-Spline mesh is associated with which lagrange mesh
        tParameterList.insert( "lagrange_to_bspline", std::string( "0" ) );

        //! output severity level for moris
        tParameterList.insert( "severity_level", 1 );

        //! boolean for truncated B-Splines
        tParameterList.insert( "truncate_bsplines", 1 );

        //! boolean for multigrid
        tParameterList.insert( "use_multigrid", 0 );

        //! boolean for numbering of aura
        tParameterList.insert( "use_number_aura", 1 );

        //! initial refinement level
        tParameterList.insert( "initial_refinement", 0 );

        //! add comment by the person who implemented this
        tParameterList.insert( "max_refinement_level", -1 );

        //! legacy functions
        tParameterList.insert( "additional_lagrange_refinement", 0 );
        tParameterList.insert( "adaptive_refinement_level", 0 );
        tParameterList.insert( "use_refinement_interrelation", 0 );
        tParameterList.insert( "renumber_lagrange_nodes", 0 );

        return tParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_ */
