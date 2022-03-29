/*
 * fn_PRM_FEM_Parameters.hpp
 *
 *  Created on: Feb 6, 2020
 *      Author: schmidt
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_HMR_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_HMR_PARAMETERS_HPP_

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
        inline ParameterList
        create_hmr_parameter_list()
        {
            ParameterList tParameterList;

            //! number of elements per direction in overall mesh, without aura
            //! 2D or 3D is determined by length of this vector
            tParameterList.insert( "number_of_elements_per_dimension", std::string( "2, 2" ) );

            //! Processor Decomposition Method (0=user defined; 1=min proc interface; 2=min mesh interface)
            tParameterList.insert( "processor_decomposition_method", 1 );

            //! User defined processor grid.  Decomp method must = 0.  Product of array must match number of processors used
            tParameterList.insert( "processor_dimensions", std::string( "2, 2" ) );

            //! width, height and depth of domain (without aura)
            tParameterList.insert( "domain_dimensions", std::string( "1, 1" ) );
            //! offset from the origin
            tParameterList.insert( "domain_offset", std::string( "0, 0 " ) );
            //! sidesets which should be built
            tParameterList.insert( "domain_sidesets", std::string( "" ) );

            //! Lagrange Meshes that are used as output meshes
            tParameterList.insert( "lagrange_output_meshes", std::string( "" ) );
            //! Lagrange Meshe Names of output meshes
            tParameterList.insert( "lagrange_output_meshe_names", std::string( "" ) );
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
            tParameterList.insert( "initial_refinement", "0" );

            //! initial refinement level
            tParameterList.insert( "initial_refinement_pattern", "0" );

            //! label of background mesh output file
            tParameterList.insert( "write_background_mesh", std::string( "" ) );

            //! label of lagrange mesh output file (VTK)
            tParameterList.insert( "write_lagrange_output_mesh", std::string( "" ) );

            //! label of lagrange mesh output file (Exodus)
            tParameterList.insert( "write_lagrange_output_mesh_to_exodus", std::string( "" ) );

            //! name of restart file - write
            tParameterList.insert( "write_refinement_pattern_file", true );

            //! name of restart file - load
            tParameterList.insert( "restart_refinement_pattern_file", std::string( "" ) );

            //! add comment by the person who implemented this
            tParameterList.insert( "max_refinement_level", -1 );

            //! legacy functions
            tParameterList.insert( "additional_lagrange_refinement", 0 );
            tParameterList.insert( "adaptive_refinement_level", 0 );
            tParameterList.insert( "use_refinement_interrelation", 0 );
            tParameterList.insert( "renumber_lagrange_nodes", 0 );
            tParameterList.insert( "use_advanced_T_matrix_scheme", 0 );

            tParameterList.insert( "refinement_function_names", "" );

            //! Expert functionality. This function is only for developers. It is not tested in all use cases and will not work in all use cases. T
            //! When using this function the user has to know the limitations and unexpected behaviors
            tParameterList.insert( "use_refine_low_level_elements", false );

            return tParameterList;
        }
        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif /* PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_ */
