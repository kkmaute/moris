/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_HMR_Parameters.hpp
 *
 */

#pragma once

#include "cl_Submodule_Parameter_Lists.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    inline Parameter_List
    create_hmr_parameter_list()
    {
        Parameter_List tParameterList( "General" );

        // number of elements per direction in overall mesh, without aura
        // 2D or 3D is determined by length of this vector
        tParameterList.insert( "number_of_elements_per_dimension", "2, 2" );

        // Processor Decomposition Method (0=user defined; 1=min proc interface; 2=min mesh interface)
        tParameterList.insert( "processor_decomposition_method", 1 );

        // User defined processor grid.  Decomp method must = 0.  Product of array must match number of processors used
        tParameterList.insert( "processor_dimensions", "2, 2" );

        // width, height and depth of domain (without aura)
        tParameterList.insert( "domain_dimensions", "1, 1" );
        // offset from the origin
        tParameterList.insert( "domain_offset", "0, 0 " );
        // sidesets which should be built
        tParameterList.insert( "domain_sidesets", "" );

        // Lagrange Meshes that are used as output meshes
        tParameterList.insert( "lagrange_output_meshes", "" );

        // Lagrange Mesh Names of output meshes
        tParameterList.insert( "lagrange_output_mesh_names", "" );
        // Lagrange Meshes that are used as input meshes
        tParameterList.insert( "lagrange_input_meshes", "" );

        // size of refinement buffer
        tParameterList.insert( "refinement_buffer", 0 );
        // size of staircase buffer
        tParameterList.insert( "staircase_buffer", 0 );

        // Lagrange orders
        tParameterList.insert( "lagrange_orders", "1" );
        // Lagrange pattern
        tParameterList.insert( "lagrange_pattern", "0" );

        // B-Spline orders
        tParameterList.insert( "bspline_orders", "1" );
        // B-Spline orders
        tParameterList.insert( "bspline_pattern", "0" );

        tParameterList.insert( "union_pattern", 6 );
        tParameterList.insert( "working_pattern", 7 );

        // defines which B-Spline mesh is associated with which lagrange mesh
        tParameterList.insert( "lagrange_to_bspline", "0" );

        // output severity level for moris
        tParameterList.insert( "severity_level", 0 );

        // boolean for truncated B-Splines
        tParameterList.insert( "truncate_bsplines", 1 );

        // boolean for multigrid
        tParameterList.insert( "use_multigrid", 0 );

        // boolean for numbering of aura
        tParameterList.insert( "use_number_aura", 1 );

        // initial refinement level
        tParameterList.insert( "initial_refinement", "0" );

        // initial refinement level
        tParameterList.insert( "initial_refinement_pattern", "0" );

        // label of background mesh output file
        tParameterList.insert( "write_background_mesh", "" );

        // label of lagrange mesh output file (VTK)
        tParameterList.insert( "write_lagrange_output_mesh", "" );

        // label of lagrange mesh output file (Exodus)
        tParameterList.insert( "write_lagrange_output_mesh_to_exodus", "" );

        // name of restart file - write
        tParameterList.insert( "write_refinement_pattern_file", false );

        // name of restart file - load
        tParameterList.insert( "restart_refinement_pattern_file", "" );

        // name of vtk file for writing basis function locations
        tParameterList.insert( "basis_function_vtk_file", "" );

        // add comment by the person who implemented this
        tParameterList.insert( "max_refinement_level", -1 );

        // legacy functions
        tParameterList.insert( "additional_lagrange_refinement", 0 );
        tParameterList.insert( "adaptive_refinement_level", 0 );
        tParameterList.insert( "use_refinement_interrelation", 0 );
        tParameterList.insert( "renumber_lagrange_nodes", 0 );
        tParameterList.insert( "use_advanced_T_matrix_scheme", 0 );

        tParameterList.insert( "refinement_function_names", "" );

        // Expert functionality. This function is only for developers. It is not tested in all use cases and will not work in all use cases. T
        // When using this function the user has to know the limitations and unexpected behaviors
        tParameterList.insert( "use_refine_low_level_elements", false );

        return tParameterList;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::prm
