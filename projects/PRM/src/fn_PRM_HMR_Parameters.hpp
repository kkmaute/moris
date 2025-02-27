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

#include "cl_Parameter_List.hpp"

namespace moris::prm
{

    //------------------------------------------------------------------------------

    // creates a parameter list with default inputs
    /**
     * Creates an HMR general parameter list with default inputs.
     *
     * @return HMR general parameter list
     */
    inline Parameter_List
    create_hmr_parameter_list()
    {
        Parameter_List tParameterList( "General" );

        // width, height and depth of domain (without aura)
        tParameterList.insert( "domain_dimensions", Vector< real >{ 1.0, 1.0 } );

        // offset from the origin
        tParameterList.insert( "domain_offset", Vector< real >() );

        // number of elements per direction in overall mesh, without aura
        // 2D or 3D is determined by length of this vector
        tParameterList.insert( "number_of_elements_per_dimension", Vector< uint >() );

        // Processor Decomposition Method (0=user defined; 1=min proc interface; 2=min mesh interface)
        tParameterList.insert( "processor_decomposition_method", 1 );

        // User defined processor grid.  Decomp method must = 0.  Product of array must match number of processors used
        tParameterList.insert( "processor_dimensions", Vector< uint >() );

        // Lagrange Meshes that are used as output meshes
        tParameterList.insert( "lagrange_output_meshes", "" );

        // Lagrange Mesh Names of output meshes
        tParameterList.insert( "lagrange_output_mesh_names", "" );

        // size of refinement buffer
        tParameterList.insert( "refinement_buffer", 0u );
        // size of staircase buffer
        tParameterList.insert( "staircase_buffer", 0u );

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

        // output severity level for moris TODO this should not be a part of HMR
        tParameterList.insert( "severity_level", 0 );

        // boolean for truncated B-Splines
        tParameterList.insert( "truncate_bsplines", true );

        // boolean for multigrid
        tParameterList.insert( "use_multigrid", false );

        // boolean for numbering of aura
        tParameterList.insert( "use_number_aura", true );

        // initial refinement level per pattern. Each entry creates a new pattern
        tParameterList.insert( "initial_refinement", Vector< uint >{ 0 } );

        // label of background mesh output file
        tParameterList.insert( "write_background_mesh", "" );

        // label of lagrange mesh output file (VTK/Exodus)
        tParameterList.insert( "lagrange_mesh_output_file_name", "" );

        // name of restart file - write
        tParameterList.insert( "write_refinement_pattern_file", false );

        // name of restart file - load
        tParameterList.insert( "restart_refinement_pattern_file", "" );

        // name of vtk file for writing basis function locations
        tParameterList.insert( "basis_function_vtk_file", "" );

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

    /**
     * Creates an HMR Lagrange mesh parameter list.
     *
     * @return B-spline mesh parameter list
     */
    inline Parameter_List
    create_lagrange_mesh_parameter_list()
    {
        Parameter_List tParameterList( "Lagrange Mesh" );

        tParameterList.insert( "pattern_index", 0u ); // Pattern to use for this Lagrange mesh
        tParameterList.insert( "order", 1u ); // Lagrange order
        tParameterList.insert( "is_output_mesh", true ); // If this is an output mesh or not
        tParameterList.insert( "output_mesh_name", "" ); // Custom name for this output mesh

        return tParameterList;
    }

    /**
     * Creates an HMR B-spline mesh parameter list.
     *
     * @return B-spline mesh parameter list
     */
    inline Parameter_List
    create_bspline_mesh_parameter_list()
    {
        Parameter_List tParameterList( "B-spline Mesh" );

        tParameterList.insert( "pattern_index", 0u ); // Pattern to use for this B-spline mesh
        tParameterList.insert( "orders", Vector< uint >{ 1 } ); // B-spline orders (single order, or order per dimension)
        tParameterList.insert( "paired_lagrange_mesh_index", 0u ); // Lagrange mesh index to assign this B-spline mesh to

        return tParameterList;
    }

}    // namespace moris::prm
