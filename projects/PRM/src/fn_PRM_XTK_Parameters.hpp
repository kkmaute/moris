/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_PRM_XTK_Parameters.hpp
 *
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_XTK_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_XTK_PARAMETERS_HPP_

#include "cl_Param_List.hpp"

namespace moris
{
    namespace prm
    {

        //------------------------------------------------------------------------------

        // creates a parameter list with default inputs
        inline ParameterList
        create_xtk_parameter_list()
        {
            ParameterList tParameterList;

            // decomposition and decomposition related parameters
            tParameterList.insert( "decompose", true );
            tParameterList.insert( "decomposition_type", "conformal" );
            tParameterList.insert( "octree_refinement_level", "-1" );
            tParameterList.insert( "triangulate_all", false );    // NOTE: this option does fail if the Lagrange mesh is not uniformly refined
            tParameterList.insert( "ig_element_order", moris::uint( 1 ) );

            // cleanup
            tParameterList.insert( "cleanup_cut_mesh", false );

            // enrichment and enrichment related parameters
            tParameterList.insert( "enrich", false );
            tParameterList.insert( "use_SPG_based_enrichment", false );
            tParameterList.insert( "basis_rank", "node" );
            tParameterList.insert( "enrich_mesh_indices", "0" );
            tParameterList.insert( "sort_basis_enrichment_levels", false );
            tParameterList.insert( "unenriched_mesh_indices", "" );

            // ghost stabilization and ghost related parameters
            tParameterList.insert( "ghost_stab", false );         // Perform ghost stabilization
            tParameterList.insert( "visualize_ghost", false );    // writes the ghost blocks, on the XTK output mesh

            // multigrid
            tParameterList.insert( "multigrid", false );

            // contact sandbox
            tParameterList.insert( "contact_sandbox", false );
            tParameterList.insert( "potential_phases_in_contact", "" );
            tParameterList.insert( "bb_epsilon", 0.1 );

            // verbose - should be replaced by the severity level of the logger
            tParameterList.insert( "verbose", false );
            tParameterList.insert( "verbose_level", moris::uint( 0 ) );    // 0 - basic outputs // 1- lots of outputs
            tParameterList.insert( "diagnostics", false );
            tParameterList.insert( "diagnostics_id", "" );
            tParameterList.insert( "diagnostics_path", "" );

            // if to deactivate empty sets - used only if outputting ig mesh as well, set to true only for debugging
            tParameterList.insert( "deactivate_empty_sets", false );

            // deactivate all but selected
            tParameterList.insert( "deactivate_all_but_blocks", "" );
            tParameterList.insert( "deactivate_all_but_side_sets", "" );

            // Write enrichment fields on mesh (only recommended on very small meshes)
            tParameterList.insert( "write_enrichment_fields", false );
            tParameterList.insert( "write_basis_functions", false );

            // request B-spline cluster information (SPG IDs, indices, etc.) to be written to the xtk_temp.exo
            tParameterList.insert( "write_bspline_cluster_info", false );

            // a sphere where I write enrichment fields locations (r,xc,yv,zc)
            tParameterList.insert( "write_enrichment_fields_probe_spheres", "" );

            // T-Matrix output if needed
            tParameterList.insert( "global_T_matrix_output_file", "" );
            tParameterList.insert( "elemental_T_matrix_output_file", "" );
            tParameterList.insert( "MPC_output_file", "" );

            // triangulate mesh at the end of constructing, this allows non-uniform refinement in the Lagrange mesh
            // NOTE: this option is only intended for mesh output, the produced mesh may miss some domain boundary side sets
            tParameterList.insert( "triangulate_all_in_post", false );

            // kill workflow after outputting XTK mesh, used for quickly testing the geometry and assign material phases
            tParameterList.insert( "only_generate_xtk_temp", false );

            // path to folder for XTK output
            tParameterList.insert( "output_path", "./" );
            tParameterList.insert( "output_file", "xtk_temp.exo" );

            // print enriched integration meshes
            tParameterList.insert( "print_enriched_ig_mesh", false );

            // save enriched integration mesh for every optimization iteration
            tParameterList.insert( "keep_all_opt_iters", false );

            // write XTK exodus mesh
            tParameterList.insert( "exodus_output_XTK_ig_mesh", false );
            tParameterList.insert( "exodus_output_XTK_ip_mesh", false );

            // write Cut_Integration_Mesh
            tParameterList.insert( "output_cut_ig_mesh", false );

            // write Intersection_Mesh
            tParameterList.insert( "output_intersection_mesh", false );

            // enriched integration mesh options
            tParameterList.insert( "high_to_low_dbl_side_sets", false );

            // print memory usage
            tParameterList.insert( "print_memory", false );
            tParameterList.insert( "low_memory", true );

            // probe a cell - Debug
            tParameterList.insert( "probe_bg_cells", "" );

            // union
            tParameterList.insert( "union_blocks", "" );
            tParameterList.insert( "union_block_names", "" );
            tParameterList.insert( "union_block_colors", "" );

            // union
            tParameterList.insert( "union_side_sets", "" );
            tParameterList.insert( "union_side_set_names", "" );
            tParameterList.insert( "union_side_set_colors", "" );

            tParameterList.insert( "identify_hanging_nodes", false );

            tParameterList.insert( "delete_xtk_after_generation", true );

            tParameterList.insert( "activate_basis_agglomeration", false );
            tParameterList.insert( "volume_fraction", 1.0 );    // by default all the cut cells are considered bad
            tParameterList.insert( "activate_cell_agglomeration", false );
            tParameterList.insert( "visualize_cell_association", true );

            return tParameterList;
        }
        //------------------------------------------------------------------------------

    }    // namespace prm
}    // namespace moris

#endif /* PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_ */
