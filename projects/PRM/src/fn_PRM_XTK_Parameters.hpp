/*
 * fn_PRM_XTK_Parameters.hpp
 *
 *  Created on: March 10, 2020
 *      Author: doble
 */

#ifndef PROJECTS_PRM_SRC_FN_PRM_XTK_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_FN_PRM_XTK_PARAMETERS_HPP_

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
    create_xtk_parameter_list()
    {
        ParameterList tParameterList;

        // decomposition and decomposition related parameters
        tParameterList.insert( "decompose", true );
        tParameterList.insert( "decomposition_type", "conformal" );
        tParameterList.insert( "octree_refinement_level", "-1" );
        tParameterList.insert( "triangulate_all", false );

        // cleanup
        tParameterList.insert( "cleanup_cut_mesh", false );

        // enrichment and enrichment related parameters
        tParameterList.insert( "enrich", false );
        tParameterList.insert( "basis_rank", "node" );
        tParameterList.insert( "enrich_mesh_indices", "0" );

        // ghost stabilization and ghost related parameters
        tParameterList.insert( "ghost_stab", false );// Perform ghost stabilization
        tParameterList.insert( "visualize_ghost", false );// writes the ghost blocks, on the XTK output mesh

        // multigrid
        tParameterList.insert( "multigrid", false );

        // contact sandbox
        tParameterList.insert( "contact_sandbox", false );
        tParameterList.insert( "potential_phases_in_contact", "" );
        tParameterList.insert( "bb_epsilon", 0.1 );

        // verbose - should be replaced by the severity level of the logger
        tParameterList.insert( "verbose", false );
        tParameterList.insert( "diagnostics", false );
        tParameterList.insert( "diagnostics_id", "" );
        tParameterList.insert( "diagnostics_path", "" );

        // if to deactivate empty sets - used only if outputting ig mesh as well, set to true only for debugging
        tParameterList.insert( "deactivate_empty_sets", false );

        // deactivate all but selected
        tParameterList.insert( "deactivate_all_but_blocks", "" );
        tParameterList.insert( "deactivate_all_but_side_sets", "" );

        // Write enrichement fields on mesh (only recommended on very small meshes)
        tParameterList.insert( "write_enrichment_fields", false );

        // path to folder for XTK output
        tParameterList.insert( "output_path", "./" );
        tParameterList.insert( "output_file", "xtk_temp.exo" );

        // print enriched integration meshes
        tParameterList.insert( "print_enriched_ig_mesh", false );

        // save enriched integration mesh for every optimization iteration
        tParameterList.insert( "keep_all_opt_iters", false );

        // write XTK exodus mesh
        tParameterList.insert( "exodus_output_XTK_ig_mesh", false );

        // enriched integration mesh options
        tParameterList.insert( "high_to_low_dbl_side_sets", false );

        // print memory usage
        tParameterList.insert( "print_memory", false );

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

        //Periodic environment sidesets
        tParameterList.insert( "periodic_side_set_pair", "" );

        tParameterList.insert( "identify_hanging_nodes", false );

        return tParameterList;
    }
    //------------------------------------------------------------------------------

}// namespace prm
}// namespace moris

#endif /* PROJECTS_PRM_SRC_FN_PRM_MSI_PARAMETERS_HPP_ */
