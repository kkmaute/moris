/*
 * cl_PRM_XTK_Parameters.hpp
 *
 *  Created on: March 10, 2020
 *      Author: doble
 */

#ifndef PROJECTS_PRM_SRC_CL_PRM_XTK_PARAMETERS_HPP_
#define PROJECTS_PRM_SRC_CL_PRM_XTK_PARAMETERS_HPP_

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
    ParameterList create_xtk_parameter_list()
    {
        ParameterList tParameterList;

        // decomposition and decomposition related parameters
        tParameterList.insert( "decompose", true );
        tParameterList.insert( "decomposition_type", std::string("conformal") );

        // enrichment and enrichment related parameters
        tParameterList.insert( "enrich", false );
        tParameterList.insert( "basis_rank",std::string("node") );
        tParameterList.insert( "enrich_mesh_indices",std::string("0"));

        // ghost stabilization and ghost related parameters
        tParameterList.insert( "ghost_stab", false );

        // multigrid
        tParameterList.insert( "multigrid", false );

        // verbose - should be replaced by the severity level of the logger
        tParameterList.insert( "verbose", false );

        tParameterList.insert( "print_enriched_ig_mesh", false );

        tParameterList.insert( "exodus_output_XTK_ig_mesh", false );

        return tParameterList;
    }
//------------------------------------------------------------------------------

    }/* end_namespace_prm */
}/* end_namespace_moris */

#endif /* PROJECTS_PRM_SRC_CL_PRM_MSI_PARAMETERS_HPP_ */
