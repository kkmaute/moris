/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_GEN_create_simple_mesh.cpp
 *
 */

#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::gen
{
    mtk::Interpolation_Mesh* create_simple_mesh(
            uint aNumXElements,
            uint aNumYElements,
            uint aLagrangeOrder,
            uint aBSplineOrder,
            uint aRefinement)
    {
        ParameterList tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension",
                         std::to_string(aNumXElements) + ", " + std::to_string(aNumYElements));
        tParameters.set( "domain_dimensions", "2, 2");
        tParameters.set( "domain_offset", "-1.0, -1.0");
        tParameters.set( "domain_sidesets", "1,2,3,4");
        tParameters.set( "lagrange_output_meshes", "0");

        tParameters.set( "lagrange_orders", std::to_string(aLagrangeOrder));
        tParameters.set( "lagrange_pattern", "0");
        tParameters.set( "bspline_orders", std::to_string(aBSplineOrder));
        tParameters.set( "bspline_pattern", "0");
        tParameters.set( "lagrange_to_bspline", "0");

        tParameters.set( "initial_refinement", std::to_string(aRefinement) );
        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );

        tParameters.set( "severity_level", 2 );

        hmr::HMR tHMR( tParameters );

        // initial refinement
        tHMR.perform_initial_refinement();
        tHMR.finalize();

        return tHMR.create_interpolation_mesh(0);
    }
}
