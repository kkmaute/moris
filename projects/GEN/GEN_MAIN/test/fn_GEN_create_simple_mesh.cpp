#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_PRM_HMR_Parameters.hpp"

namespace moris
{
    namespace ge
    {
        mtk::Interpolation_Mesh* create_simple_mesh()
        {
            ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 2");
            tParameters.set( "domain_dimensions", "2, 2");
            tParameters.set( "domain_offset", "-1.0, -1.0");
            tParameters.set( "domain_sidesets", "1,2,3,4");
            tParameters.set( "lagrange_output_meshes", "0");

            tParameters.set( "lagrange_orders", "1");
            tParameters.set( "lagrange_pattern", "0");
            tParameters.set( "bspline_orders", "2");
            tParameters.set( "bspline_pattern", "0");

            tParameters.set( "lagrange_to_bspline", "0");

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 3 );
            tParameters.set( "staircase_buffer", 3 );

            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );

            // initial refinement
            tHMR.perform_initial_refinement( 0 );
            tHMR.finalize();

            return tHMR.create_interpolation_mesh(0);
        }
    }
}
