/*
 * cl_HMR_Integration_Mesh.cpp
 *
 *  Created on: May 30, 2019
 *      Author: doble
 */
#include "catch.hpp"

#include "cl_HMR.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_HMR_Field.hpp"
#include "fn_norm.hpp"

namespace moris
{
moris::real
LevelSetFunction( const moris::Matrix< moris::DDRMat > & aPoint )
{
    return norm( aPoint ) - 0.5;
}
TEST_CASE( "HMR Integration Mesh" , "[IG_Mesh]")
{
    //------------------------------------------------------------------------------

    moris::uint tLagrangeOrder = 1;

    hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

    tParameters.set( "number_of_elements_per_dimension", "10, 4, 4" );
    tParameters.set( "domain_dimensions", "10, 4, 4" );
    tParameters.set( "domain_offset", "-10.0, -2.0, -2.0" );
    tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");
    tParameters.set( "verbose", 0 );
    tParameters.set( "truncate_bsplines", 1 );
    tParameters.set( "bspline_orders", "1" );
    tParameters.set( "lagrange_orders", "1" );

    tParameters.set( "use_multigrid", 0 );

    tParameters.set( "refinement_buffer", 1 );
    tParameters.set( "staircase_buffer", 1 );

    hmr::HMR tHMR( tParameters );

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

    // create field
    std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tLagrangeOrder );

    for( uint k=0; k<3; ++k )
    {
        tField->evaluate_scalar_function( LevelSetFunction );
        tHMR.flag_surface_elements( tField );
        tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
        tHMR.update_refinement_pattern();
    }

    tHMR.finalize();

    // create pointer to output mesh
    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tOutputInterpMesh = tHMR.create_interpolation_mesh( tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern() );

    // create pointer to output mesh
    std::shared_ptr< hmr::Integration_Mesh_HMR > tOutputIntegMesh = tHMR.create_integration_mesh( tLagrangeOrder, tHMR.mParameters->get_lagrange_output_pattern(), *tOutputInterpMesh );

    moris::Cell<std::string> tBlockNames = tOutputIntegMesh->get_block_set_names();

    moris::Cell<moris::mtk::Cluster const *> tCellClustersInBlock = tOutputIntegMesh->get_cell_clusters_in_set(0);

    CHECK(tBlockNames.size() == 1);
    CHECK(tBlockNames(0).compare("HMR_dummy")  == 0);
    CHECK(tCellClustersInBlock.size() == tOutputIntegMesh->get_num_elems());
    CHECK(tOutputInterpMesh->get_num_elems() == tOutputIntegMesh->get_num_elems());

    uint tSideNames = tOutputIntegMesh->get_num_side_sets();
    CHECK(tSideNames == 6);

}
}

