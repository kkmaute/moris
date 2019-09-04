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

    moris::uint tLagrangeMeshIndex = 0;
    moris::uint tBSplineMeshIndex = 0;

    moris::hmr::Parameters tParameters;

    tParameters.set_number_of_elements_per_dimension( { {10}, {4}, {4} } );
    tParameters.set_domain_dimensions({ {10}, {4}, {4} });
    tParameters.set_domain_offset({ {-10.0}, {-2.0}, {-2.0} });
    tParameters.set_bspline_truncation( true );
    tParameters.set_side_sets({ {1}, {6}, {3}, {4}, {5}, {2} });

    tParameters.set_output_meshes( { {0} } );

    tParameters.set_lagrange_orders  ( { {1} });
    tParameters.set_lagrange_patterns({ {0} });

    tParameters.set_bspline_orders   ( { {1} } );
    tParameters.set_bspline_patterns ( { {0} } );

    tParameters.set_union_pattern( 2 );
    tParameters.set_working_pattern( 3 );

    tParameters.set_refinement_buffer( 1 );
    tParameters.set_staircase_buffer( 1 );

    Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
    tLagrangeToBSplineMesh( 0 ) = { {0} };

    tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

    hmr::HMR tHMR( tParameters );

    std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

    // create field
    std::shared_ptr< moris::hmr::Field > tField = tMesh->create_field( "Circle", tBSplineMeshIndex );

    for( uint k=0; k<3; ++k )
    {
        tField->evaluate_scalar_function( LevelSetFunction );
        tHMR.flag_surface_elements( tField );
        tHMR.perform_refinement( 0 );
        tHMR.update_refinement_pattern( 0 );
    }

    tHMR.finalize();

    // create pointer to output mesh
    std::shared_ptr< hmr::Interpolation_Mesh_HMR > tOutputInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

    // create pointer to output mesh
    std::shared_ptr< hmr::Integration_Mesh_HMR > tOutputIntegMesh = tHMR.create_integration_mesh( 1, 0, *tOutputInterpMesh );

    moris::Cell<std::string> tBlockNames = tOutputIntegMesh->get_block_set_names();

    moris::Cell<moris::mtk::Cluster const *> tCellClustersInBlock = tOutputIntegMesh->get_cell_clusters_in_set(0);

    CHECK(tBlockNames.size() == 1);
    CHECK(tBlockNames(0).compare("HMR_dummy")  == 0);
    CHECK(tCellClustersInBlock.size() == tOutputIntegMesh->get_num_elems());
    CHECK(tOutputInterpMesh->get_num_elems() == tOutputIntegMesh->get_num_elems());

    uint tSideNames = tOutputIntegMesh->get_num_side_sets();

    CHECK(tSideNames == 6);
}

TEST_CASE( "HMR_Basis_Support" , "[hmr][HMR_Basis_Support]")
{
    //------------------------------------------------------------------------------

    if(par_size() == 1)
    {
        moris::uint tLagrangeMeshIndex = 0;
        moris::uint tBSplineMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {2}, {2} } );
        tParameters.set_bspline_truncation( true );

        tParameters.set_output_meshes( { {0} } );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {1} } );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        Cell< Matrix< DDUMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        // std::shared_ptr< Database >
        auto tDatabase = tHMR.get_database();

        // manually select output pattern
        tDatabase->set_activation_pattern( 0 );

        for( uint tLevel = 0; tLevel < 1; ++tLevel )
        {
            tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

            // manually refine, do not reset pattern
            tDatabase->get_background_mesh()->perform_refinement(0);
        }

        tDatabase->get_background_mesh()->save_to_vtk("Basis_support.vtk");

        tHMR.finalize();

        auto tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        Matrix< IndexMat > tElementIdnInBasisSupport_1;
        Matrix< IndexMat > tElementIdnInBasisSupport_2;
        Matrix< IndexMat > tElementIdnInBasisSupport_3;

        tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 1, tElementIdnInBasisSupport_1 );
        tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 4, tElementIdnInBasisSupport_2 );
        tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 8, tElementIdnInBasisSupport_3 );

        CHECK(tElementIdnInBasisSupport_1( 0 ) == 0 );    CHECK(tElementIdnInBasisSupport_1( 1 ) == 1 );
        CHECK(tElementIdnInBasisSupport_1( 2 ) == 2 );    CHECK(tElementIdnInBasisSupport_1( 3 ) == 3 );
        CHECK(tElementIdnInBasisSupport_1( 4 ) == 4 );

        CHECK(tElementIdnInBasisSupport_2( 0 ) == 0 );    CHECK(tElementIdnInBasisSupport_2( 1 ) == 1 );
        CHECK(tElementIdnInBasisSupport_2( 2 ) == 2 );    CHECK(tElementIdnInBasisSupport_2( 3 ) == 3 );
        CHECK(tElementIdnInBasisSupport_2( 4 ) == 4 );    CHECK(tElementIdnInBasisSupport_2( 5 ) == 5 );
        CHECK(tElementIdnInBasisSupport_2( 6 ) == 6 );

        CHECK(tElementIdnInBasisSupport_3( 0 ) == 6 );
    }
}

}

