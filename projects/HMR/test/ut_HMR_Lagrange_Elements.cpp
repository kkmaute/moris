/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Lagrange_Elements.cpp
 *
 */

#include "catch.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Hex27.hpp"
#include "cl_MTK_Cell_Info_Hex64.hpp"

namespace moris
{
TEST_CASE("Single Hex 8 Lagrange Mesh","[Lag_Hex8]")
{
    if(par_size() == 1)
    {
        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {1} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {1} });
        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-1.0} });
        tParameters.set_bspline_truncation( true );
//        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {2} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 1 );

        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//        tHMR.save_to_exodus( tLagrangeMeshIndex, "hex8_hmr.exo" );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);

        //Check that the vertices on side are correct
        moris::mtk::Cell_Info_Hex8 tConn;
        moris::Matrix<moris::IndexMat> tVertexToSideMap = tConn.get_node_to_face_map();

        moris::Cell< mtk::Vertex* > tCellVertices = tCell.get_vertex_pointers();

        // verify vertices on side ordinals
        for(moris::uint  i = 0; i <6; i++)
        {
            moris::Cell< mtk::Vertex const * > tVertsOnSide = tCell.get_vertices_on_side_ordinal(i);

            for(moris::uint j = 0; j <4; j ++)
            {
                //Expected vertex
                mtk::Vertex* tExpVert = tCellVertices(tVertexToSideMap(i,j));

                // received vertex
                mtk::Vertex const * tRecVert = tVertsOnSide(j);

                CHECK(tExpVert->get_id() == tRecVert->get_id());
            }
        }

        tHMR.finalize();
    }
}

TEST_CASE("Single Hex 27 Lagrange Mesh","[Lag_Hex27]")
{
    if(par_size() == 1)
    {
        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {1} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {1} });
        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-1.0} });
        tParameters.set_bspline_truncation( true );
//        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set_lagrange_orders  ( { {2} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {2} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 1 );

        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

//        tHMR.save_to_exodus( tLagrangeMeshIndex, "hex27_hmr.exo" );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);

        //Check that the vertices on side are correct
        moris::mtk::Cell_Info_Hex27 tConn;
        moris::Matrix<moris::IndexMat> tVertexToSideMap = tConn.get_node_to_face_map();

        moris::Cell< mtk::Vertex* > tCellVertices = tCell.get_vertex_pointers();

        // verify vertices on side ordinals
        for(moris::uint  i = 0; i <6; i++)
        {
            moris::Cell< mtk::Vertex const * > tVertsOnSide = tCell.get_vertices_on_side_ordinal(i);

            CHECK(tVertsOnSide.size() == 9);

            for(moris::uint j = 0; j <9; j ++)
            {
                //Expected vertex
                mtk::Vertex* tExpVert = tCellVertices(tVertexToSideMap(i,j));

                // received vertex
                mtk::Vertex const * tRecVert = tVertsOnSide(j);

                CHECK(tExpVert->get_id() == tRecVert->get_id());
            }
        }

        tHMR.finalize();
    }
}

TEST_CASE("Single Hex 64 Lagrange Mesh","[Lag_Hex64]")
{
    if(par_size() == 1)
    {
        moris::uint tLagrangeMeshIndex = 0;

        moris::hmr::Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { {1}, {1}, {1} } );
        tParameters.set_domain_dimensions({ {1}, {1}, {1} });
        tParameters.set_domain_offset({ {-1.0}, {-1.0}, {-1.0} });
        tParameters.set_bspline_truncation( true );
//        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set_lagrange_orders  ( { {3} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {2} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 1 );

        Cell< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);

        //Check that the vertices on side are correct
        moris::mtk::Cell_Info_Hex64 tConn;
        moris::Matrix<moris::IndexMat> tVertexToSideMap = tConn.get_node_to_face_map();

        moris::Cell< mtk::Vertex* > tCellVertices = tCell.get_vertex_pointers();

        // verify vertices on side ordinals
        for(moris::uint  i = 0; i <6; i++)
        {
            moris::Cell< mtk::Vertex const * > tVertsOnSide = tCell.get_vertices_on_side_ordinal(i);

            CHECK(tVertsOnSide.size() == 16);

            for(moris::uint j = 0; j <  16; j ++)
            {
                //Expected vertex
                mtk::Vertex* tExpVert = tCellVertices(tVertexToSideMap(i,j));

                // received vertex
                mtk::Vertex const * tRecVert = tVertsOnSide(j);

                CHECK(tExpVert->get_id() == tRecVert->get_id());
            }
        }

        tHMR.finalize();
    }
}

}

