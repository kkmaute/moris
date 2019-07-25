/*
 * cl_HMR_Lagrange_Elements.cpp
 *
 *  Created on: May 31, 2019
 *      Author: doble
 */


#include "catch.hpp"

#include "cl_HMR.hpp"
#include "cl_MTK_Hex8_Connectivity.hpp"
#include "cl_MTK_Hex27_Connectivity.hpp"
#include "cl_MTK_Hex64_Connectivity.hpp"

namespace moris
{

TEST_CASE("Single Hex 8 Lagrange Mesh","[Lag_Hex8]")
        {
    if(par_size() == 1)
    {
        moris::uint tLagrangeOrder = 1;

        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "1, 1, 1" );
        tParameters.set( "domain_dimensions", "1, 1, 1" );
        tParameters.set( "domain_offset", "-1.0, -1.0, -1.0" );
        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "bspline_orders", "2" );
        tParameters.set( "lagrange_orders", "1" );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 2 );
        tParameters.set( "staircase_buffer", 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

        tHMR.save_to_exodus( "hex8_hmr.exo", 0, tLagrangeOrder );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);
        std::cout<<tCell.get_id()<<std::endl;
        std::cout<<tCell.get_index()<<std::endl;

        //Check that the vertices on side are correct
        moris::Matrix<moris::IndexMat> tVertexToSideMap = Hex8::get_node_to_face_map();

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
        moris::uint tLagrangeOrder = 2;

        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "1, 1, 1" );
        tParameters.set( "domain_dimensions", "1, 1, 1" );
        tParameters.set( "domain_offset", "-1.0, -1.0, -1.0" );
        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "bspline_orders", "2" );
        tParameters.set( "lagrange_orders", "2" );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 2 );
        tParameters.set( "staircase_buffer", 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

        tHMR.save_to_exodus( "hex27_hmr.exo", 0, tLagrangeOrder );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);
        std::cout<<tCell.get_id()<<std::endl;
        std::cout<<tCell.get_index()<<std::endl;

        //Check that the vertices on side are correct
        moris::Matrix<moris::IndexMat> tVertexToSideMap = Hex27::get_node_to_face_map();

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
        moris::uint tLagrangeOrder = 3;

        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "1, 1, 1" );
        tParameters.set( "domain_dimensions", "1, 1, 1" );
        tParameters.set( "domain_offset", "-1.0, -1.0, -1.0" );
        tParameters.set( "domain_sidesets", "1, 6, 3, 4, 5, 2");

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "bspline_orders", "2" );
        tParameters.set( "lagrange_orders", "3" );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 2 );
        tParameters.set( "staircase_buffer", 1 );

        hmr::HMR tHMR( tParameters );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeOrder );

        // get the cells
        moris::mtk::Cell & tCell = tMesh->get_mtk_cell(0);
        std::cout<<tCell.get_id()<<std::endl;
        std::cout<<tCell.get_index()<<std::endl;

        //Check that the vertices on side are correct
        moris::Matrix<moris::IndexMat> tVertexToSideMap = Hex64::get_node_to_face_map();

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
