/*
 * cl_MTK_Blocks.cpp
 *
 *  Created on: May 28, 2019
 *      Author: schmidt
 */

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "cl_Param_List.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Set.hpp" //MTK/src
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "catch.hpp"
#include "paths.hpp"

// implementations to test
#include "cl_MTK_Mesh_Factory.hpp"

namespace moris
{
    namespace mtk
    {
        TEST_CASE("MTK Periodic","[MTK],[MTK_Periodic]")
        {
            if(par_size() ==1)
            {
                //generate a cubic mesh
                std::string tInterpString = "generated:2x1x1|sideset:xX";

                //create interpolation integration mesh
                mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tInterpString );
                mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );

                // create a mesh manager and register the mesh pair
                auto tMeshManager = std::make_shared< mtk::Mesh_Manager >();
                tMeshManager->register_mesh_pair( tInterpMesh, tIntegMesh );

                //parameter list input for the surfaces that are going to be periodic
                moris::ParameterList tParameterList;
                tParameterList.insert( "periodic_side_set_pair", "surface_1,surface_2");

                //construct the boundary condition helper
                mtk::Periodic_Boundary_Condition_Helper tPBCHelper(tMeshManager,0, tParameterList);

                //perform the periodic boundary condition
                tPBCHelper.setup_periodic_boundary_conditions();

                //check number of double sided sets
                REQUIRE( tIntegMesh->get_num_double_side_set() == 1 );

                //FIXME -- implemented in STK
                //REQUIRE( tIntegMesh->get_num_double_sided_sets() == 1 );

                // get the periodic side set by index
                mtk::Set * tSideSet1 = tIntegMesh->get_set_by_index( 5 );

                // testing get set by name function
                REQUIRE( tSideSet1->get_set_name() == "Periodic0" );

                // recover index of the periodic dbl sided set
                REQUIRE(tIntegMesh->get_set_index_by_name("Periodic0") == 5 );

                //get the side ordinal that has PBC
                Matrix< IndexMat > tSideOrdinal1 = tSideSet1->get_clusters_by_index( 0 )
                                                                   ->get_cell_side_ordinals( mtk::Master_Slave::MASTER );
                //get the vertices on the PBC master side
                moris::Cell< moris::mtk::Vertex const* > tVertex1 = tSideSet1->get_clusters_by_index( 0 )
                                                                      ->get_primary_cells_in_cluster( mtk::Master_Slave::MASTER )( 0 )->get_vertices_on_side_ordinal( tSideOrdinal1( 0, 0) );

                //matrix to store IDs of master side set
                Matrix<IdMat> tVertex1ID = Matrix<IdMat> ( 1, 4);

                //fill in values of IDs
                for(uint j = 0; j < tVertex1.size(); j++)
                {
                    tVertex1ID( 0, j) = tVertex1( j )->get_id();
                }

                //get the side ordinal that has PBC
                Matrix< IndexMat > tSideOrdinal2 = tSideSet1->get_clusters_by_index( 0 )->get_cell_side_ordinals(mtk::Master_Slave::SLAVE);

                //get the vertices on the PBC master side
                moris::Cell< moris::mtk::Vertex const* > tVertex2 = tSideSet1->get_clusters_by_index( 0 )
                                                                                     ->get_primary_cells_in_cluster( mtk::Master_Slave::SLAVE )( 0 )->get_vertices_on_side_ordinal( tSideOrdinal2( 0, 0));

                //matrix to store IDs of master side set
                Matrix<IdMat> tVertex2ID = Matrix<IdMat> (1,4);

                //fill in values of IDs
                for(uint j=0; j<tVertex2.size(); j++ )
                {
                    tVertex2ID(0,j)=tVertex2(j)->get_id();
                }

                REQUIRE( tVertex1ID( 0, 0) == 1 );                   REQUIRE( tVertex2ID( 0, 0) == 3 );
                REQUIRE( tVertex1ID( 0, 1) == 7 );                   REQUIRE( tVertex2ID( 0, 1) == 6 );
                REQUIRE( tVertex1ID( 0, 2) == 10);                   REQUIRE( tVertex2ID( 0, 2) == 12 );
                REQUIRE( tVertex1ID( 0, 3) == 4 );                   REQUIRE( tVertex2ID( 0, 3) == 9 );

                //matrix to store IDs of master side set
                Matrix<IdMat> tVertexPairID = Matrix<IdMat> (  1, 4);

                for(uint k = 0;  k < tVertex2.size(); k++)
                {
                    moris::mtk::Vertex const* tVertex = tSideSet1->get_clusters_by_index( 0 )->get_master_vertex_pair( tVertex1(k) );
                    tVertexPairID( 0, k ) = tVertex->get_id();
                }

                REQUIRE( tVertexPairID( 0, 0 ) == 3 );
                REQUIRE( tVertexPairID( 0, 1 ) == 9 );
                REQUIRE( tVertexPairID( 0, 2 ) == 12 );
                REQUIRE( tVertexPairID( 0, 3 ) == 6 );

                //delete tInterpMesh;
                delete tIntegMesh;
            }
        }
    }
}
