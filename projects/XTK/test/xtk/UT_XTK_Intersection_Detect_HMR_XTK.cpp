/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Intersection_Detect_HMR_XTK.cpp
 *
 */

#include "cl_Communication_Tools.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "cl_Param_List.hpp"

#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"

#include "cl_MTK_Mesh_Manager.hpp"

#//include "cl_XTK_Periodic_Boundary_Condition_Helper.hpp"
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

#define protected public
#define private public
#include "cl_XTK_Model.hpp"
#undef protected
#undef private

#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Plane.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

#include "cl_MTK_Intersection_Mesh.hpp"

namespace xtk
{
    TEST_CASE("XTK Intersect","[XTK],[XTK_Intersect_HMR_XTK]")
                {
        if(par_size() ==1)
        {
            uint tLagrangeMeshIndex = 0;

            //HMR parameter list
            moris::ParameterList tParameters = moris::prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 2, 2");
            tParameters.set( "domain_dimensions", "1, 1, 1" );
            tParameters.set( "domain_offset", "0.0, 0.0, 0.0" );
            tParameters.set( "domain_sidesets", "1,2,3,4,5,6" );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );
            tParameters.set( "initial_refinement", "0" );
            tParameters.set( "initial_refinement_pattern", "0" );

            tParameters.set( "use_multigrid", 0 );
            tParameters.set( "severity_level", 0 );

            //create the hmr mesh
            hmr::HMR tHMR( tParameters );

            tHMR.finalize();

            moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);

            //XTK parameter list
            std::string tPeriodicSidePairs= "4,2; 1,3; 5,6";

            moris::ParameterList tXTKParameters = moris::prm::create_xtk_parameter_list();
            tXTKParameters.set( "decompose",                 true );
            tXTKParameters.set( "decomposition_type",        "conformal") ;
            tXTKParameters.set( "enrich",                    true );
            tXTKParameters.set( "basis_rank",                "bspline") ;
            tXTKParameters.set( "enrich_mesh_indices",       "0") ;
            tXTKParameters.set( "ghost_stab",                true );
            tXTKParameters.set( "multigrid",                 false );
            tXTKParameters.set( "verbose",                   true );
            tXTKParameters.set( "print_enriched_ig_mesh",    false );
            tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
            tXTKParameters.set( "high_to_low_dbl_side_sets", true );
            tXTKParameters.set( "periodic_side_set_pair", tPeriodicSidePairs );

            //define the sphere such that it is non interacting
            Vector<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
            tGeometry(0) = std::make_shared<moris::ge::Sphere>(3,3,3,0.1);

            //define ge engine
            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGeometryEngine(tInterpolationMesh, tGeometryEngineParameters);

            size_t tModelDimension = 3;

            Model tXTKModel(tModelDimension, tInterpolationMesh, &tGeometryEngine);
            tXTKModel.mVerbose  =  false;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE,0);

            tXTKModel.construct_face_oriented_ghost_penalization_cells();

            // get meshes
            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tMeshManager->register_mesh_pair(&tEnrInterpMesh, &tEnrIntegMesh);

            // Apply the intersection algorithm
            tXTKModel.mIntersectionDetect = new moris::mtk::Intersection_Detect( tMeshManager, 0, tXTKParameters, 2 ) ;
            tXTKModel.mIntersectionDetect->perform();
            tXTKModel.mIntersectionDetect2D=nullptr;

            // Output the IG mesh
            moris::mtk::Integration_Mesh* tIntegrationMesh2 = tMeshManager->get_integration_mesh( 0 );

            // Get the periodic set
            moris_index tPeriodicSetIndex = tIntegrationMesh2->get_set_index_by_name("P11");
            moris::mtk::Set* tPeriodicSet = tIntegrationMesh2->get_set_by_index( tPeriodicSetIndex );

            // Get clusters on the set
            Vector<moris::mtk::Cluster const *> tClusters = tPeriodicSet->get_clusters_on_set();

            // Cells to populate and check later
            Vector<moris_id > tMorisIdCellDiff;
            Vector<moris_index > tMorisIndexCellDiff;

            //loop over all double sided clusters
            for(uint i = 0 ; i < 12 ; i++)
            {
                // leader side measures
                real tVolumeM = tClusters(i)->compute_cluster_cell_measure(mtk::Primary_Void::PRIMARY , mtk::Leader_Follower::LEADER) ;
                real tLengthM = tClusters(i)->compute_cluster_cell_side_measure(mtk::Primary_Void::PRIMARY , mtk::Leader_Follower::LEADER);

                //follower side measure
                moris::real tVolumeS = tClusters(i)->compute_cluster_cell_measure(mtk::Primary_Void::PRIMARY , mtk::Leader_Follower::FOLLOWER) ;
                moris::real tLengthS = tClusters(i)->compute_cluster_cell_side_measure(mtk::Primary_Void::PRIMARY , mtk::Leader_Follower::FOLLOWER);

                // Measures should be the same
                REQUIRE( tVolumeM == tVolumeS );
                REQUIRE( tLengthM == tLengthS );

                //loop over the integration cells in each double sided cluster
                for( uint k = 0; k < 4 ; k++ )
                {
                    // Get Vertices on the leader and follower
                    Vector< moris::mtk::Vertex const * > tVertex1 = tPeriodicSet->get_clusters_by_index( i )->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER )( k )->get_vertices_on_side_ordinal( 3 );
                    Vector< moris::mtk::Vertex const * > tVertex2 = tPeriodicSet->get_clusters_by_index( i )->get_primary_cells_in_cluster( mtk::Leader_Follower::FOLLOWER )( k )->get_vertices_on_side_ordinal( 3 );

                    // Side ordinal  of the leader and follower cell
                    Matrix< IndexMat > tSideOrdinal2 = tPeriodicSet->get_clusters_by_index( i )->get_cell_side_ordinals(mtk::Leader_Follower::FOLLOWER);
                    Matrix< IndexMat > tSideOrdinal1 = tPeriodicSet->get_clusters_by_index( i )->get_cell_side_ordinals(mtk::Leader_Follower::LEADER);

                    //matrix to store IDs of leader side set
                    Matrix<IdMat> tVertex1ID = Matrix<IdMat> ( 1, 3);
                    Matrix<IdMat> tVertex2ID = Matrix<IdMat> ( 1, 3);

                    //fill in values of IDs
                    for(uint j = 0; j < tVertex1.size(); j++)
                    {
                        tVertex1ID( 0, j) = tVertex1( j )->get_id();
                        tVertex2ID( 0, j) = tVertex2( j )->get_id();
                    }

                    //loop over the vertices on the dide
                    for(uint j = 0 ; j < 3 ; j++ )
                    {
                        moris::mtk::Vertex const *
                        tLeader_Vertex = tClusters(i)->get_leader_vertex_pair(tVertex1(j));

                        moris_id tIdDiff  = tLeader_Vertex->get_id() -  tVertex1(j)->get_id();
                        tMorisIdCellDiff.push_back(tIdDiff);

                        moris_index tVertexOrd = tClusters(i)->get_follower_vertex_ord_on_facet(k,tVertex2(j));

                        if( j == 0 )
                        {
                            tMorisIndexCellDiff.push_back(tVertexOrd == 0 );
                        }
                        else if (j ==1 )
                        {
                            tMorisIndexCellDiff.push_back(tVertexOrd == 1 );
                        }
                        else
                        {
                            tMorisIndexCellDiff.push_back( tVertexOrd == 2 );
                        }
                    }

                }
            }

            for( uint i = 0 ; i < tMorisIndexCellDiff.size() ; i++ )
            {
                REQUIRE( tMorisIndexCellDiff(i) == 1) ;
            }

            for( uint i = 0 ; i < tMorisIdCellDiff.size() ; i++ )
            {

                REQUIRE( tMorisIdCellDiff(i) == 11) ;
            }

        }
                }

}

