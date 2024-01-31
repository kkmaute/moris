/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Enrichment.cpp
 *
 */

#include <memory>
#include <mpi.h>

#include "catch.hpp"

// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
#include "fn_verify_tet_topology.hpp"
#include "fn_write_element_ownership_as_field.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes
#include "cl_MTK_Mesh_Checker.hpp"

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_GEN_Mesh_Field.hpp"
#include "cl_GEN_Plane.hpp"

namespace xtk
{

TEST_CASE( "Enrichment Example 1", "[ENRICH_1]" )
{
    if ( par_size() == 1 || par_size() == 1 )
    {
        bool tOutputEnrichmentFields = true;

        // Generate mesh from string
        std::string tMeshFileName = "generated:1x1x3|sideset:XzZ";

        // Add level set field to add onto file

        // Specify field parameters
        moris::mtk::Scalar_Field_Info< DDRMat > tLSF;
        std::string                             tLSFName = "lsf1";
        tLSF.set_field_name( tLSFName );
        tLSF.set_field_entity_rank( mtk::EntityRank::NODE );

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input( &tLSF, tFieldsInfo );

        // add to mesh data input container
        moris::mtk::MtkMeshData tSuppMeshData;
        tSuppMeshData.FieldsInfo = &tFieldsInfo;

        // Create mesh with supplementary data
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tSuppMeshData );

        xtk::size_t tNumNodes = tMeshData->get_num_entities( mtk::EntityRank::NODE );

        moris::Matrix< moris::DDRMat > tLevelsetVal( tNumNodes, 1, -1.3 );

        moris_id tIndexOfNodeId6 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 6, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId3 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3, mtk::EntityRank::NODE );

        // Bottom face
        tLevelsetVal( tIndexOfNodeId3 ) = 1;
        tLevelsetVal( tIndexOfNodeId6 ) = 1;

        tMeshData->add_mesh_field_real_scalar_data_loc_inds( tLSFName, mtk::EntityRank::NODE, tLevelsetVal );
        tMeshData->mVerbose          = false;
        std::string tMeshOutputFile2 = "./xtk_exo/unit_enrichment_1_background.e";
        tMeshData->create_output_mesh( tMeshOutputFile2 );

        auto tField = std::make_shared< moris::ge::Mesh_Field >( tMeshData, tLSFName );
        Cell< std::shared_ptr< ge::Geometry > > tGeometry = { std::make_shared< ge::Level_Set_Geometry >( tField ) };

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

        /*
         * Setup XTK Model and tell it how to cut
         */
        size_t                          tModelDimension       = 3;
        Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
        Model                           tXTKModel( tModelDimension, tMeshData, &tGeometryEngine );
        tXTKModel.mVerbose = false;
        /*
         * Decompose
         */
        tXTKModel.decompose( tDecompositionMethods );

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE );

        Enrichment const& tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
        Cell< std::string > tEnrichmentFieldNames;
        if ( tOutputEnrichmentFields )
        {
            tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        }

        delete tMeshData;
    }
}

TEST_CASE( "8 Element 10 enrichment Levels", "[ENRICH_10_EL_CLUSTER]" )
{
    if ( par_size() == 1 )
    {
        bool tOutputEnrichmentFields = true;

        // Generate mesh from string
        std::string tMeshFileName = "generated:2x2x2";

        // Add level set field to add onto file

        // Specify field parameters
        moris::mtk::Scalar_Field_Info< DDRMat > tLSF;
        std::string                             tLSFName = "lsf1";
        tLSF.set_field_name( tLSFName );
        tLSF.set_field_entity_rank( mtk::EntityRank::NODE );

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input( &tLSF, tFieldsInfo );

        // add to mesh data input container
        moris::mtk::MtkMeshData tSuppMeshData;
        tSuppMeshData.FieldsInfo = &tFieldsInfo;

        // Create mesh with supplementary data
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tSuppMeshData );

        xtk::size_t tNumNodes = tMeshData->get_num_entities( mtk::EntityRank::NODE );

        moris::Matrix< moris::DDRMat > tLevelsetVal( tNumNodes, 1, -1.2 );

        moris_id tIndexOfNodeId1  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId3  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId5  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 5, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId7  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 7, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId9  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 9, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId11 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 11, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId13 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 13, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId15 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 15, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId17 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 17, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId19 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 19, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId21 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 21, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId23 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 23, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId25 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 25, mtk::EntityRank::NODE );
        moris_id tIndexOfNodeId27 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 27, mtk::EntityRank::NODE );

        // Bottom face
        tLevelsetVal( tIndexOfNodeId1 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId3 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId7 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId9 ) = 1.1;

        // Top Face
        tLevelsetVal( tIndexOfNodeId19 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId21 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId25 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId27 ) = 1.1;

        tLevelsetVal( tIndexOfNodeId5 )  = 1.1;
        tLevelsetVal( tIndexOfNodeId11 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId17 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId23 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId15 ) = 1.1;
        tLevelsetVal( tIndexOfNodeId13 ) = 1.1;

        tMeshData->add_mesh_field_real_scalar_data_loc_inds( tLSFName, mtk::EntityRank::NODE, tLevelsetVal );
        tMeshData->mVerbose          = false;
        std::string tMeshOutputFile2 = "./xtk_exo/enrichment_test_10_cluster_background.e";
        tMeshData->create_output_mesh( tMeshOutputFile2 );

        auto tField = std::make_shared< moris::ge::Mesh_Field >( tMeshData, tLSFName );
        Cell< std::shared_ptr< ge::Geometry > > tGeometry = { std::make_shared< ge::Level_Set_Geometry >( tField ) };

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

        /*
         * Setup XTK Model and tell it how to cut
         */
        size_t                          tModelDimension       = 3;
        Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
        Model                           tXTKModel( tModelDimension, tMeshData, &tGeometryEngine );
        tXTKModel.mVerbose = false;

        /*
         * Decompose
         */
        tXTKModel.decompose( tDecompositionMethods );

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE );

        Enrichment const& tEnrichment = tXTKModel.get_basis_enrichment();

        // declare cell enrichment fields in output mesh
        Cell< std::string > tEnrichmentFieldNames;
        if ( tOutputEnrichmentFields )
        {
            tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        }

        Output_Options tOutputOptions;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        Enriched_Integration_Mesh&   tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        tEnrIntegMesh.create_dbl_sided_interface_sets( { 1 }, { 0 } );

        moris::Cell< mtk::Cluster const* > tDoubleSideCluster = tEnrIntegMesh.get_double_side_set_cluster( 0 );
        moris::real                        tGoldVolume        = 2;
        moris::real                        tGoldSurface       = 0.6862003781;

        for ( moris::uint i = 0; i < tDoubleSideCluster.size(); i++ )
        {
            moris::real tLeaderPrimaryVolume = tDoubleSideCluster( i )->compute_cluster_cell_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER );
            moris::real tFollowerPrimaryVolume  = tDoubleSideCluster( i )->compute_cluster_cell_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::FOLLOWER );
            moris::real tLeaderVoidVolume    = tDoubleSideCluster( i )->compute_cluster_cell_measure( mtk::Primary_Void::VOID, mtk::Leader_Follower::LEADER );
            moris::real tFollowerVoidVolume     = tDoubleSideCluster( i )->compute_cluster_cell_measure( mtk::Primary_Void::VOID, mtk::Leader_Follower::FOLLOWER );
            CHECK( std::abs( tLeaderPrimaryVolume + tFollowerPrimaryVolume + tLeaderVoidVolume + tFollowerVoidVolume - tGoldVolume ) < 1e-8 );

            moris::real tLeaderSurfaceArea = tDoubleSideCluster( i )->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER );
            CHECK( std::abs( tGoldSurface - tLeaderSurfaceArea ) < 1e-8 );

            moris::real tFollowerSurfaceArea = tDoubleSideCluster( i )->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::FOLLOWER );
            CHECK( std::abs( tGoldSurface - tFollowerSurfaceArea ) < 1e-8 );
        }

        // Enriched_Interpolation_Mesh& tEnrIpMesh    = tXTKModel.get_enriched_interp_mesh();
        // std::cout < "Check Mesh" << std::endl;
        // mtk::Mesh_Checker tMeshChecker( 0, &tEnrIpMesh, &tEnrIntegMesh );
        // tMeshChecker.perform();
        // tMeshChecker.print_diagnostics();
        // std::cout << "Basis indexing:" << tMeshChecker.bool_to_string( tMeshChecker.verify_basis_indexing( &tEnrIpMesh ) ) << std::endl;
        // std::cout << "Side Clustering:" << tMeshChecker.bool_to_string( tMeshChecker.verify_double_side_sets( &tEnrIntegMesh ) ) << std::endl;
        delete tMeshData;
    }
}

}// namespace xtk

