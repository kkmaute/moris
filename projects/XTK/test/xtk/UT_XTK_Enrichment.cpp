/*
 * cl_XTK_Enrichment.cpp
 *
 *  Created on: Dec 1, 2017
 *      Author: doble
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

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_GEN_Mesh_Field_Geometry.hpp"
#include "cl_GEN_Plane.hpp"

namespace xtk
{



TEST_CASE("Enrichment Example 1","[ENRICH_1]")
{
    if(par_size() == 1 || par_size() == 1)
    {
        bool tOutputEnrichmentFields = true;

        // Generate mesh from string
        std::string tMeshFileName     = "generated:1x1x3|sideset:XzZ";

        // Add level set field to add onto file

        // Specify field parameters
        moris::mtk::Scalar_Field_Info<DDRMat> tLSF;
        std::string tLSFName = "lsf1";
        tLSF.set_field_name(tLSFName);
        tLSF.set_field_entity_rank(moris::EntityRank::NODE);

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input(&tLSF,tFieldsInfo);

        // add to mesh data input container
        moris::mtk::MtkMeshData tSuppMeshData;
        tSuppMeshData.FieldsInfo = &tFieldsInfo;

        // Create mesh with supplementary data
        moris::mtk::Interpolation_Mesh* tMeshData   = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tSuppMeshData );


        xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);

        moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.3);

        moris_id tIndexOfNodeId6 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 6,EntityRank::NODE);
        moris_id tIndexOfNodeId3 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3,EntityRank::NODE);

        // Bottom face
        tLevelsetVal(tIndexOfNodeId3) = 1;
        tLevelsetVal(tIndexOfNodeId6) = 1;



        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
        tMeshData->mVerbose = false;
        std::string tMeshOutputFile2 = "./xtk_exo/unit_enrichment_1_background.e";
        tMeshData->create_output_mesh(tMeshOutputFile2);

        moris::Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
        tGeometry(0) = std::make_shared<moris::ge::Mesh_Field_Geometry>(tMeshData, tLSFName);

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

        /*
         * Setup XTK Model and tell it how to cut
         */
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
        tXTKModel.mVerbose  =  false;
        /*
         * Decompose
         */
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(EntityRank::NODE);


        Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();

        // Declare the fields related to enrichment strategy in output options
        Cell<std::string> tEnrichmentFieldNames;
        if(tOutputEnrichmentFields)
        {
            tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        }

        // TODO: run some FEM Temperature problem perturbing an enrichment level and checking whether other disconnected subdomains are heated up.

        // setup output mesh options with cell enrichment fields
        Output_Options tOutputOptions;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;


        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);

        std::string tMeshOutputFile = "./xtk_exo/unit_enrichment_1.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tCutMeshData;
    }

}

TEST_CASE("8 Element 10 enrichment Levels","[ENRICH_10_EL_CLUSTER]")
{
    if(par_size() == 1)
    {
        bool tOutputEnrichmentFields = true;

        // Generate mesh from string
        std::string tMeshFileName     = "generated:2x2x2";

        // Add level set field to add onto file

        // Specify field parameters
        moris::mtk::Scalar_Field_Info<DDRMat> tLSF;
        std::string tLSFName = "lsf1";
        tLSF.set_field_name(tLSFName);
        tLSF.set_field_entity_rank(moris::EntityRank::NODE);

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input(&tLSF,tFieldsInfo);

        // add to mesh data input container
        moris::mtk::MtkMeshData tSuppMeshData;
        tSuppMeshData.FieldsInfo = &tFieldsInfo;

        // Create mesh with supplementary data
        moris::mtk::Interpolation_Mesh* tMeshData   = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tSuppMeshData );


        xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);

        moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.2);

        moris_id tIndexOfNodeId1  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1,EntityRank::NODE);
        moris_id tIndexOfNodeId3  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3,EntityRank::NODE);
        moris_id tIndexOfNodeId5  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 5,EntityRank::NODE);
        moris_id tIndexOfNodeId7  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 7,EntityRank::NODE);
        moris_id tIndexOfNodeId9  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 9,EntityRank::NODE);
        moris_id tIndexOfNodeId11 = tMeshData->get_loc_entity_ind_from_entity_glb_id(11,EntityRank::NODE);
        moris_id tIndexOfNodeId13 = tMeshData->get_loc_entity_ind_from_entity_glb_id(13,EntityRank::NODE);
        moris_id tIndexOfNodeId15 = tMeshData->get_loc_entity_ind_from_entity_glb_id(15,EntityRank::NODE);
        moris_id tIndexOfNodeId17 = tMeshData->get_loc_entity_ind_from_entity_glb_id(17,EntityRank::NODE);
        moris_id tIndexOfNodeId19 = tMeshData->get_loc_entity_ind_from_entity_glb_id(19,EntityRank::NODE);
        moris_id tIndexOfNodeId21 = tMeshData->get_loc_entity_ind_from_entity_glb_id(21,EntityRank::NODE);
        moris_id tIndexOfNodeId23 = tMeshData->get_loc_entity_ind_from_entity_glb_id(23,EntityRank::NODE);
        moris_id tIndexOfNodeId25 = tMeshData->get_loc_entity_ind_from_entity_glb_id(25,EntityRank::NODE);
        moris_id tIndexOfNodeId27 = tMeshData->get_loc_entity_ind_from_entity_glb_id(27,EntityRank::NODE);


        // Bottom face
        tLevelsetVal(tIndexOfNodeId1) = 1.1;
        tLevelsetVal(tIndexOfNodeId3) = 1.1;
        tLevelsetVal(tIndexOfNodeId7) = 1.1;
        tLevelsetVal(tIndexOfNodeId9) = 1.1;

        // Top Face
        tLevelsetVal(tIndexOfNodeId19) = 1.1;
        tLevelsetVal(tIndexOfNodeId21) = 1.1;
        tLevelsetVal(tIndexOfNodeId25) = 1.1;
        tLevelsetVal(tIndexOfNodeId27) = 1.1;

        tLevelsetVal(tIndexOfNodeId5) = 1.1;
        tLevelsetVal(tIndexOfNodeId11) = 1.1;
        tLevelsetVal(tIndexOfNodeId17) = 1.1;
        tLevelsetVal(tIndexOfNodeId23) = 1.1;
        tLevelsetVal(tIndexOfNodeId15) = 1.1;
        tLevelsetVal(tIndexOfNodeId13) = 1.1;


        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
        tMeshData->mVerbose = false;
        std::string tMeshOutputFile2 = "./xtk_exo/enrichment_test_10_cluster_background.e";
        tMeshData->create_output_mesh(tMeshOutputFile2);

        moris::Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
        tGeometry(0) = std::make_shared<moris::ge::Mesh_Field_Geometry>(tMeshData, tLSFName);

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

        /*
         * Setup XTK Model and tell it how to cut
         */
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
        tXTKModel.mVerbose  = false;

        /*
         * Decompose
         */
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(EntityRank::NODE);

        Enrichment const & tEnrichment = tXTKModel.get_basis_enrichment();


        // declare cell enrichment fields in output mesh
        Cell<std::string> tEnrichmentFieldNames;
        if(tOutputEnrichmentFields)
        {
            tEnrichmentFieldNames = tEnrichment.get_cell_enrichment_field_names();
        }

        Output_Options tOutputOptions;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);


        std::string tMeshOutputFile = "./xtk_exo/enrichment_test_10_cluster.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);


        Enriched_Integration_Mesh & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

        tEnrIntegMesh.create_dbl_sided_interface_set(1,0);

        moris::Cell<mtk::Cluster const*> tDoubleSideCluster = tEnrIntegMesh.get_double_side_set_cluster(0);
        moris::real tGoldVolume = 2;
        moris::real tGoldSurface = 0.6862003781;

        for(moris::uint i = 0; i < tDoubleSideCluster.size(); i++)
        {
            moris::real tMasterPrimaryVolume = tDoubleSideCluster(i)->compute_cluster_cell_measure(mtk::Primary_Void::PRIMARY, mtk::Master_Slave::MASTER);
            moris::real tSlavePrimaryVolume = tDoubleSideCluster(i)->compute_cluster_cell_measure(mtk::Primary_Void::PRIMARY, mtk::Master_Slave::SLAVE);
            moris::real tMasterVoidVolume = tDoubleSideCluster(i)->compute_cluster_cell_measure(mtk::Primary_Void::VOID, mtk::Master_Slave::MASTER);
            moris::real tSlaveVoidVolume = tDoubleSideCluster(i)->compute_cluster_cell_measure(mtk::Primary_Void::VOID, mtk::Master_Slave::SLAVE);
            CHECK(std::abs(tMasterPrimaryVolume + tSlavePrimaryVolume + tMasterVoidVolume + tSlaveVoidVolume - tGoldVolume) < 1e-8);

            moris::real tMasterSurfaceArea = tDoubleSideCluster(i)->compute_cluster_cell_side_measure(mtk::Primary_Void::PRIMARY, mtk::Master_Slave::MASTER);
            CHECK(std::abs(tGoldSurface - tMasterSurfaceArea) < 1e-8 );

            moris::real tSlaveSurfaceArea = tDoubleSideCluster(i)->compute_cluster_cell_side_measure(mtk::Primary_Void::PRIMARY, mtk::Master_Slave::SLAVE);
            CHECK(std::abs(tGoldSurface - tSlaveSurfaceArea) < 1e-8 );
        }


        delete tMeshData;
        delete tCutMeshData;
    }
}


}
