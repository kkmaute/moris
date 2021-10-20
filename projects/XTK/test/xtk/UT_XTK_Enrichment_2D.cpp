/*
 * UT_XTK_Enrichment_2D.cpp
 *
 *  Created on: Sep 13, 2019
 *      Author: doble
 */

#include <memory>
#include <mpi.h>

#include "catch.hpp"

// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Vertex_Interpolation.hpp"
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
#include "fn_all_true.hpp"
#include "fn_sort.hpp"
#include "op_equal_equal.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_XTK_Cell_Cluster.hpp"

#include "cl_GEN_Mesh_Field_Geometry.hpp"
#include "cl_GEN_Plane.hpp"

namespace xtk
{



TEST_CASE("2 Element Enrichment 2D","[ENRICH_1E_2D]")
{
    if(par_size() == 1 || par_size() == 1)
    {
        bool tOutputEnrichmentFields = true;
        // Create Mesh ---------------------------------
        // Generate data for test
        uint aNumDim = 2;
        Matrix< DDRMat >  aCoords(6,2);
        aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
        aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
        aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
        aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
        aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
        aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
        Matrix< IdMat >     aElemConn( 2, 4 );

        // 0D to 3D connectivity (node to element)
        aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
        aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

        Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

        // Create MORIS mesh using MTK database
        moris::mtk::MtkMeshData aMeshData;
        aMeshData.CreateAllEdgesAndFaces = true;
        aMeshData.SpatialDim = &aNumDim;
        aMeshData.ElemConn(0)= &aElemConn;
        aMeshData.NodeCoords = &aCoords;
        aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;
        moris::mtk::Scalar_Field_Info<DDRMat> tLSF;
        std::string tLSFName = "lsf1";
        tLSF.set_field_name(tLSFName);
        tLSF.set_field_entity_rank(moris::EntityRank::NODE);

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input(&tLSF,tFieldsInfo);
        aMeshData.FieldsInfo = &tFieldsInfo;

        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, aMeshData );

        xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);
        moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.3);

        moris_id tIndexOfNodeId1 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1, EntityRank::NODE);
        moris_id tIndexOfNodeId3 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3, EntityRank::NODE);
        moris_id tIndexOfNodeId5 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 5, EntityRank::NODE);

        tLevelsetVal(tIndexOfNodeId1) = 1.0;
        tLevelsetVal(tIndexOfNodeId3) = 1.0;
        tLevelsetVal(tIndexOfNodeId5) = 1.0;

        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);

        Cell<std::shared_ptr<ge::Geometry>> tGeometry(1);
        tGeometry(0) = std::make_shared<moris::ge::Mesh_Field_Geometry>(tMeshData, tLSFName);

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 2;
        Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
        tXTKModel.mVerbose  =  false;

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
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

        // get the enriched meshes
        Enriched_Integration_Mesh   const & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh();
        Enriched_Interpolation_Mesh const & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();

        // Check the basis to enriched basis table
        moris::Cell<moris::Matrix<moris::IndexMat>> tGoldCoeffToEnrCoeff(6);
        tGoldCoeffToEnrCoeff(0) = {{0,1,2}};
        tGoldCoeffToEnrCoeff(1) = {{3,4,5,6}};
        tGoldCoeffToEnrCoeff(2) = {{7,8,9,10}};
        tGoldCoeffToEnrCoeff(3) = {{11,12,13}};
        tGoldCoeffToEnrCoeff(4) = {{14,15,16}};
        tGoldCoeffToEnrCoeff(5) = {{17,18,19}};

        for(moris::uint i = 0; i < tEnrInterpMesh.get_num_background_coefficients(0); i++)
        {
            Matrix<IndexMat> const & tEnrichedCoeffs = tEnrInterpMesh.get_enriched_coefficients_at_background_coefficient(0,(moris_index)i);
            CHECK(all_true(tEnrichedCoeffs == tGoldCoeffToEnrCoeff(i)));
        }

        // verify there is one interpolation cell per subphase
        moris::uint tExpectedNumSubphase = 6;
        Cell<Interpolation_Cell_Unzipped*> const & tCells = tEnrInterpMesh.get_enriched_interpolation_cells();
        CHECK(tCells.size() == tExpectedNumSubphase);

        // Expected primary cells in each cluster
        moris::Cell<moris::Matrix<moris::IndexMat>> tGoldPrimaryCellsInClusters(6);
        tGoldPrimaryCellsInClusters(0) = {{8,18}};
        tGoldPrimaryCellsInClusters(1) = {{9,11,12,13,15,16,17,19}};
        tGoldPrimaryCellsInClusters(2) = {{10,14}};
        tGoldPrimaryCellsInClusters(3) = {{20,22,25,26,27,29,30,31}};
        tGoldPrimaryCellsInClusters(4) = {{21,24}};
        tGoldPrimaryCellsInClusters(5) = {{23,28}};

        // Expected void cells in each cluster
        moris::Cell<moris::Matrix<moris::IndexMat>> tGoldVoidCellsInClusters(6);
        tGoldVoidCellsInClusters(0) = {{9,10,11,12,13,14,15,16,17,19}};
        tGoldVoidCellsInClusters(1) = {{8,10,14,18}};
        tGoldVoidCellsInClusters(2) = {{8,9,11,12,13,15,16,17,18,19}};
        tGoldVoidCellsInClusters(3) = {{21,23,24,28}};
        tGoldVoidCellsInClusters(4) = {{20,22,23,25,26,27,28,29,30,31}};
        tGoldVoidCellsInClusters(5) = {{20,21,22,24,25,26,27,29,30,31}};

//        // Expected interpolation vertices
//        moris::Cell<moris::Matrix<moris::IndexMat>> tGoldInterpCoeff(6);
//        tGoldInterpCoeff(0) = {{0,6,9,15}};
//        tGoldInterpCoeff(1) = {{1,7,10,16}};
//        tGoldInterpCoeff(2) = {{2,8,11,17}};
//        tGoldInterpCoeff(3) = {{3,18,21,12}};
//        tGoldInterpCoeff(4) = {{4,19,22,13}};
//        tGoldInterpCoeff(5) = {{ 5,20,23,14}};
//

        // Expected interpolation vertices
        moris::Cell<moris::Matrix<moris::IndexMat>> tGoldInterpCoeff(6);
        tGoldInterpCoeff(0) = {{0,6,7,11}};
        tGoldInterpCoeff(1) = {{1,3,8,12}};
        tGoldInterpCoeff(2) = {{2,5,9,13}};
        tGoldInterpCoeff(3) = {{3,14,17,8}};
        tGoldInterpCoeff(4) = {{4,15,18,10}};
        tGoldInterpCoeff(5) = {{ 5,16,19,9}};

        // iterate through cells and get cell clusters
        for(moris::uint i = 0; i < tCells.size() ; i++ )
        {
            moris::mtk::Cell const & tInterpCell = tEnrInterpMesh.get_mtk_cell((moris_index)i);

            xtk::Cell_Cluster const &           tXTKCellCluster = tEnrIntegMesh.get_xtk_cell_cluster(tInterpCell);
            Interpolation_Cell_Unzipped const * tXTKInterpCell  = tXTKCellCluster.get_xtk_interpolation_cell();

            // verify the interpolation cells are the same
            CHECK(tInterpCell.get_id() == tXTKInterpCell->get_id());

            // check the primary ids
            Matrix<IdMat> tPrimaryCellIds = tXTKCellCluster.get_primary_cell_ids_in_cluster();


            CHECK(all_true(tPrimaryCellIds == tGoldPrimaryCellsInClusters(i)));

            // check the void ids
            Matrix<IdMat> tVoidCellIds    = tXTKCellCluster.get_void_cell_ids_in_cluster();

            Matrix<IdMat> tSortedVoidCellIds;
            sort( tVoidCellIds, tSortedVoidCellIds, "ascend", 1 );

            CHECK(all_true(tSortedVoidCellIds == tGoldVoidCellsInClusters(i)));

            // check the vertices of the interpolation cell
            Cell<moris::mtk::Vertex *> tVertices = tXTKInterpCell->get_vertex_pointers();

            Matrix<IndexMat> tVertexInterpInds(1,4);
            for(moris::uint j = 0; j < tVertices.size(); j++)
            {
                mtk::Vertex_Interpolation* tInterp = tVertices(j)->get_interpolation(0);

                Matrix< IndexMat > tIndices  = tInterp->get_indices();

                // since this is a lagrange mesh there should be a 1-1 relationship between interp vert and coefficient
                // this does not necessarily imply they have the same index though
                CHECK(tIndices.numel() == 1);

                tVertexInterpInds(j) = tIndices(0);
            }

            CHECK(all_true(tVertexInterpInds == tGoldInterpCoeff(i)));
        }


        // outputting for visual check
        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        tEnrichment.write_cell_enrichment_to_fields(tEnrichmentFieldNames,tCutMeshData);

        std::string tMeshOutputFile = "./xtk_exo/enrich_2_element_ig.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tCutMeshData;
    }

}
}
