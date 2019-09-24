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


#include "geometry/cl_Discrete_Level_Set.hpp"
#include "geometry/cl_Plane.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

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

        Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});


        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable,2);

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 2;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(EntityRank::NODE);

        // get the enriched interpolation mesh
        Enriched_Interpolation_Mesh const & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();

        Cell<Interpolation_Cell_Unzipped> const & tCells = tEnrInterpMesh.get_enriched_interpolation_cells();


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

        tXTKModel.print_subphase_neighborhood();

        std::string tMeshOutputFile = "./xtk_exo/enrich_2_element_ig.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tMeshData;
        delete tCutMeshData;
    }

}
}
