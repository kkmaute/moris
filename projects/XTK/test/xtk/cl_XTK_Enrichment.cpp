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
#include "mesh/fn_verify_tet_topology.hpp"

// XTKL: Geometry  Include
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"


#include "geometry/cl_Discrete_Level_Set.hpp"
#include "geometry/cl_Plane.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"

namespace xtk
{

void
create_checkerboard_pattern(size_t const &     aNumX,
                            size_t const &     aNumY,
                            size_t const &     aNumZ,
                            std::string        aFieldName,
                            moris::mtk::Mesh & aMesh)
{
    // Get information about number of nodes and their coordinates
    // Split into two loops to avoid rewriting add_mesh_field_data function and to collect all field data first then apply to mesh
    size_t tNumNodes = aMesh.get_num_entities(moris::EntityRank::NODE);
    moris::Matrix< moris::DDRMat > tFieldData(tNumNodes,1);

    bool tOn = true;
    size_t tCount = 0;
    for(size_t iZ = 0; iZ<aNumZ+1; iZ++)
    {
        for(size_t iY = 0; iY<aNumY+1; iY++)
        {
            for(size_t iX = 0; iX<aNumX+1; iX++)
            {
                if(tOn)
                {
                    tFieldData(tCount) = 1.1;
                    tOn = false;
                }
                else
                {
                    tFieldData(tCount) = -0.5;
                    tOn = true;
                }
                tCount++;
            }
            tOn = false;
        }
        tOn = true;
    }

    aMesh.add_mesh_field_real_scalar_data_loc_inds(aFieldName, moris::EntityRank::NODE, tFieldData);

}



Cell<std::string>
declare_enrichment_fields_in_output_options(size_t aNumBasis)
{
    // declare  enrichment fields
    Cell<std::string> tEnrichmentFields(aNumBasis + 1);
    std::string tBaseEnrich = "subphase_";
    for(size_t i = 0; i<aNumBasis; i++)
    {
        tEnrichmentFields(i) = tBaseEnrich + std::to_string(i);
    }


    // Add local floodfill field to the output mesh
    std::string tLocalFFStr = "local_ff";
    tEnrichmentFields(aNumBasis) = tLocalFFStr;
    return tEnrichmentFields;
}

void
write_enrichment_data_to_fields(size_t               aNumBasis,
                                Cut_Mesh           & aCutMesh,
                                moris::mtk::Mesh   & aOutputMesh,
                                Enrichment   const & aEnrichment,
                                Cell<std::string>  & aEnrichmentFields)
{
    // Local subphas bins
    moris::Matrix<moris::DDRMat> tLocalSubphaseVal(aOutputMesh.get_num_entities(moris::EntityRank::ELEMENT),1);
    for(size_t i = 0; i<aCutMesh.get_num_child_meshes(); i++)
    {

        Child_Mesh & tChildMesh = aCutMesh.get_child_mesh(i);

        moris::Matrix< moris::IndexMat > const & tElementSubphases = tChildMesh.get_elemental_subphase_bin_membership();

        moris::Matrix< moris::IdMat > const & tChildElementIds = tChildMesh.get_element_ids();


        for(size_t j = 0; j<tChildElementIds.n_cols(); j++)
        {
            moris::moris_index tElementInd = aOutputMesh.get_loc_entity_ind_from_entity_glb_id(tChildElementIds(0,j),moris::EntityRank::ELEMENT);
            moris::moris_index tBulkPhaseInd = tChildMesh.get_element_phase_index(j);
            tLocalSubphaseVal(tElementInd) = (real)(tElementSubphases(j))+10;

        }
    }

    std::string tLocalFFStr = "local_ff";
    aOutputMesh.add_mesh_field_real_scalar_data_loc_inds(tLocalFFStr, moris::EntityRank::ELEMENT, tLocalSubphaseVal);

    // Enrichment values
    Cell<moris::Matrix< moris::IndexMat >> const & tElementIdsInBasis = aEnrichment.get_element_ids_in_basis_support();
    Cell<moris::Matrix< moris::IndexMat >> const & tElementEnrichmentInBasis = aEnrichment.get_element_enrichment_levels_in_basis_support();

    for(size_t i = 0; i<aNumBasis; i++)
    {
        moris::Matrix<moris::DDRMat> tEnrichmentLevels(aOutputMesh.get_num_entities(moris::EntityRank::ELEMENT),1);
        tEnrichmentLevels.fill(0);

        for(size_t j = 0; j<tElementIdsInBasis(i).n_cols(); j++)
        {
            size_t tElementId = (tElementIdsInBasis(i))(0,j);
            size_t tElementInd = aOutputMesh.get_loc_entity_ind_from_entity_glb_id(tElementId,moris::EntityRank::ELEMENT);
            tEnrichmentLevels(tElementInd) = (real)(((tElementEnrichmentInBasis(i)))(0,j));

        }

        aOutputMesh.add_mesh_field_real_scalar_data_loc_inds(aEnrichmentFields(i), moris::EntityRank::ELEMENT, tEnrichmentLevels);

    }
}

TEST_CASE("Enrichment Example 1","[ENRICH_1]")
        {

    // This problem has all background elements intersected
    //TODO: ADD declare field to STK create_mesh

//    bool tOutputEnrichmentFields = true;
//     // Load Mesh
//     xtk::size_t tNumX = 3;
//     xtk::size_t tNumY = 3;
//     xtk::size_t tNumZ = 1;
//
//     std::string tMeshFileName = get_generated_mesh_string(tNumX,tNumY,tNumZ);
//     Cell<std::string> tFieldNames = {"lsf1"};
//     moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );
//     size_t tNumNodes = tMeshData->get_num_nodes();
//     create_checkerboard_pattern(tNumX,tNumY,tNumZ,tFieldNames(0),*tMeshData);
//
//     Discrete_Level_Set tLevelSetMesh(tMeshData,tFieldNames);
//
//     Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
//     Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);
//
//     tGeometryEngine.mThresholdValue = 0.0;
//     tGeometryEngine.mComputeDxDp = false;
//
//     /*
//      * Setup XTK Model and tell it how to cut
//      */
//     size_t tModelDimension = 3;
//     Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//     Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
//     tXTKModel.mSameMesh = true;
//
//
//     /*
//      * Decompose
//      */
//     tXTKModel.decompose(tDecompositionMethods);
//
//
//     // Create output mesh
//     Output_Options<size_t> tOutputOptions;
//     tOutputOptions.mInternalUseFlag = true;
//     tOutputOptions.mAddPhaseField = true;
//
//
//     // Declare the fields related to enrichment strategy in output options
//     if(tOutputEnrichmentFields)
//     {
//         declare_enrichment_fields_in_output_options(tNumNodes, tOutputOptions.mIntElementExternalFieldNames);
//     }
//
////     moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);
//
//     // Perform enrichment
//     Enrichment tEnrichment(2,&tXTKModel.get_cut_mesh(),&tXTKModel.get_background_mesh());
//     tEnrichment.perform_enrichment();

//     Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();

//     if(tOutputEnrichmentFields)
//     {
//         write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tOutputOptions.mIntElementExternalFieldNames);
//     }



//
//     std::string tPrefix = std::getenv("XTKOUTPUT");
//     std::string tMeshOutputFile = tPrefix + "/enrichment_test_1.e";
//     tCutMeshData->write_output_mesh(tMeshOutputFile,{},{},tOutputOptions.mIntElementExternalFieldNames,{},{});

//     delete tMeshData;
        }


TEST_CASE("8 Element 10 enrichment Levels","[ENRICH_10_EL_CLUSTER]")
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
    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tSuppMeshData );


    xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);

    moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.0);

    // Bottom face
    tLevelsetVal(0) = 1;
    tLevelsetVal(2) = 1;
    tLevelsetVal(6) = 1;
    tLevelsetVal(8) = 1;

    // Top Face
    tLevelsetVal(18) = 1;
    tLevelsetVal(20) = 1;
    tLevelsetVal(24) = 1;
    tLevelsetVal(26) = 1;

    tLevelsetVal(4) = 1;
    tLevelsetVal(10) = 1;
    tLevelsetVal(16) = 1;
    tLevelsetVal(22) = 1;
    tLevelsetVal(14) = 1;
    tLevelsetVal(12) = 1;


    tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
    tMeshData->mVerbose = true;
    std::string tPrefix2 = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile2 = tPrefix2 + "/enrichment_test_10_cluster_background.e";
    tMeshData->create_output_mesh(tMeshOutputFile2);

    Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});

    Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);

    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp = false;

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mSameMesh = true;

    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);


      // Declare the fields related to enrichment strategy in output options
    Cell<std::string> tEnrichmentFieldNames;
    if(tOutputEnrichmentFields)
    {
        tEnrichmentFieldNames = declare_enrichment_fields_in_output_options(tNumNodes);
    }

    // Perform the enrichment
    Enrichment tEnrichment(2,&tXTKModel.get_cut_mesh(),&tXTKModel.get_background_mesh());
    tEnrichment.perform_enrichment();


    // Create output mesh
    Output_Options tOutputOptions;
    tOutputOptions.mInternalUseFlag = true;
    tOutputOptions.mAddPhaseField = true;
    tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    if(tOutputEnrichmentFields)
    {
        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();
        write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tEnrichmentFieldNames);
    }

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/enrichment_test_10_cluster.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);



    delete tMeshData;
    delete tCutMeshData;
}



TEST_CASE("Mixed Unintersected and Intersected Parent Elements","[Enrich_2]")
        {

    // This problem has a mix between intersected background elements and non-intersected background elements
    bool tOutputEnrichmentFields = true;

    // Generate mesh from string and then adding a level set field
    xtk::size_t tNumX = 2;
    xtk::size_t tNumY = 2;
    xtk::size_t tNumZ = 2;
    std::string tMeshFileName ="generated:2x2x2";
    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

    xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);


    Plane tPlane1( 0.0, 0.0, 0.1, 0.0, 0.0, 1.0);
    Plane tPlane2( 0.0, 0.0, 0.1, 1.0, 0.0, 0.0);
    Plane tPlane3( 0.0, 0.0, 1.4, 1.0, 0.0, 0.0);
    Plane tPlane4( 0.0, 0.0, 1.5, 0.0, 1.0, 0.0);

    Phase_Table tPhaseTable (4, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine({&tPlane1,&tPlane2,&tPlane3,&tPlane4},tPhaseTable);

    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp = false;

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mSameMesh = true;
    tXTKModel.mVerbose = true;


    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);


    // Create output mesh
    Output_Options tOutputOptions;
    tOutputOptions.mInternalUseFlag = true;
    tOutputOptions.mAddPhaseField = true;


    // Declare the fields related to enrichment strategy in output options
    if(tOutputEnrichmentFields)
    {
        tOutputOptions.mRealElementExternalFieldNames = declare_enrichment_fields_in_output_options(tNumNodes);
    }

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    // Perform enrichment
    Enrichment tEnrichment(8,&tXTKModel.get_cut_mesh(),&tXTKModel.get_background_mesh());
    tEnrichment.perform_enrichment();


    if(tOutputEnrichmentFields)
    {
        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();
        write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tOutputOptions.mRealElementExternalFieldNames);
    }

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/enrichment_test_2.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);

    delete tCutMeshData;
    delete tMeshData;
        }

}
