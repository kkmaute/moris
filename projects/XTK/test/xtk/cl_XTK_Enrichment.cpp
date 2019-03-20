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

std::string
get_generated_mesh_string(size_t const & aNumX,
                          size_t const & aNumY,
                          size_t const & aNumZ)
{
 return "generated:" + std::to_string(aNumX) + "x" + std::to_string(aNumY) + "x" + std::to_string(aNumZ);
}


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



void
declare_enrichment_fields_in_output_options(size_t aNumBasis,
                                            Cell<std::string> & aOutputElementIntFields)
{
    // declare  enrichment fields
    Cell<std::string> tEnrichmentFields(aNumBasis);
    std::string tBaseEnrich = "subphase_";
    for(size_t i = 0; i<aNumBasis; i++)
    {
        tEnrichmentFields(i) = tBaseEnrich + std::to_string(i);
    }

    aOutputElementIntFields = tEnrichmentFields;

    // Add local floodfill field to the output mesh
    std::string tLocalFFStr = "local_ff";
    aOutputElementIntFields.push_back(tLocalFFStr);
}

void
write_enrichment_data_to_fields(size_t aNumBasis,
                                Cut_Mesh           & aCutMesh,
                                moris::mtk::Mesh                                                           & aOutputMesh,
                                Enrichment   const & aEnrichment,
                                Cell<std::string> aEnrichmentFieldStrs)
{

    // Local subphas bins
    moris::Matrix<moris::DDRMat> tLocalSubphaseVal(aOutputMesh.get_num_entities(moris::EntityRank::ELEMENT),0);
    for(size_t i = 0; i<aCutMesh.get_num_simple_meshes(); i++)
    {

        Child_Mesh_Test & tChildMesh = aCutMesh.get_child_mesh(i);

        moris::Matrix< moris::IndexMat > const & tElementSubphases = tChildMesh.get_elemental_subphase_bin_membership();

        moris::Matrix< moris::IdMat > const & tChildElementIds = tChildMesh.get_element_ids();


        for(size_t j = 0; j<tChildElementIds.n_cols(); j++)
        {
            moris::moris_index tElementInd = aOutputMesh.get_loc_entity_ind_from_entity_glb_id(tChildElementIds(0,j),moris::EntityRank::ELEMENT);
            moris::moris_index tBulkPhaseInd = tChildMesh.get_element_phase_index(j);
            tLocalSubphaseVal(tElementInd) = (real)(tElementSubphases(0,j));

        }
    }
    std::string tLocalFFStr = "local_ff";
    aOutputMesh.add_mesh_field_real_scalar_data_loc_inds(tLocalFFStr, moris::EntityRank::ELEMENT, tLocalSubphaseVal);



    // Enrichment values
    Cell<moris::Matrix< moris::IndexMat >> const & tElementIdsInBasis = aEnrichment.get_element_ids_in_basis_support();
    Cell<moris::Matrix< moris::IndexMat >> const & tElementEnrichmentInBasis = aEnrichment.get_element_enrichment_levels_in_basis_support();

    for(size_t i = 0; i<aNumBasis; i++)
    {
        moris::Matrix<moris::DDRMat> tEnrichmentLevels(aOutputMesh.get_num_entities(moris::EntityRank::ELEMENT),10);

        for(size_t j = 0; j<tElementIdsInBasis(i).n_cols(); j++)
        {
            size_t tElementId = (tElementIdsInBasis(i))(0,j);
            size_t tElementInd = aOutputMesh.get_loc_entity_ind_from_entity_glb_id(tElementId,moris::EntityRank::ELEMENT);
            tEnrichmentLevels(tElementInd) = (real)(((tElementEnrichmentInBasis(i)))(0,j));

        }

        aOutputMesh.add_mesh_field_real_scalar_data_loc_inds(aEnrichmentFieldStrs(i), moris::EntityRank::ELEMENT, tEnrichmentLevels);

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
//     Enrichment tEnrichment(2,&tXTKModel.get_cut_mesh(),&tXTKModel.get_xtk_mesh());
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

    //TODO: ADD declare field to STK create_mesh
    // This problem has a mix between intersected background elements and non-intersected background elements
//    bool tOutputEnrichmentFields = true;
//
//    // Generate mesh from string and then adding a level set field
//    xtk::size_t tNumX = 2;
//    xtk::size_t tNumY = 2;
//    xtk::size_t tNumZ = 2;
//    std::string tMeshFileName     = get_generated_mesh_string(tNumX,tNumY,tNumZ);
//    Cell<std::string> tFieldNames = {"lsf1"};
//    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );
//
//
//    xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);
//
//    moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1);
//
//    // Bottom face
//    tLevelsetVal(0) = 1;
//    tLevelsetVal(2) = 1;
//    tLevelsetVal(6) = 1;
//    tLevelsetVal(8) = 1;
//
//    // Center Node
//    tLevelsetVal(13) = 1;
//
//    // Top Face
//    tLevelsetVal(18) = 1;
//    tLevelsetVal(20) = 1;
//    tLevelsetVal(24) = 1;
//    tLevelsetVal(26) = 1;
//
//
//    tMeshData->add_mesh_field_real_scalar_data_loc_inds(tFieldNames(0), moris::EntityRank::NODE, tLevelsetVal);
//
//    Discrete_Level_Set tLevelSetMesh(tMeshData,tFieldNames);
//
//    Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
//    Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);
//
//    tGeometryEngine.mThresholdValue = 0.0;
//    tGeometryEngine.mComputeDxDp = false;
//
//    /*
//     * Setup XTK Model and tell it how to cut
//     */
//    size_t tModelDimension = 3;
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
//    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
//    tXTKModel.mSameMesh = true;
//
//
//    /*
//     * Decompose
//     */
//    tXTKModel.decompose(tDecompositionMethods);
//
//
//    // Create output mesh
//    Output_Options<size_t> tOutputOptions;
//    tOutputOptions.mInternalUseFlag = true;
//    tOutputOptions.mAddPhaseField = true;
//
//
////    // Declare the fields related to enrichment strategy in output options
////    if(tOutputEnrichmentFields)
////    {
////        declare_enrichment_fields_in_output_options(tNumNodes, tOutputOptions.mIntElementExternalFieldNames);
////    }
//
//    // Perform the enrichment
//    Enrichment tEnrichment(2,&tXTKModel.get_cut_mesh(),&tXTKModel.get_xtk_mesh());
//    tEnrichment.perform_enrichment();
//
////    if(tOutputEnrichmentFields)
////    {
////        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();
////        write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tOutputOptions.mIntElementExternalFieldNames);
////    }
//
//
//
//
////    std::string tPrefix = std::getenv("XTKOUTPUT");
////    std::string tMeshOutputFile = tPrefix + "/enrichment_test_10_cluster.e";
////    tCutMeshData->write_output_mesh(tMeshOutputFile,{},{},tOutputOptions.mIntElementExternalFieldNames,{},{});
//
//    delete tMeshData;
////    delete tCutMeshData;
        }



TEST_CASE("Mixed Unintersected and Intersected Parent Elements","[Enrich_2]")
        {

    // This problem has a mix between intersected background elements and non-intersected background elements
    bool tOutputEnrichmentFields = true;

    // Generate mesh from string and then adding a level set field
    xtk::size_t tNumX = 2;
    xtk::size_t tNumY = 2;
    xtk::size_t tNumZ = 2;
    std::string tMeshFileName = get_generated_mesh_string(tNumX,tNumY,tNumZ);
    Cell<std::string> tFieldNames = {"lsf"};
    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

    xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);


    Plane tPlane1( 0.0, 0.0, 0.1, 0.0, 0.0, 1.0);
    Plane tPlane2( 0.0, 0.0, 0.7, 0.0, 0.0, 1.0);



    Phase_Table tPhaseTable (2, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine({&tPlane1,&tPlane2},tPhaseTable);

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


    // Create output mesh
    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mInternalUseFlag = true;
    tOutputOptions.mAddPhaseField = true;


    // Declare the fields related to enrichment strategy in output options
    if(tOutputEnrichmentFields)
    {
        declare_enrichment_fields_in_output_options(tNumNodes, tOutputOptions.mIntElementExternalFieldNames);
    }
//    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    // Perform enrichment
    Enrichment tEnrichment(8,&tXTKModel.get_cut_mesh(),&tXTKModel.get_xtk_mesh());
    tEnrichment.perform_enrichment();


//    if(tOutputEnrichmentFields)
//
//    {
//        Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();
//        write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tOutputOptions.mIntElementExternalFieldNames);
//    }
//
//
//
//
//    std::string tPrefix = std::getenv("XTKOUTPUT");
//    std::string tMeshOutputFile = tPrefix + "/enrichment_test_2.e";
//    tCutMeshData->write_output_mesh(tMeshOutputFile,{},{},tOutputOptions.mIntElementExternalFieldNames,{},{});
        }

}
