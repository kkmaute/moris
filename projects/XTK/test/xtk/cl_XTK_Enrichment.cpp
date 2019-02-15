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


#include "geometry/cl_Discrete_Level_Set.hpp"
#include "geometry/cl_Plane.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"

namespace xtk
{

class Vertex_Enrichment
{
public:
    Vertex_Enrichment(){}
private:
    moris::moris_index               mNodeIndex;
    moris::Matrix< moris::IndexMat > mBasisIndices;
    moris::Matrix< moris::IndexMat > mBasisEnrichmentLevel;
    moris::Matrix< moris::DDRMat >   mWeights;
};






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


TEST_CASE("Enrichment Example 1","[ENRICH_1]")
{
    if(par_size() == 1 || par_size() == 2)
    {
        bool tOutputEnrichmentFields = false;

        // Generate mesh from string
        std::string tMeshFileName     = "generated:1x1x1";

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

        moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.3);

        moris_id tIndexOfNodeId7  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 7,EntityRank::NODE);
        moris_id tIndexOfNodeId10 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1,EntityRank::NODE);

        // Bottom face
        tLevelsetVal(tIndexOfNodeId10) = 1;
        tLevelsetVal(tIndexOfNodeId7)  = 1;



        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
        tMeshData->mVerbose = true;
        std::string tPrefix2 = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile2 = tPrefix2 + "/unit_enrichment_1_background.e";
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
        tXTKModel.mVerbose = true;
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

        tXTKModel.unzip_interface();

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment();

        //
        Enrichment const & tEnrichment  = tXTKModel.get_basis_enrichment();

        Cell<moris::Matrix< moris::IdMat >> const & tElementIdsInBasisSupport = tEnrichment.get_element_inds_in_basis_support();
        Cell<moris::Matrix< moris::IdMat >> const & tELs = tEnrichment.get_element_enrichment_levels_in_basis_support();


        std::cout<<"tElementIdsInBasisSupport size"<<tElementIdsInBasisSupport.size()<<std::endl;
        std::cout<<"tELs size " <<tELs.size()<<std::endl;

        // Create output mesh
        Output_Options tOutputOptions;
        tOutputOptions.mInternalUseFlag = true;
        tOutputOptions.mAddPhaseField = true;
        tOutputOptions.mRealElementExternalFieldNames = tEnrichmentFieldNames;

        std::string tOwnerFieldName = "par_owner";
        tOutputOptions.mRealElementExternalFieldNames.push_back(tOwnerFieldName);

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        if(tOutputEnrichmentFields)
        {
            Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();
//            write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tEnrichmentFieldNames);
            write_element_ownership_as_field(tOwnerFieldName,
                                             tXTKModel.get_background_mesh(),
                                             tXTKModel.get_cut_mesh(),
                                             *tCutMeshData);

        }

        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/unit_enrichment_1.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);



        delete tMeshData;
        delete tCutMeshData;
    }

}

TEST_CASE("8 Element 10 enrichment Levels","[ENRICH_10_EL_CLUSTER]")
{
    if(par_size() == 1 || par_size() == 2)
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
        tLevelsetVal(tIndexOfNodeId1) = 1;
        tLevelsetVal(tIndexOfNodeId3) = 1;
        tLevelsetVal(tIndexOfNodeId7) = 1;
        tLevelsetVal(tIndexOfNodeId9) = 1;

        // Top Face
        tLevelsetVal(tIndexOfNodeId19) = 1;
        tLevelsetVal(tIndexOfNodeId21) = 1;
        tLevelsetVal(tIndexOfNodeId25) = 1;
        tLevelsetVal(tIndexOfNodeId27) = 1;

        tLevelsetVal(tIndexOfNodeId5) = 1;
        tLevelsetVal(tIndexOfNodeId11) = 1;
        tLevelsetVal(tIndexOfNodeId17) = 1;
        tLevelsetVal(tIndexOfNodeId23) = 1;
        tLevelsetVal(tIndexOfNodeId15) = 1;
        tLevelsetVal(tIndexOfNodeId13) = 1;


        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
        tMeshData->mVerbose = true;
        std::string tPrefix2 = std::getenv("MORISOUTPUT");
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
        tXTKModel.unzip_interface();


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

        std::string tOwnerFieldName = "par_owner";
        tOutputOptions.mRealElementExternalFieldNames.push_back(tOwnerFieldName);

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        if(tOutputEnrichmentFields)
        {
            Cut_Mesh & tCutMesh = tXTKModel.get_cut_mesh();

            // TODO:REWRITE THIS FUNCTION
//            write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tEnrichmentFieldNames);
            write_element_ownership_as_field(tOwnerFieldName,
                                             tXTKModel.get_background_mesh(),
                                             tXTKModel.get_cut_mesh(),
                                             *tCutMeshData);

        }

        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/enrichment_test_10_cluster.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);



        delete tMeshData;
        delete tCutMeshData;
    }
}



TEST_CASE("Mixed Unintersected and Intersected Parent Elements","[Enrich_2]")
{
    if(par_size() == 1 || par_size() == 2)
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
            // TODO:REWRITE THIS FUNCTION
//            write_enrichment_data_to_fields(tNumNodes,tCutMesh,*tCutMeshData,tEnrichment,tOutputOptions.mRealElementExternalFieldNames);
        }

        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/enrichment_test_2.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tCutMeshData;
        delete tMeshData;
    }
}

}
