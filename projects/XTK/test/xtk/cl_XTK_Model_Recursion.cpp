


/*
 * cl_XTK_Model_Recursion.hpp
 *
 *  Created on: Oct 20, 2017
 *      Author: ktdoble
 */

#include <memory>
#include <mpi.h>
#include "catch.hpp"
#include "ios/cl_Logger.hpp"

#include <cmath>

#include "xtk/cl_XTK_Child_Mesh.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geometry/cl_Plane.hpp"
// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/fn_verify_tet_topology.hpp"

// XTKL: Geometry  Include
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Phase_Table.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Face_Registry.hpp"

#include "xtk/fn_local_child_mesh_flood_fill.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "xtk/fn_generate_element_to_element.hpp"


namespace xtk
{
TEST_CASE("Phase Table","[Phase_Table]")
{

    moris::Matrix< Default_Matrix_Integer > tPhaseTableData (
            {{0,0},
             {0,1},
             {1,0},
             {1,1}});

    Cell<std::string> tPhaseNames = {"m0","m1","m2"};

    // TODO: figure out how to check a throw in this constructor
//    CHECK_THROWS(Phase_Table<size_t, Default_Matrix_Integer>(tPhaseTableData,tPhaseNames));

    tPhaseNames = {"m0","m1","m2","m3"};

    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (tPhaseTableData,tPhaseNames);
    moris::Matrix< Default_Matrix_Integer > tRow(0,0);

    size_t tIndex = 0;
    for(size_t iR = 0; iR<tPhaseTableData.n_rows(); iR++ )
    {
        tRow = tPhaseTableData.get_row(iR);

        tIndex = tPhaseTable.get_phase_index(tRow);

        CHECK(tIndex == iR);
    }

    (tRow)(0,0) = 2;
    CHECK_THROWS(tPhaseTable.get_phase_index(tRow));
    (tRow)(0,0) = 0;
    (tRow)(0,1) = 2;

    CHECK_THROWS(tPhaseTable.get_phase_index(tRow));


    // Check a 3 phase problem

    tPhaseTableData = moris::Matrix< Default_Matrix_Integer >(
                {{0,0,0},
                 {0,0,1},
                 {0,1,0},
                 {0,1,1},
                 {1,0,0},
                 {1,0,1},
                 {1,1,0},
                 {1,1,1}});


        // TODO: figure out how to check a throw in this constructor
    //    CHECK_THROWS(Phase_Table<size_t, Default_Matrix_Integer>(tPhaseTableData,tPhaseNames));

        tPhaseNames = {"m0","m1","m2","m3","m4","m5","m6","m7"};

        Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable2 (tPhaseTableData,tPhaseNames);

        //Check indices are correct
        tRow = moris::Matrix< Default_Matrix_Integer >(1,tPhaseTableData.n_cols());
        for(size_t iR = 0; iR<tPhaseTableData.n_rows(); iR++ )
        {
            tRow = tPhaseTableData.get_row(iR);

            tIndex = tPhaseTable2.get_phase_index(tRow);

            CHECK(tIndex == iR);
        }

        (tRow)(0,0) = 2;
        CHECK_THROWS(tPhaseTable2.get_phase_index(tRow));
        (tRow)(0,0) = 0;
        (tRow)(0,1) = 2;

        CHECK_THROWS(tPhaseTable2.get_phase_index(tRow));



}


TEST_CASE("Autogenerate Exponential Base 2 Table","[AUTO_PHASE_TABLE]")
{
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (3, Phase_Table_Structure::EXP_BASE_2);


    moris::Matrix< Default_Matrix_Integer > tExpectedPhaseTableData (
                {{0,0,0},
                 {0,0,1},
                 {0,1,0},
                 {0,1,1},
                 {1,0,0},
                 {1,0,1},
                 {1,1,0},
                 {1,1,1}});

    CHECK(xtk::equal_to(tExpectedPhaseTableData,tPhaseTable.get_phase_table_data()));
}
TEST_CASE("2 Nonintersecting geometries","[2_Phase],[NO_OVER]")
{

    moris::Matrix< Default_Matrix_Integer > tPhaseTableData (
            {{0,0},
             {0,1},
             {1,0},
             {1,1}});

    Cell<std::string> tPhaseNames = {"m0","m1","m2","m3"};

    // Geometry Engine Setup -----------------------
    // Using a 2 Levelset Spheres as the Geometry
    real tRadius1  = 2.70;
    real tXCenter1 = 4.51;
    real tYCenter1 = 4.51;
    real tZCenter1 = 4.51;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere1(tRadius1, tXCenter1, tYCenter1, tZCenter1);

    real tRadius2  = 2.8;
    real tXCenter2 = 4.51;
    real tYCenter2 = 4.51;
    real tZCenter2 = 5.51;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere2(tRadius2, tXCenter2, tYCenter2, tZCenter2);

//    real tRadius3  = 3.25;
//    real tXCenter3 = 10.0;
//    real tYCenter3 = 10.0;
//    real tZCenter3 = 10.0;
//    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere3(tRadius3, tXCenter3, tYCenter3, tZCenter3);


    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (2,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(
        {&tLevelsetSphere1,
         &tLevelsetSphere2},
        tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:10x10x10";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    Output_Options<size_t> tOutputOptions;

    // Write output mesh
    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData =
            tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/unit_recursive_no_intersect.e";
    tCutMeshData->write_output_mesh(tMeshOutputFile);



}


TEST_CASE("2 Intersecting Geometries","[2_Phase][OVER]")
{


    moris::Matrix< Default_Matrix_Integer > tPhaseTableData (
            {{0,0},
             {0,1},
             {1,0},
             {1,1}});

    Cell<std::string> tPhaseNames = {"m0","m1","m2","m3"};

    // Geometry Engine Setup -----------------------
    real tXc1 = 0.4;
    real tYc1 = 0.4;
    real tZc1 = 0.4;
    real tXn1 = 0.0;
    real tYn1 = 1.0;
    real tZn1 = 0.0;

    Plane<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tPlane1(tXc1,tYc1,tZc1,tXn1,tYn1,tZn1);


    real tXc2 = 0.55;
    real tYc2 = 0.55;
    real tZc2 = 0.55;
    real tXn2 = 0.0;
    real tYn2 = 0.0;
    real tZn2 = 1.0;

    Plane<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tPlane2(tXc2,tYc2,tZc2,tXn2,tYn2,tZn2);

    real tXc3 = 0.7;
    real tYc3 = 0.7;
    real tZc3 = 0.7;
    real tXn3 = 0.0;
    real tYn3 = 0.0;
    real tZn3 = 1.0;

    Plane<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tPlane3(tXc3,tYc3,tZc3,tXn3,tYn3,tZn3);


    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (3,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine({&tPlane1, &tPlane2,&tPlane3},tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:2x2x2";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string( tMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);


    // Verify that the tets created have correct topology
    Cut_Mesh<real, size_t, Default_Matrix_Real,Default_Matrix_Integer> & tCutMesh          = tXTKModel.get_cut_mesh();
    Child_Mesh_Test<real, size_t, Default_Matrix_Real,Default_Matrix_Integer> & tChildMesh = tCutMesh.get_child_mesh(0);

    bool tValidTopo = verify_tet4_topology(tChildMesh.get_element_to_node(),
                                           tChildMesh.get_element_to_edge(),
                                           tChildMesh.get_element_to_face(),
                                           tChildMesh.get_edge_to_node(),
                                           tChildMesh.get_face_to_node());

    CHECK(tValidTopo);


    Output_Options<size_t> tOutputOptions;
    tOutputOptions.mAddNodeSets = false;
    tOutputOptions.mAddSideSets = false;
    tOutputOptions.mInternalUseFlag = true;
    tOutputOptions.mRealNodeExternalFieldNames = {"lsf1","lsf2"};
    tOutputOptions.mIntElementExternalFieldNames= {"loc_ff","loc_ff2"};
    std::shared_ptr<mesh::Mesh_Data<real, size_t,Default_Matrix_Real,Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder,tOutputOptions);

//    size_t tNumNodes = tCutMeshData->get_num_entities(EntityRank::NODE);
//    size_t tNumElems = tCutMeshData->get_num_entities(EntityRank::ELEMENT);
//    Cell<real> tLSF1Vals(tNumNodes);
//    Cell<real> tLSF2Vals(tNumNodes);
//
//    for(size_t i = 0; i<tNumNodes; i++)
//    {
//        real tLSV = tXTKModel.get_geom_engine().get_entity_phase_val(i,0);
//
//        size_t tNodeId = tMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//
//        size_t tNodeIndex = tCutMeshData->get_loc_entity_index_from_entity_glb_id(tNodeId,EntityRank::NODE);
//
//        tLSF1Vals(tNodeIndex) = tLSV;
//    }
//    tCutMeshData->add_mesh_field_data_loc_indices(tOutputOptions.mRealNodeExternalFieldNames(0), EntityRank::NODE, tLSF1Vals);
//
//
//    for(size_t i = 0; i<tNumNodes; i++)
//    {
//        real tLSV = tXTKModel.get_geom_engine().get_entity_phase_val(i,1);
//
//        size_t tNodeId = tMeshData->get_glb_entity_id_from_entity_loc_index(i,EntityRank::NODE);
//
//        size_t tNodeIndex = tCutMeshData->get_loc_entity_index_from_entity_glb_id(tNodeId,EntityRank::NODE);
//
//        tLSF2Vals(tNodeIndex) = tLSV;
//    }
//    tCutMeshData->add_mesh_field_data_loc_indices(tOutputOptions.mRealNodeExternalFieldNames(1), EntityRank::NODE, tLSF2Vals);
//
//
//    std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tElementSubphase = local_child_mesh_flood_fill(tChildMesh, tMatrixFactory);
//    Cell<real> tSubphase(tNumElems);
//    moris::Matrix< Default_Matrix_Integer > const & tElemIds = tChildMesh.get_element_ids();
//    for(size_t i = 0; i<tNumElems; i++)
//    {
//        size_t tElemIndex = tCutMeshData->get_loc_entity_index_from_entity_glb_id(tElemIds(0,i),EntityRank::ELEMENT);
//        tSubphase(tElemIndex) = (real)(tElementSubphase)(0,i);
//    }
//    tCutMeshData->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(0), EntityRank::ELEMENT, tSubphase);
//
//    moris::Matrix< Default_Matrix_Integer > tElementToNode = tChildMesh.get_element_to_node();
//
//
//    size_t tNumElem = tElementToNode.n_rows();
//    size_t tNumFacePerElem = 4;
//    size_t tMax = std::numeric_limits<size_t>::max();
//    std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tFacAncInds (1,1,0);
//    std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tFacAncRank (1,1,3);
//    moris::Matrix< Default_Matrix_Integer > tActiveElements (1,tNumElem);
//    moris::Matrix< Default_Matrix_Integer > tElementPhase (1,tNumElem);
//
//    size_t t0 = 0;
//    std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tExpectedFaceToNode = construct_expected_face_to_node_tet4(t0,tElementToNode);
//    (tActiveElements)(0,0) = 0;
//    (tElementPhase)(0,0) = tChildMesh.get_element_phase_index( t0);
//    std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tFaceToElement ({{0},{0},{0},{0}});
//
//    Face_Registry<real, size_t, Default_Matrix_Real,Default_Matrix_Integer> tFaceRegistry(tNumElem,tNumFacePerElem,tExpectedFaceToNode,tFaceToElement,tFacAncInds,tFacAncRank);
//
//    for(size_t i = 1; i<tNumElem; i++)
//    {
//        std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tExpectedFaceToNode = construct_expected_face_to_node_tet4(i,tElementToNode);
//
//        std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tFaceInds = tFaceRegistry.get_face_indices(tExpectedFaceToNode,true);
//        print(tFaceInds,"tFaceInds");
//
//        (tActiveElements)(0,i) = i;
//        (tElementPhase)(0,i) = tChildMesh.get_element_phase_index( i);
//        for(size_t j = 0; j<tNumFacePerElem; j++)
//        {
//            if((tFaceInds)(0,j) == std::numeric_limits<size_t>::max())
//            {
//                (tFaceInds)(0,j) = tFaceRegistry.append_face(j, tExpectedFaceToNode);
//            }
//        }
//
//        print(tFaceInds,"tFaceInds end");
//        std::shared_ptr<Matrix_Base<size_t,Default_Matrix_Integer>> tElementInd (1,1);
//        (tElementInd)(0,0) = i;
//        tFaceRegistry.set_face_to_element(tElementInd,tFaceInds);
//    }
//
//    tFaceToElement = tFaceRegistry.get_face_to_element();
//
//    print(tFaceToElement,"tFaceToElement");
//
//    moris::Matrix< Default_Matrix_Integer > tElementToElement = generate_element_to_element(tFaceToElement,tNumElem,tNumFacePerElem,tMax);
//
//    moris::Matrix< Default_Matrix_Integer > tIncludedElementMarker (1,tNumElem);
//    tIncludedElementMarker->fill(1); // Mark all elements as included
//
//
//
//
//    // Run flood fill Algorithm
//    size_t tNumPhases = 4;
//     tElementSubphase = flood_fill( tElementToElement,
//                                    tElementPhase,
//                                    tActiveElements,
//                                    tIncludedElementMarker,
//                                    tNumPhases,
//                                    tMax,
//                                    true);
//
//
//     for(size_t i = 0; i<tNumElems; i++)
//     {
//         size_t tElemIndex = tCutMeshData->get_loc_entity_index_from_entity_glb_id(tElemIds(0,i),EntityRank::ELEMENT);
//         tSubphase(tElemIndex) = (real)(tElementSubphase)(0,i);
//     }
//
//     tCutMeshData->add_mesh_field_data_loc_indices(tOutputOptions.mIntElementExternalFieldNames(1), EntityRank::ELEMENT, tSubphase);



     std::string tPrefix = std::getenv("XTKOUTPUT");
     std::string tMeshOutputFile = tPrefix + "/unit_recursive_intersect.e";


     tCutMeshData->write_output_mesh(tMeshOutputFile,tOutputOptions.mRealNodeExternalFieldNames,{},tOutputOptions.mIntElementExternalFieldNames);



}
}
