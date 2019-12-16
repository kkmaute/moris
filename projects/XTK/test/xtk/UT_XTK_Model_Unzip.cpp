
#include "catch.hpp"

// Xtk includes
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Background_Mesh.hpp"
#include "cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Interface_Element.hpp"

// Geometry
#include "cl_Sphere.hpp"
#include "cl_Plane.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

// Linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "op_plus.hpp"
#include "op_equal_equal.hpp"

using namespace moris;
using namespace xtk;


// Returns a matrix of node coordinates which have the interface coordinates shitfted by the interface element
// normals (this is not a fast code but helpful for debugging small problems)
moris::Matrix<moris::DDRMat>
extrude_interface_coordinates(xtk::Model  & aModel,
                              moris::real   aShiftDist)
{
    Background_Mesh & tBM = aModel.get_background_mesh();

    Cut_Mesh        & tCM = aModel.get_cut_mesh();

    // Allocate extruded coordinates
    moris::Matrix<moris::DDRMat> tCoords     = tBM.get_all_node_coordinates_loc_inds();
    moris::Matrix<moris::DDRMat> tCoordShift(tCoords.n_rows(),tCoords.n_cols(), 0.0);


    // Interface nodes
    moris::Matrix< moris::IndexMat > tIntNodes = tBM.get_interface_nodes_loc_inds(0);


    // Interface element objects
    moris::Cell<xtk::Interface_Element> & tIntElems = tCM.get_interface_elements();

    // interface element extracted
    moris::Matrix<moris::IndexMat> tInterfaceElementsInds = tCM.get_extracted_interface_elements_loc_inds();

    //
    moris::Matrix<moris::IndexMat> tNodeToInterfaceElem(tCoords.n_rows(),tIntElems.size(),MORIS_INDEX_MAX);
    moris::Matrix<moris::IndexMat> tNodeCount(tCoords.n_rows(),1,0);

    // loop over interface elements to construct a node to element index
    for(moris::uint i = 0; i <tIntElems.size(); i++)
    {

        for(moris::uint  j = 0 ; j< tInterfaceElementsInds.n_cols(); j++)
        {
            moris::moris_index tNodeIndex = tInterfaceElementsInds(i,j);
            moris::moris_index tCount     = tNodeCount(tNodeIndex);
            tNodeToInterfaceElem(tNodeIndex,tNodeCount(tNodeIndex)) = i;
            tNodeCount(tNodeIndex) = tNodeCount(tNodeIndex) + 1;
        }
    }



    // Loop over nodes and compute average normal of all elements connected to it
    for(moris::uint i = 0; i < tIntNodes.numel(); i++)
    {
        moris::moris_index tNodeIndex = tIntNodes(i);

        moris::Matrix<moris::F31RMat> tNodeNormal(3,1);


        moris::uint tCount = 0;
        for(moris::uint j = 0; j < tNodeToInterfaceElem.n_cols(); j++)
        {
            moris::moris_index tIntElemIndex = tNodeToInterfaceElem(tNodeIndex,j);

            if(tIntElemIndex == MORIS_INDEX_MAX)
            {
                break;
            }

            tNodeNormal = tNodeNormal + tIntElems(tIntElemIndex).get_outward_normal(0);
            tCount++;
        }

        tNodeNormal = tNodeNormal/tCount;

        tCoordShift.get_row(tNodeIndex) = aShiftDist*moris::trans(tNodeNormal);
    }


    return tCoords+tCoordShift;


}


TEST_CASE("XTK Model Unzip Interface","[unzip_xtk]")
{
    if(par_size() == 1)
    {
//        real tRadius  = 0.75;
        real tXCenter = 0.6;
        real tYCenter = 0.6;
        real tZCenter = 0.6;
//        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

        Plane tLevelsetSphere(tXCenter,tYCenter,tZCenter,0.0,1.0,0.0);
        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

        // Create Mesh --------------------------------------------------------------------
        std::string tMeshFileName = "generated:1x1x2";
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose  =  false;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        moris::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                      Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        // Test interface information prior to unzipping
        Background_Mesh const & tBM = tXTKModel.get_background_mesh();

        // get the interface nodes
        moris::Matrix< moris::IndexMat > tGoldInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54}});
        moris::Matrix< moris::IndexMat > tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);


        CHECK(all_true(tGoldInterfaceNodes == tInterfaceNodes));
        // Unzip interface
        tXTKModel.unzip_interface();

        // Now there should be twice as many nodes
        moris::Matrix< moris::IndexMat > tGoldUnzippedInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}});
        tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);

        // interface node pairs
        moris::Matrix<moris::IndexMat> tGoldInterfaceNodePairs({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54},
                                                                {55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}});

        moris::Matrix<moris::DDRMat> tCoords  = tBM.get_all_node_coordinates_loc_inds();
        moris::real tTol = 1e-8;
        // verify pairs have the same coordinates
        for(moris::uint i = 0; i < tGoldInterfaceNodePairs.n_cols(); i++)
        {
            CHECK( moris::norm(tCoords.get_row(tGoldInterfaceNodePairs(0,i)) - tCoords.get_row(tGoldInterfaceNodePairs(1,i)))<tTol);
        }



//
//        // Create a mesh with only the interface elements but extruded
//        // get the interface elements local node index element connectivity
//        xtk::Cut_Mesh & tCM     = tXTKModel.get_cut_mesh();
//        moris::mtk::MtkMeshData   tMeshDataInput;
//        moris::mtk::MtkSetsInfo   tMtkMeshSets;
//
//        moris::Matrix<moris::IdMat>  tLocalToGlobalNodeMap = tBM.get_local_to_global_map(EntityRank::NODE);
//        moris::Matrix<moris::DDRMat> tNodeCoordinates      = extrude_interface_coordinates(tXTKModel,1.5);
//        moris::Matrix<moris::IdMat> tNodeOwner(1,tBM.get_num_entities(EntityRank::NODE),moris::par_rank());
//
//        // Setup the interface element block set
//        moris::mtk::MtkBlockSetInfo tUnzippedInterfaceBlockSet;
//        moris::Matrix<moris::IndexMat> tInterfaceElementsInds = tCM.get_extracted_interface_elements_loc_inds();
//        tBM.convert_loc_entity_ind_to_glb_entity_ids(EntityRank::NODE,tInterfaceElementsInds);
//        moris::Matrix<moris::IndexMat> tInterfaceElementIds = tCM.get_interface_element_ids();
//        tUnzippedInterfaceBlockSet.mCellIdsInSet = &tInterfaceElementIds;
//        tUnzippedInterfaceBlockSet.mBlockSetName = "interface";
//        tUnzippedInterfaceBlockSet.mBlockSetTopo = CellTopology::PRISM6;
//        tMtkMeshSets.add_block_set(&tUnzippedInterfaceBlockSet);
//
//        moris::uint tSpatialDim = 3;
//
//        tMeshDataInput.ElemConn                = moris::Cell<moris::Matrix<IdMat>*>(1);
//        tMeshDataInput.LocaltoGlobalElemMap    = moris::Cell<moris::Matrix<IdMat>*>(1);
//        tMeshDataInput.CreateAllEdgesAndFaces  = false;
//        tMeshDataInput.Verbose                 = true;
//        tMeshDataInput.SpatialDim              = &tSpatialDim;
//        tMeshDataInput.ElemConn(0)             = &tInterfaceElementsInds;
//        tMeshDataInput.LocaltoGlobalElemMap(0) = &tInterfaceElementIds;
//        tMeshDataInput.NodeCoords              = &tNodeCoordinates;
//        tMeshDataInput.NodeProcOwner           = &tNodeOwner;
////        tMeshDataInput.LocaltoGlobalElemMap    = &aElemLocaltoGlobal;
//        tMeshDataInput.LocaltoGlobalNodeMap    = &tLocalToGlobalNodeMap;
//        tMeshDataInput.SetsInfo                = &tMtkMeshSets;
//        moris::mtk::Mesh* tMeshDataInt = moris::mtk::create_mesh( MeshType::STK, tMeshDataInput );
//        std::string tMeshOutputFile = tPrefix + "/xtk_test_unzip_int.e";
//        tMeshDataInt->create_output_mesh(tMeshOutputFile);
//

        // Output to file
        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/xtk_test_unzip.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);
        delete tCutMeshData;
        delete tMeshData;
//        delete tMeshDataInt;
    }
}
TEST_CASE("Unzipped element","[interface_element]")
{

    //    /**
    //     * This tests the functionality of the Interface_Element class.
    //     */
    //
    //    // initialize the element
    //    Interface_Element tInterfaceElement;
    //
    //    // Specify that element with index 10 and element with index 14 share an
    //    // interface  with their side ordinals 0 and 2 respectively
    //    moris_index tElementIndex0       = 10;
    //    moris_index tElementIndex1       = 14;
    //    moris_index tElement0SideOrdinal = 0;
    //    moris_index tElement1SideOrdinal = 2;
    //
    //    Matrix<IndexMat> tElementIndexPair({{tElementIndex0,tElementIndex1}});
    //    Matrix<IndexMat> tElementPairSideOrdinals({{tElement0SideOrdinal,tElement1SideOrdinal}});
    //
    //    // Set this information in the interface element
    //    tInterfaceElement.set_element_pair_and_side_ordinal(tElementIndexPair,tElementPairSideOrdinals);
    //
    //#ifdef DEBUG
    //    // Check the asserts work
    //    CHECK_THROWS(tInterfaceElement.set_element_pair_and_side_ordinal(tElementIndexPair,tElementPairSideOrdinals)); /*Trying to set the function twice*/
    //
    //    Interface_Element tInterfaceElement2;
    //    CHECK_THROWS(tInterfaceElement2.set_element_pair_and_side_ordinal({{0}},tElementPairSideOrdinals));/* incorrectly size vector */
    //    CHECK_THROWS(tInterfaceElement2.set_element_pair_and_side_ordinal(tElementIndexPair,{{0}}));/* incorrectly sized vector */
    //#endif
    //
    //    // Verify the element copied things correctly
    //    CHECK(all_true(tInterfaceElement.get_element_indices_pair() == tElementIndexPair));
    //    CHECK(all_true(tInterfaceElement.get_element_pair_side_ordinals() == tElementPairSideOrdinals));
    //


}


