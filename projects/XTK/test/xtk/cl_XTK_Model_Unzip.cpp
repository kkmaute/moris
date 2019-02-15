
#include "catch.hpp"

// Xtk includes
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Interface_Element.hpp"

// Geometry
#include "cl_Sphere.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

// Linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

using namespace moris;
using namespace xtk;

TEST_CASE("XTK Model Unzip Interface","[unzip_xtk]")
{
    if(par_size() == 1 || par_size() ==2)
    {
        real tRadius  = 0.75;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 1.0;
        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

        // Create Mesh --------------------------------------------------------------------
        std::string tMeshFileName = "generated:1x1x2";
        moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName );

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        moris::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                      Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        // Test interface information prior to unzipping
        Background_Mesh const & tBM = tXTKModel.get_background_mesh();

        // get the interface nodes
        moris::Matrix< moris::IndexMat > tGoldInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35}});
        moris::Matrix< moris::IndexMat > tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);
        CHECK(all_true(tGoldInterfaceNodes == tInterfaceNodes));

        // Unzip interface
        tXTKModel.unzip_interface();

        // Now there should be twice as many nodes
        moris::Matrix< moris::IndexMat > tGoldUnzippedInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
            36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46}});
        tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);
        CHECK(all_true(tGoldUnzippedInterfaceNodes == tInterfaceNodes));


        // Output to file
        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/xtk_test_unzip.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);
        delete tCutMeshData;
        delete tMeshData;
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


