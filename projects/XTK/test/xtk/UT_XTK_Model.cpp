/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/fn_write_element_ownership_as_field.hpp"

#include "topology/cl_XTK_Hexahedron_8_Basis_Function.hpp"

#include "cl_MTK_Visualization_STK.hpp"

#include "cl_GEN_Sphere.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
#include "geomeng/fn_Triangle_Geometry.hpp" // For surface normals

// Linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"

#include "cl_Profiler.hpp"
#include "Child_Mesh_Verification_Utilities.hpp"


#include "fn_compute_interface_surface_area.hpp"

namespace xtk
{


TEST_CASE("Regular Subdivision Method","[XTK] [REGULAR_SUBDIVISION_MODEL]")
        {
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize==1)
    {
        // Geometry Engine Setup -----------------------
        // Using a Levelset Sphere as the Geometry

        real tRadius = 0.7;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 0;
        moris::ge::Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

        moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:1x1x1";
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
        tXTKModel.decompose(tDecompositionMethods);

        // Access the decomposed XTK Mesh
        Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

        // Do some testing
        size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
        size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

        // 1 element was subdivided
        CHECK(tNumNodesAfterDecompositionXTK == 15);
        CHECK(tNumElementsAfterDecompositionXTK == 24);

        moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

        moris::Matrix< moris::DDRMat > tExpectedNodeCoordinates(
           {{0, 0, 0},
            {1, 0, 0},
            {0, 1, 0},
            {1, 1, 0},
            {0, 0, 1},
            {1, 0, 1},
            {0, 1, 1},
            {1, 1, 1},
            {0.5, 0, 0.5},
            {1, 0.5, 0.5},
            {0.5, 1, 0.5},
            {0, 0.5, 0.5},
            {0.5, 0.5, 0},
            {0.5, 0.5, 1},
            {0.5, 0.5, 0.5}});
        CHECK(equal_to(tNodeCoordinates,tExpectedNodeCoordinates));

        verify_child_mesh_ancestry(tXTKModel.get_background_mesh(),
                                   tCutMesh);

        // Verify parametric coordinates reproduce the same results as tabulated in global node coordinates
        // Basis function to use for hex8
        Hexahedron_8_Basis_Function tHex8Basis;

        // Allocate a basis function weight matrix
        moris::Matrix< moris::DDRMat > tBasisWeights(1,8);

        // tolerance for difference between coordinates
        real tTol = 1e-12;

        // Iterate over child meshes
        for(size_t iCM = 0; iCM < tCutMesh.get_num_child_meshes(); iCM++)
        {
            // Get reference to child mesh
            Child_Mesh const & tChildMesh = tCutMesh.get_child_mesh(iCM);

            // iterate over nodes
            size_t tNumNodes = tChildMesh.get_num_entities(EntityRank::NODE);

            // Get child node indices from the child mesh (note these aren't guaranteed to be monotonically increasing)
            moris::Matrix<moris::IndexMat> const & tNodeIndicesOfCM = tChildMesh.get_node_indices();

            // Parent element index
            moris::moris_index tParentIndex = tChildMesh.get_parent_element_index();

            // Nodes attached to parent element
            moris::Matrix< moris::IndexMat > tNodesAttachedToParentElem = tMeshData->get_entity_connected_to_entity_loc_inds(tParentIndex,moris::EntityRank::ELEMENT,moris::EntityRank::NODE);

            // Get the node coordinates
            moris::Matrix< moris::DDRMat > tHex8NodeCoords = tXTKModel.get_background_mesh().get_selected_node_coordinates_loc_inds(tNodesAttachedToParentElem);

            for(size_t i= 0; i<tNumNodes; i++)
            {
                // Node index
                moris::moris_index tNodeIndex = tNodeIndicesOfCM(i);

                // Get the nodes parametric coordinate
                moris::Matrix<moris::DDRMat> tNodeParamCoord = tChildMesh.get_parametric_coordinates(tNodeIndex);

                // Get the basis function values at this point
                tHex8Basis.evaluate_basis_function(tNodeParamCoord,tBasisWeights);

                // Evaluate the nodes global coordinate from the basis weights
                moris::Matrix<moris::DDRMat> tInterpNodeCoord = tBasisWeights*tHex8NodeCoords;

                // Verify the interpolated coordinate is equal to the node coordinate row
                CHECK(norm(tInterpNodeCoord - tNodeCoordinates.get_row(tNodeIndex)) < tTol);
            }
        }

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

        std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_regular_subdivision.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);
        delete tMeshData;
        delete tCutMeshData;
    }
}
TEST_CASE("Regular Subdivision and Nodal Hierarchy Subdivision","[XTK] [CONFORMAL_MODEL]")
{
    int tProcRank = moris::par_rank();
    int tProcSize = moris::par_size();

    if(tProcSize<=2)
    {
            // Geometry Engine Setup ---------------------------------------------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius  = 0.25;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0.0;
            moris::ge::Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
            moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
            moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:1x1x4";
            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );
            std::string tBMOutputFile ="./xtk_exo/xtk_test_output_conformal_bm.e";
            tMeshData->create_output_mesh(tBMOutputFile);

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
            tXTKModel.mVerbose  =  false;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            // output to exodus file ----------------------------------------------------------
            Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = true;
            tOutputOptions.mAddSideSets = true;

            moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

            std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_conformal.e";
            tCutMeshData->create_output_mesh(tMeshOutputFile);

            // Access the Cut Mesh-------------------------------------------------------------
            Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

            // Do some testing
            size_t tNumNodesAfterDecompositionSTK    = tXTKModel.get_background_mesh().get_num_entities(EntityRank::NODE);
            size_t tNumElementsAfterDecompositionSTK = tXTKModel.get_background_mesh().get_num_entities(EntityRank::ELEMENT);
            size_t tNumNodesAfterDecompositionXTK    = tCutMesh.get_num_entities(EntityRank::NODE);
            size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

            if(par_size() == 1)
            {
                CHECK(tNumNodesAfterDecompositionXTK == 22 );
                CHECK(tNumElementsAfterDecompositionXTK == 42);

                CHECK(tNumElementsAfterDecompositionSTK == 46);
                // But it does add nodes because this is where the coordinates are stored
                CHECK(tNumNodesAfterDecompositionSTK == 34);
            }

            moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

            // verify ancestry of child mesh
            verify_child_mesh_ancestry(tXTKModel.get_background_mesh(),
                                       tCutMesh);

            tXTKModel.perform_basis_enrichment(EntityRank::NODE);

            // Verify parametric coordinates reproduce the same results as tabulated in global node coordinates

            // Basis function to use for hex8
            Hexahedron_8_Basis_Function tHex8Basis;

            // Allocate a basis function weight matrix
              moris::Matrix< moris::DDRMat > tBasisWeights(1,8);

             // tolerance for difference between coordinates
              real tTol = 1e-12;

             // Iterate over child meshes
            for(size_t iCM = 0; iCM < tCutMesh.get_num_child_meshes(); iCM++)
            {
                // Get reference to child mesh
                Child_Mesh const & tChildMesh = tCutMesh.get_child_mesh(iCM);

                // iterate over nodes
                size_t tNumNodes = tChildMesh.get_num_entities(EntityRank::NODE);

                // Get child node indices from the child mesh (note these aren't guaranteed to be monotonically increasing)
                 moris::Matrix<moris::IndexMat> const & tNodeIndicesOfCM = tChildMesh.get_node_indices();

                  // Parent element index
                moris::moris_index tParentIndex = tChildMesh.get_parent_element_index();

                // Nodes attached to parent element
                moris::Matrix< moris::IndexMat > tNodesAttachedToParentElem = tMeshData->get_entity_connected_to_entity_loc_inds(tParentIndex,moris::EntityRank::ELEMENT,moris::EntityRank::NODE);

                // Get the node coordinates
                moris::Matrix< moris::DDRMat > tHex8NodeCoords = tXTKModel.get_background_mesh().get_selected_node_coordinates_loc_inds(tNodesAttachedToParentElem);

               for(size_t i= 0; i<tNumNodes; i++)
               {
                   // Node index
                   moris::moris_index tNodeIndex = tNodeIndicesOfCM(i);

                   // Get the nodes parametric coordinate
                   moris::Matrix<moris::DDRMat> tNodeParamCoord = tChildMesh.get_parametric_coordinates(tNodeIndex);

                   // Get the basis function values at this point
                   tHex8Basis.evaluate_basis_function(tNodeParamCoord,tBasisWeights);

                   // Evaluate the nodes global coordinate from the basis weights
                   moris::Matrix<moris::DDRMat> tInterpNodeCoord = tBasisWeights*tHex8NodeCoords;

                   // Verify the interpolated coordinate is equal to the node coordinate row
                   CHECK(norm(tInterpNodeCoord - tNodeCoordinates.get_row(tNodeIndex)) < tTol);
               }

            }

            // verify the interface nodes have level set values sufficiently close to 0
            // when interpolated to
            moris::Matrix<moris::IndexMat> tInterfaceNodes = tXTKModel.get_background_mesh().get_interface_nodes_loc_inds(0);

            // Iterate over interface nodes
            for(size_t iIN = 0; iIN < tInterfaceNodes.numel(); iIN++)
            {
                // Get the child meshes these nodes belong to
                moris::Matrix<moris::IndexMat> tCMIndices = tXTKModel.get_background_mesh().get_node_child_mesh_assocation(tInterfaceNodes(iIN));

                // Iterate over child meshes this node belongs to
                for(size_t iCM = 0; iCM<tCMIndices.numel(); iCM++)
                {
                    // Get reference to child mesh
                    Child_Mesh const & tChildMesh = tCutMesh.get_child_mesh(tCMIndices(iCM));

                    // Parent element index
                    moris::moris_index tParentIndex = tChildMesh.get_parent_element_index();

                    // Nodes attached to parent element
                    moris::Matrix< moris::IndexMat > tNodesAttachedToParentElem = tMeshData->get_entity_connected_to_entity_loc_inds(tParentIndex,moris::EntityRank::ELEMENT,moris::EntityRank::NODE);

                    // Get the level set values of each node on the hex8
                    moris::Matrix< moris::DDRMat > tHex8LSVs(8,1);
                    for(size_t iNode = 0; iNode<8; iNode++)
                    {
                        tHex8LSVs(iNode)  = tXTKModel.get_geom_engine().get_entity_phase_val(tNodesAttachedToParentElem(iNode),0);
                    }

                    // Get the interface node parametric coordinate in iCM
                    moris::Matrix<moris::DDRMat> tNodeParamCoord = tChildMesh.get_parametric_coordinates(tInterfaceNodes(iIN));

                    // Get the basis function values at this point
                    tHex8Basis.evaluate_basis_function(tNodeParamCoord,tBasisWeights);

                    // Evaluate the nodes global coordinate from the basis weights
                    moris::Matrix<moris::DDRMat> tInterfaceLSV = tBasisWeights*tHex8LSVs;

                    // Verify it is  approximately 0.0
                    CHECK(approximate(tInterfaceLSV(0),0.0) );

                }
            }


            delete tCutMeshData;
            delete tMeshData;
        }
    }

TEST_CASE("Propagate Mesh Sets","[SET_PROPOGATION]")
{
    /*
     * Loads an exodus file with a Block Set and Side Set already populated
     * Performs regular subdivison method and then checks to see if the
     * children XTK mesh entities have the same as their parents
     */

    /*
     * Set up:
     * Geometry,
     * Geometry Engine,
     * Mesh
     */

    real tRadius = 5.1111;
    real tXCenter = 0.0;
    real tYCenter = 0.0;
    real tZCenter = 0.0;
    moris::ge::Sphere tLevelSetSphere(tRadius,tXCenter,tYCenter,tZCenter);
    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tLevelSetSphere,tPhaseTable);

    tGeometryEngine.mThresholdValue = 0.0;
    tGeometryEngine.mComputeDxDp = false;

    /*
     * Load Mesh which has 3 block sets. These blocks are named:
     *  - top_bread
     *  - meat
     *  - bottom_bread
     *
     * Side Sets will eventually be named
     *  - top_crust
     *  - bottom_crust
     */
    std::string tPrefix;
    tPrefix = std::getenv("MORISROOT");
    std::string tMeshFileName = tPrefix + "/projects/XTK/test/test_exodus_files/sandwich.e";
    moris::Cell<std::string> tFieldNames;


    // add parallel fields to the mesh
    moris::mtk::Visualization_STK tVizTool;
    moris::mtk::MtkFieldsInfo* tFieldsInfo = tVizTool.setup_parallel_cell_fields_for_declaration();

    // Declare some supplementary fields
    moris::mtk::MtkMeshData tMeshDataInput;
    tMeshDataInput.FieldsInfo = tFieldsInfo;

    // fill in the parallel fields
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName,  &tMeshDataInput  );

    tVizTool.populate_parallel_cell_fields_on_mesh(tMeshData);

    std::string tBackgroundfile = "./xtk_exo/xtk_ut_set_propogation_bm.e";
    tMeshData->create_output_mesh(tBackgroundfile);


    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};

    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mVerbose  =  false;
    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);

//    tXTKModel.unzip_interface();

    moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();
    moris::real tMySurfaceArea = compute_interface_surface_area(tNodeCoords,tXTKModel);
    // Collect all volumes
    moris::real tGlbSurf = 0.0;
    sum_all(tMySurfaceArea,tGlbSurf);


    tXTKModel.perform_basis_enrichment(EntityRank::NODE);

    /*
     * Get the output mesh and write to exodus file
     */

    Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    // Set the sensitivity field names
    tOutputOptions.mPackageDxDpSparsely = true;

    // Add field for enrichment
    tOutputOptions.mIntElementExternalFieldNames = {"owner"};
    tOutputOptions.mInternalUseFlag = true;

    tOutputOptions.mAddParallelFields = true;

    moris::mtk::Mesh* tOutputMeshData = tXTKModel.get_output_mesh(tOutputOptions);

   std::string tMeshOutputFile = "./xtk_exo/xtk_ut_set_propogation.e";
   tOutputMeshData->create_output_mesh(tMeshOutputFile);

   delete tMeshData;
   delete tOutputMeshData;

}

}
