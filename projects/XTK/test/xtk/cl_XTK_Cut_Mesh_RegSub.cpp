/*
 * cl_XTK_Cut_Mesh_RegSub.cpp
 *
 *  Created on: Feb 23, 2018
 *      Author: ktdoble
 */
#include "catch.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"
#include "fn_print.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Child_Mesh.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/fn_local_child_mesh_flood_fill.hpp"
#include "topology/cl_XTK_Hexahedron_8_Topology.hpp"


// DEBUGGING UTILITY INCLUDES
#include "tools/fn_tet_volume.hpp"
#include "geomeng/fn_Triangle_Geometry.hpp" // For surface normals
#include "mesh/fn_verify_tet_topology.hpp"

namespace xtk
{
TEST_CASE("Direct Testing of the regular subdivision","[NEW_REG_SUB_TEMPLATE]")
{
    // Model dimension
    size_t tModelDim = 3;

    // Set up global coordinates
    moris::Matrix< Default_Matrix_Real > tNodeCoords(15,3);
    tNodeCoords(0,0)  = 0.0; tNodeCoords(0,1)  = 0.0; tNodeCoords(0,2)  = 0.0;
    tNodeCoords(1,0)  = 1.0; tNodeCoords(1,1)  = 0.0; tNodeCoords(1,2)  = 0.0;
    tNodeCoords(2,0)  = 0.0; tNodeCoords(2,1)  = 1.0; tNodeCoords(2,2)  = 0.0;
    tNodeCoords(3,0)  = 1.0; tNodeCoords(3,1)  = 1.0; tNodeCoords(3,2)  = 0.0;
    tNodeCoords(4,0)  = 0.0; tNodeCoords(4,1)  = 0.0; tNodeCoords(4,2)  = 1.0;
    tNodeCoords(5,0)  = 1.0; tNodeCoords(5,1)  = 0.0; tNodeCoords(5,2)  = 1.0;
    tNodeCoords(6,0)  = 0.0; tNodeCoords(6,1)  = 1.0; tNodeCoords(6,2)  = 1.0;
    tNodeCoords(7,0)  = 1.0; tNodeCoords(7,1)  = 1.0; tNodeCoords(7,2)  = 1.0;
    tNodeCoords(8,0)  = 0.5; tNodeCoords(8,1)  = 0.0; tNodeCoords(8,2)  = 0.5;
    tNodeCoords(9,0)  = 1.0; tNodeCoords(9,1)  = 0.5; tNodeCoords(9,2)  = 0.5;
    tNodeCoords(10,0) = 0.5; tNodeCoords(10,1) = 1.0; tNodeCoords(10,2) = 0.5;
    tNodeCoords(11,0) = 0.0; tNodeCoords(11,1) = 0.5; tNodeCoords(11,2) = 0.5;
    tNodeCoords(12,0) = 0.5; tNodeCoords(12,1) = 0.5; tNodeCoords(12,2) = 0.0;
    tNodeCoords(13,0) = 0.5; tNodeCoords(13,1) = 0.5; tNodeCoords(13,2) = 1.0;
    tNodeCoords(14,0) = 0.5; tNodeCoords(14,1) = 0.5; tNodeCoords(14,2) = 0.5;

    // Initialize the Node Indices
    moris::Matrix< moris::IndexMat > tNodeIndex({{0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13, 14}});

    // Intialize Node Ids
    moris::Matrix< moris::IdMat > tNodeId({{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}});

    // Initialize the ancestry
    // Setup the parent tet ancestry (this should be 1 to 1)
    moris::Matrix< moris::IndexMat > tNodesAncestry({{0}});
    moris::Matrix< moris::IndexMat > tParentEdgeInds({{0,1,2,3,4,5,6,7,8,9,10,11}});
    moris::Matrix< Default_Matrix_Integer > tParentEdgeRanks(1,12,1);
    moris::Matrix< moris::IndexMat > tParentFaceInds({{0,1,2,3,4,5}});
    moris::Matrix< Default_Matrix_Integer > tParentFaceRanks(1,6,2);
    moris::Matrix< moris::IndexMat > tElementsAncestry({{0}});

    // Initialize Template
    Mesh_Modification_Template<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tRegSubTemplate(tElementsAncestry(0,0),
                                                                                                       0,
                                                                                                       tNodeIndex,
                                                                                                       tParentEdgeInds,
                                                                                                       tParentEdgeRanks,
                                                                                                       tParentFaceInds,
                                                                                                       tParentFaceRanks,
                                                                                                       TemplateType::REGULAR_SUBDIVISION_HEX8);

    // Initialize child mesh with template
    Child_Mesh_Test<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tRegSubChildMesh(tRegSubTemplate);


    // Check the volume
    moris::Matrix< moris::IndexMat > const & tElemToNode = tRegSubChildMesh.get_element_to_node();
    real tVolume = compute_volume_for_multiple_tets(tNodeCoords,tElemToNode);
    CHECK(approximate(tVolume,1.0));


    //
    moris::Matrix< moris::IndexMat > tElementPhase(1,24,0);

    moris::moris_index tMax = std::numeric_limits<moris::moris_index>::max();
    size_t tNumPhases = 2;

    moris::Matrix< moris::IndexMat > tActiveElements({{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}});
    moris::Matrix< moris::IndexMat > tIncludedElementMarker(1,24,1);

    // Run flood fill Algorithm to ensure that the floodfill can traverse the mesh
    moris::Matrix< moris::IndexMat > tElementSubphase = flood_fill( tRegSubChildMesh.get_element_to_element(),
                                                                    tElementPhase,
                                                                    tActiveElements,
                                                                    tIncludedElementMarker,
                                                                    tNumPhases,
                                                                    tMax,
                                                                    true);


    moris::Matrix< moris::IndexMat > tExpElementSubphase(1,24,0);
    CHECK(equal_to(tExpElementSubphase,tElementSubphase));

    // Make sure the tet4 topology is valid
    bool tValidTopo = verify_tet4_topology(tRegSubChildMesh.get_element_to_node(),
                                           tRegSubChildMesh.get_element_to_edge(),
                                           tRegSubChildMesh.get_element_to_face(),
                                           tRegSubChildMesh.get_edge_to_node(),
                                           tRegSubChildMesh.get_face_to_node());

    CHECK(tValidTopo);

    // Check that the volume is 1
    //TODO:

    // Parametric Coordinates (zeta, eta, xsi)
    // NOTE: these are ordered based on {0,1,3,2,4,6,5,7}
    moris::Matrix< moris::DDRMat > tParamCoords(15,3);

    // Base hex
    tParamCoords( 0,0) = -1.0; tParamCoords( 0,1) = -1.0; tParamCoords( 0,2) = -1.0;
    tParamCoords( 1,0) =  1.0; tParamCoords( 1,1) = -1.0; tParamCoords( 1,2) = -1.0;
    tParamCoords( 2,0) =  1.0; tParamCoords( 2,1) =  1.0; tParamCoords( 2,2) = -1.0;
    tParamCoords( 3,0) = -1.0; tParamCoords( 3,1) =  1.0; tParamCoords( 3,2) = -1.0;
    tParamCoords( 4,0) = -1.0; tParamCoords( 4,1) = -1.0; tParamCoords( 4,2) =  1.0;
    tParamCoords( 5,0) =  1.0; tParamCoords( 5,1) = -1.0; tParamCoords( 5,2) =  1.0;
    tParamCoords( 6,0) =  1.0; tParamCoords( 6,1) =  1.0; tParamCoords( 6,2) =  1.0;
    tParamCoords( 7,0) = -1.0; tParamCoords( 7,1) =  1.0; tParamCoords( 7,2) =  1.0;

    // Nodes at midside of hex faces
    tParamCoords(8,0)  =  0.0; tParamCoords(8,1)  = -1.0; tParamCoords(8,2)  = 0.0;
    tParamCoords(9,0)  =  1.0; tParamCoords(9,1)  =  0.0; tParamCoords(9,2)  = 0.0;
    tParamCoords(10,0) =  0.0; tParamCoords(10,1) =  1.0; tParamCoords(10,2) = 0.0;
    tParamCoords(11,0) = -1.0; tParamCoords(11,1) =  0.0; tParamCoords(11,2) = 0.0;
    tParamCoords(12,0) =  0.0; tParamCoords(12,1) =  0.0; tParamCoords(12,2) = -1.0;
    tParamCoords(13,0) =  0.0; tParamCoords(13,1) =  0.0; tParamCoords(13,2) = 1.0;

    // Nodes at center of hex element
    tParamCoords(14,0) = 0.0; tParamCoords(14,1) = 0.0; tParamCoords(14,2) = 0.0;

    // allocate space
    tRegSubChildMesh.allocate_parametric_coordinates(15);

    // Add parametric coordinate
    tRegSubChildMesh.add_node_parametric_coordinate(tRegSubChildMesh.get_node_indices(),tParamCoords);

    // Check the parametric coordinates are added as expected
    moris::Matrix<moris::DDRMat>   const & tCMParamCoords = tRegSubChildMesh.get_parametric_coordinates();
    CHECK(all_true(tParamCoords == tCMParamCoords));

    // Verify we can retrieve the correct global coordinate
    // by interpolating from the base element and parametric
    // coordinate of a node

    // Get child node indices from the child mesh (note these aren't guaranteed to be monotonically increasing)
     moris::Matrix<moris::IndexMat> const & tNodeIndicesOfCM = tRegSubChildMesh.get_node_indices();

     // Topology of the base hex
     Hexahedron_8_Topology<real,size_t,Default_Matrix_Real,Default_Matrix_Integer>
     tHex8Topo({{0,1,2,3,4,5,6,7}});

     // Get the basis of the hex8
     Basis_Function<real, Default_Matrix_Real> const &
     tHex8Basis = tHex8Topo.get_basis_function();

     // Coordinates of the base hex8 (note these are in a different order from the tNodeCoords)
     moris::Matrix<moris::DDRMat> tHex8Coords(8,3);
     tHex8Coords.set_row(0,tNodeCoords.get_row(tNodeIndicesOfCM(0)));
     tHex8Coords.set_row(1,tNodeCoords.get_row(tNodeIndicesOfCM(1)));
     tHex8Coords.set_row(2,tNodeCoords.get_row(tNodeIndicesOfCM(2)));
     tHex8Coords.set_row(3,tNodeCoords.get_row(tNodeIndicesOfCM(3)));
     tHex8Coords.set_row(4,tNodeCoords.get_row(tNodeIndicesOfCM(4)));
     tHex8Coords.set_row(5,tNodeCoords.get_row(tNodeIndicesOfCM(5)));
     tHex8Coords.set_row(6,tNodeCoords.get_row(tNodeIndicesOfCM(6)));
     tHex8Coords.set_row(7,tNodeCoords.get_row(tNodeIndicesOfCM(7)));

     // iterate over nodes
     size_t tNumNodes = tRegSubChildMesh.get_num_entities(EntityRank::NODE);

     // Allocate a basis function weight matrix
     moris::Matrix<moris::DDRMat> tBasisWeights(1,8);

     // tolerance for difference between coordinates
      real tTol = 1e-12;

    for(size_t i= 0; i<tNumNodes; i++)
    {
        // Node index
        moris::moris_index tNodeIndex = tNodeIndicesOfCM(i);

        // Get the nodes parametric coordinate
        moris::Matrix<moris::DDRMat> tNodeParamCoord = tRegSubChildMesh.get_parametric_coordinates(tNodeIndex);

        // Get the basis function values at this point
        tHex8Basis.evaluate_basis_function(tNodeParamCoord,tBasisWeights);

        // Evaluate the nodes global coordinate from the basis weights
        moris::Matrix<moris::DDRMat> tInterpNodeCoord = tBasisWeights*tHex8Coords;

        // Verify the interpolated coordinate is equal to the node coordinate row
        CHECK(moris::norm(tInterpNodeCoord - tNodeCoords.get_row(tNodeIndex)) < tTol);
    }
}

}


