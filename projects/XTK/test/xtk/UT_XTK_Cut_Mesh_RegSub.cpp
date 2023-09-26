/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Mesh_RegSub.cpp
 *
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
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "cl_XTK_Output_Options.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_generate_element_to_element.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"
#include "cl_XTK_Hexahedron_8_Topology.hpp"

// DEBUGGING UTILITY INCLUDES
#include "fn_tet_volume.hpp"
#include "fn_verify_tet_topology.hpp"
#include "fn_GEN_Triangle_Geometry.hpp"

namespace xtk
{

// note not written to be very efficient
bool
node_is_on_face(moris_index     aFaceIndex,
                Row_View_Real const & aNodeCoord)
{
    Cell<Matrix<DDRMat>> tFaceNodeCoords(6);

    // face node coordinates
    tFaceNodeCoords(0) = Matrix<DDRMat>({{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}});
    tFaceNodeCoords(1) = Matrix<DDRMat>({{1, 0, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}});
    tFaceNodeCoords(2) = Matrix<DDRMat>({{0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {0, 1, 1}});
    tFaceNodeCoords(3) = Matrix<DDRMat>({{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}});
    tFaceNodeCoords(4) = Matrix<DDRMat>({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}});
    tFaceNodeCoords(5) = Matrix<DDRMat>({{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}});

    // compute bounding box
    Matrix<DDRMat> tFaceXCoords = tFaceNodeCoords(aFaceIndex).get_column(0);
    Matrix<DDRMat> tFaceYCoords = tFaceNodeCoords(aFaceIndex).get_column(1);
    Matrix<DDRMat> tFaceZCoords = tFaceNodeCoords(aFaceIndex).get_column(2);
    real tXHigh = tFaceXCoords.max();
    real tXLow  = tFaceXCoords.min();
    real tYHigh = tFaceYCoords.max();
    real tYLow  = tFaceYCoords.min();
    real tZHigh = tFaceZCoords.max();
    real tZLow  = tFaceZCoords.min();

    // check x coordinate of node is within bounding box
    bool tXInBox = false;
    if(aNodeCoord(0) <= tXHigh && aNodeCoord(0) >= tXLow )
    {
        tXInBox = true;
    }

    // check y coordinate of node is within bounding box
    bool tYInBox = false;
    if(aNodeCoord(1) <= tYHigh && aNodeCoord(1) >= tYLow )
    {
        tYInBox = true;
    }

    // check z coordinate of node is within bounding box
    bool tZInBox = false;
    if(aNodeCoord(2) <= tZHigh && aNodeCoord(2) >= tZLow )
    {
        tZInBox = true;
    }

    return (tXInBox && tYInBox && tZInBox);
}

TEST_CASE("Direct Testing of the regular subdivision","[NEW_REG_SUB_TEMPLATE]")
{

    // Set up global coordinates
    moris::Matrix< moris::DDRMat > tNodeCoords(15,3);
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
    moris::Matrix< moris::IndexMat > tParentNodeInds({{0,1,2,3,4,5,6,7}});
    moris::Matrix< moris::DDSTMat >  tParentNodeRanks(1,8,0);
    moris::Matrix< moris::IndexMat > tParentEdgeInds({{0,1,2,3,4,5,6,7,8,9,10,11}});
    moris::Matrix< moris::DDSTMat >  tParentEdgeRanks(1,12,1);
    moris::Matrix< moris::IndexMat > tParentFaceInds({{0,1,2,3,4,5}});
    moris::Matrix< moris::DDSTMat >  tParentFaceRanks(1,6,2);
    moris::Matrix< moris::IndexMat > tElementsAncestry({{0}});

    // Initialize Template
    Mesh_Modification_Template tRegSubTemplate(tElementsAncestry(0,0),
                                               0,
                                               tNodeIndex,
                                               tParentNodeInds,
                                               tParentNodeRanks,
                                               tParentEdgeInds,
                                               tParentEdgeRanks,
                                               tParentFaceInds,
                                               tParentFaceRanks,
                                               TemplateType::REGULAR_SUBDIVISION_HEX8);

    // Initialize child mesh with template
    Child_Mesh tRegSubChildMesh(tRegSubTemplate);

    // Check the volume
    moris::Matrix< moris::IndexMat > const & tElemToNode = tRegSubChildMesh.get_element_to_node();
    real tVolume = moris::ge::compute_volume_for_multiple_tets(tNodeCoords,tElemToNode);
    CHECK(approximate(tVolume,1.0));

    //
    moris::Matrix< moris::IndexMat > tElementPhase(1,24,0);

    moris::moris_index tMax = std::numeric_limits<moris::moris_index>::max();
    size_t tNumPhases = 2;

    moris::Matrix< moris::IndexMat > tActiveElements({{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}});
    moris::Matrix< moris::IndexMat > tIncludedElementMarker(1,24,1);
    moris::moris_index tMaxFloodFill = 0;
    // Run flood fill Algorithm to ensure that the floodfill can traverse the mesh
    moris::Matrix< moris::IndexMat > tElementSubphase = flood_fill( tRegSubChildMesh.get_element_to_element(),
                                                                    tElementPhase,
                                                                    tActiveElements,
                                                                    tIncludedElementMarker,
                                                                    tNumPhases,
                                                                    tMax,
                                                                    tMaxFloodFill,
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
    tRegSubChildMesh.allocate_parametric_coordinates(15,3);

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
     Hexahedron_8_Topology tHex8Topo({{0,1,2,3,4,5,6,7}});

     // Get the basis of the hex8
     Basis_Function const & tHex8Basis = tHex8Topo.get_basis_function();

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
     size_t tNumNodes = tRegSubChildMesh.get_num_entities(mtk::EntityRank::NODE);

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

    // verify edge ancestry
    Matrix<IndexMat> const & tEdgeToNode = tRegSubChildMesh.get_edge_to_node();
    moris::Matrix< moris::IndexMat > const & tEdgeParentIndices = tRegSubChildMesh.get_edge_parent_inds();
    moris::Matrix< moris::DDSTMat >  const & tEdgeParentRanks   = tRegSubChildMesh.get_edge_parent_ranks();

    for(moris::uint i = 0; i <tRegSubChildMesh.get_num_entities(mtk::EntityRank::EDGE); i++)
    {
        if(tEdgeParentRanks(i) == 2)
        {
            moris_index tFaceIndex = tEdgeParentIndices(i);
            Row_View_Real const & tEdgeNodeCoord0 = tNodeCoords.get_row(tEdgeToNode(i,0));
            Row_View_Real const & tEdgeNodeCoord1 = tNodeCoords.get_row(tEdgeToNode(i,1));
            CHECK((node_is_on_face(tFaceIndex,tEdgeNodeCoord0) && node_is_on_face(tFaceIndex,tEdgeNodeCoord1)));
        }
    }

}

}

