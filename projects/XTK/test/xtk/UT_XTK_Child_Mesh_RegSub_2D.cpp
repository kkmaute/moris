/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Child_Mesh_RegSub_2D.cpp
 *
 */

#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"
#include "fn_local_child_mesh_flood_fill.hpp"

#include "cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "cl_XTK_Child_Mesh.hpp"
#include "fn_GEN_Triangle_Geometry.hpp"
#include "cl_XTK_Quad_4_Topology.hpp"

#include "fn_verify_tri_topology.hpp"

namespace moris::xtk
{

    TEST_CASE( "Regular Subdivision QUAD4", "[REG_SUB_TEMPLATE_QUAD4]" )
    {
        // Spatial dimension
        uint tSpatialDim = 2;

        // Number of nodes in template
        uint tNumNodesInTemplate = 5;

        // Node Coordinates
        Matrix< DDRMat > tNodeCoords( tNumNodesInTemplate, tSpatialDim );
        tNodeCoords( 0, 0 ) = 0.0;
        tNodeCoords( 0, 1 ) = 0.0;
        tNodeCoords( 1, 0 ) = 1.0;
        tNodeCoords( 1, 1 ) = 0.0;
        tNodeCoords( 2, 0 ) = 0.0;
        tNodeCoords( 2, 1 ) = 1.0;
        tNodeCoords( 3, 0 ) = 1.0;
        tNodeCoords( 3, 1 ) = 1.0;
        tNodeCoords( 4, 0 ) = 0.5;
        tNodeCoords( 4, 1 ) = 0.5;

        // All nodes in template
        Matrix< IndexMat > tNodeIndex( { { 0, 1, 3, 2, 4 } } );

        // All node ids in template
        Matrix< IndexMat > tNodeIds( { { 1, 2, 4, 3, 5 } } );

        // Add ancestry information (on board we called this template parent entity parents
        // Do for nodes, edges, cells
        Matrix< IndexMat > tParentNodeInds = { { 0, 1, 2, 3 } };
        Matrix< DDSTMat >  tParentNodeRanks( 1, 4, 0 );
        Matrix< IndexMat > tParentEdgeInds = { { 0, 1, 2, 3 } };
        Matrix< DDSTMat >  tParentEdgeRanks( 1, 4, 1 );
        Matrix< IndexMat > tElementAncestry = { { 0 } };

        // Setup mesh modification template
        // Initialize Template
        Mesh_Modification_Template tRegSubTemplate( tElementAncestry( 0, 0 ),
                0,
                tNodeIndex,
                tParentNodeInds,
                tParentNodeRanks,
                tParentEdgeInds,
                tParentEdgeRanks,
                { {} },
                { {} },
                TemplateType::REGULAR_SUBDIVISION_QUAD4 );

        // Initialize child mesh with template
        Child_Mesh tRegSubChildMesh( tRegSubTemplate );

        // Check the area
        moris::Matrix< IndexMat > const &tElemToNode = tRegSubChildMesh.get_element_to_node();
        real                             tArea       = compute_area_for_multiple_triangles( tNodeCoords, tElemToNode );
        CHECK( approximate( tArea, 1.0 ) );

        // Set up phase data
        Vector< moris_index > tElementPhase( 4, 0 );

        moris::moris_index tMax       = std::numeric_limits< moris::moris_index >::max();

        Vector< moris_index > tActiveElements = { 0, 1, 2, 3 };
        Vector< moris_index > tIncludedElementMarker( 4, 1 );
        moris::moris_index    tMaxFloodFill = 0;
        // Run flood fill Algorithm to ensure that the flood-fill can traverse the mesh
        moris::Matrix< IndexMat > tElementSubphase = mtk::flood_fill(
                tRegSubChildMesh.get_element_to_element(),
                tElementPhase,
                tActiveElements,
                tIncludedElementMarker,
                tMax,
                tMaxFloodFill,
                true );


        moris::Matrix< IndexMat > tExpectedElementSubphase( 4, 1, 0 );
        CHECK( all_true( tExpectedElementSubphase == tElementSubphase ) );

        // Make sure the tri3 topology is valid
        bool tValidTopo = verify_tri3_topology( tRegSubChildMesh.get_element_to_node(),
                tRegSubChildMesh.get_element_to_edge(),
                tRegSubChildMesh.get_edge_to_node() );

        CHECK( tValidTopo );

        // Parametric Coordinates (zeta, eta, xsi)
        moris::Matrix< DDRMat > tParamCoords( 5, 2 );

        // Base quad
        tParamCoords( 0, 0 ) = -1.0;
        tParamCoords( 0, 1 ) = -1.0;
        tParamCoords( 1, 0 ) = 1.0;
        tParamCoords( 1, 1 ) = -1.0;
        tParamCoords( 2, 0 ) = 1.0;
        tParamCoords( 2, 1 ) = 1.0;
        tParamCoords( 3, 0 ) = -1.0;
        tParamCoords( 3, 1 ) = 1.0;

        // Nodes at center of quad element
        tParamCoords( 4, 0 ) = 0.0;
        tParamCoords( 4, 1 ) = 0.0;

        // allocate space
        tRegSubChildMesh.allocate_parametric_coordinates( 5, 2 );

        // Add parametric coordinate
        tRegSubChildMesh.add_node_parametric_coordinate( tRegSubChildMesh.get_node_indices(), tParamCoords );

        // Check the parametric coordinates are added as expected
        moris::Matrix< DDRMat > const &tCMParamCoords = tRegSubChildMesh.get_parametric_coordinates();
        CHECK( all_true( tParamCoords == tCMParamCoords ) );

        // Verify we can retrieve the correct global coordinate
        // by interpolating from the base element and parametric
        // coordinate of a node

        // Get child node indices from the child mesh (note these aren't guaranteed to be monotonically increasing)
        moris::Matrix< IndexMat > const &tNodeIndicesOfCM = tRegSubChildMesh.get_node_indices();

        // Topology of the base quad
        Quad_4_Topology tQuad4Topo( { { 0, 1, 2, 3 } } );

        // Get the basis of the quad4
        Basis_Function const &tQuad4Basis = tQuad4Topo.get_basis_function();

        // Coordinates of the base quad4 (note these are in a different order from the tNodeCoords)
        moris::Matrix< DDRMat > tQuad4Coords( 4, 2 );
        tQuad4Coords.set_row( 0, tNodeCoords.get_row( tNodeIndicesOfCM( 0 ) ) );
        tQuad4Coords.set_row( 1, tNodeCoords.get_row( tNodeIndicesOfCM( 1 ) ) );
        tQuad4Coords.set_row( 2, tNodeCoords.get_row( tNodeIndicesOfCM( 2 ) ) );
        tQuad4Coords.set_row( 3, tNodeCoords.get_row( tNodeIndicesOfCM( 3 ) ) );

        // iterate over nodes
        size_t tNumNodes = tRegSubChildMesh.get_num_entities( mtk::EntityRank::NODE );

        // Allocate a basis function weight matrix
        moris::Matrix< DDRMat > tBasisWeights( 1, 4 );

        // tolerance for difference between coordinates
        real tTol = 1e-12;

        for ( size_t i = 0; i < tNumNodes; i++ )
        {
            // Node index
            moris::moris_index tNodeIndex = tNodeIndicesOfCM( i );

            // Get the nodes parametric coordinate
            moris::Matrix< DDRMat > tNodeParamCoord = tRegSubChildMesh.get_parametric_coordinates( tNodeIndex );

            // Get the basis function values at this point
            tQuad4Basis.evaluate_basis_function( tNodeParamCoord, tBasisWeights );

            // Evaluate the nodes global coordinate from the basis weights
            moris::Matrix< DDRMat > tInterpNodeCoord = tBasisWeights * tQuad4Coords;

            // Verify the interpolated coordinate is equal to the node coordinate row
            CHECK( moris::norm( tInterpNodeCoord - tNodeCoords.get_row( tNodeIndex ) ) < tTol );
        }
    }

}    // namespace moris::xtk
