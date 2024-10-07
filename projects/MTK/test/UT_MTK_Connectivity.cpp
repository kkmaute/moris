/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Connectivity.cpp
 *
 */

#include "catch.hpp"

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh_Core.hpp"

#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Cell_Info_Hex8.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_MTK_Cell_Info_Quad4.hpp"
#include "cl_MTK_Cell_Info_Tri3.hpp"
#include "fn_cross.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_sort.hpp"

namespace moris::mtk
{

Mesh* generate_single_element_hex8()
{
    // Set up a simple mtk mesh
    uint tNumDim = 3;

    // Just a single Hex8
    uint aNumElemTypes = 1;
    Matrix< IdMat > tElemConnHex8 = {{1, 2, 4, 3, 5, 6, 8, 7}};
    Matrix< IdMat > tElemLocaltoGlobalHex8 = {{1}};

    Matrix< DDRMat > tCoords = {{0.0, 0.0, 0.0},
                                {1.0, 0.0, 0.0},
                                {0.0, 1.0, 0.0},
                                {1.0, 1.0, 0.0},
                                {0.0, 0.0, 1.0},
                                {1.0, 0.0, 1.0},
                                {0.0, 1.0, 1.0},
                                {1.0, 1.0, 1.0}};

    Matrix< IdMat > tNodeLocaltoGlobal {{1, 2, 3, 4, 5, 6, 7, 8}};

    // Define block sets
    MtkSetsInfo tMtkMeshSets;

    Matrix< IdMat > tCellIdsBS({{1}});
    MtkBlockSetInfo tBlockSet;
    tBlockSet.mCellIdsInSet = &tCellIdsBS;
    tBlockSet.mBlockSetName = "blockset";
    tBlockSet.mBlockSetTopo = CellTopology::HEX8;

    tMtkMeshSets.add_block_set(&tBlockSet);

    // Define side sets
    Matrix< IdMat > tElemIdsAndSideOrdsSS = {{1, 2},
                                             {1, 3},
                                             {1, 4},
                                             {1, 5}};

    MtkSideSetInfo tSideSetStruct;
    tSideSetStruct.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS;
    tSideSetStruct.mSideSetName = "sideset";

    tMtkMeshSets.add_side_set(&tSideSetStruct);

    // Define node sets
    Matrix< IdMat > tNodeIdsNS = {{1}, {2}, {4}, {5}};

    MtkNodeSetInfo tNodeSet;
    tNodeSet.mNodeIds = &tNodeIdsNS;
    tNodeSet.mNodeSetName = "nodeset";

    tMtkMeshSets.add_node_set(&tNodeSet);

    // Create scalar field
    Scalar_Field_Info< DDRMat > tNodeField;

    Matrix < DDRMat > tNodeFieldData = {{0.0, 1.0, 0.5, 0.5, 0.5, 0.5, 1.0, 0.0}};

    tNodeField.set_field_name("nodefield");
    tNodeField.set_field_entity_rank(EntityRank::NODE);
    tNodeField.add_field_data(&tNodeLocaltoGlobal, &tNodeFieldData);

    // Collect field info
    MtkFieldsInfo tFieldsInfo;
    tFieldsInfo.mRealScalarFields.push_back(&tNodeField);

    // Coalesce mesh data
    MtkMeshData aMeshData(aNumElemTypes);
    aMeshData.SpatialDim = &tNumDim;
    aMeshData.ElemConn(0) = &tElemConnHex8;
    aMeshData.NodeCoords = &tCoords;
    aMeshData.SetsInfo = &tMtkMeshSets;
    aMeshData.LocaltoGlobalElemMap(0) = &tElemLocaltoGlobalHex8;
    aMeshData.LocaltoGlobalNodeMap = &tNodeLocaltoGlobal;
    aMeshData.FieldsInfo = &tFieldsInfo;

    // Create mtk mesh
    return create_interpolation_mesh(MeshType::STK, aMeshData);
}

Mesh* generate_single_element_tet4()
{
	// Set up a simple mtk mesh
	uint tNumDim = 3;

	// Just a single Tet4
	uint aNumElemTypes = 1;
	Matrix< IdMat > tElemConnTet4 = {{1, 2, 3, 4}};
	Matrix< IdMat > tElemLocaltoGlobalTet4 = {{1}};

	Matrix< DDRMat > tCoords = {{0.0, 0.0, 0.0},
								{1.0, 0.0, 0.0},
								{0.0, 1.0, 0.0},
								{0.0, 0.0, 1.0}};

	Matrix< IdMat > tNodeLocaltoGlobal {{1, 2, 3, 4}};

	// Define block sets
	MtkSetsInfo tMtkMeshSets;

	Matrix< IdMat > tCellIdsBS({{1}});
	MtkBlockSetInfo tBlockSet;
	tBlockSet.mCellIdsInSet = &tCellIdsBS;
	tBlockSet.mBlockSetName = "blockset";
	tBlockSet.mBlockSetTopo = CellTopology::TET4;

	tMtkMeshSets.add_block_set(&tBlockSet);

	// Define side sets
	Matrix< IdMat > tElemIdsAndSideOrdsSS = {{1, 2},
											 {1, 3}};

	MtkSideSetInfo tSideSetStruct;
	tSideSetStruct.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS;
	tSideSetStruct.mSideSetName = "sideset";

	tMtkMeshSets.add_side_set(&tSideSetStruct);

	// Define node sets
	Matrix< IdMat > tNodeIdsNS = {{1}, {2}};

	MtkNodeSetInfo tNodeSet;
	tNodeSet.mNodeIds = &tNodeIdsNS;
	tNodeSet.mNodeSetName = "nodeset";

	tMtkMeshSets.add_node_set(&tNodeSet);

	// Create scalar field
	Scalar_Field_Info< DDRMat > tNodeField;

	Matrix < DDRMat > tNodeFieldData = {{0.0, 1.0, 0.5, 0.5}};

	tNodeField.set_field_name("nodefield");
	tNodeField.set_field_entity_rank(EntityRank::NODE);
	tNodeField.add_field_data(&tNodeLocaltoGlobal, &tNodeFieldData);

	// Collect field info
	MtkFieldsInfo tFieldsInfo;
	tFieldsInfo.mRealScalarFields.push_back(&tNodeField);

	// Coalesce mesh data
	MtkMeshData aMeshData(aNumElemTypes);
	aMeshData.SpatialDim = &tNumDim;
	aMeshData.ElemConn(0) = &tElemConnTet4;
	aMeshData.NodeCoords = &tCoords;
	aMeshData.SetsInfo = &tMtkMeshSets;
	aMeshData.LocaltoGlobalElemMap(0) = &tElemLocaltoGlobalTet4;
	aMeshData.LocaltoGlobalNodeMap = &tNodeLocaltoGlobal;
	aMeshData.FieldsInfo = &tFieldsInfo;

	// Create mtk mesh
	return create_interpolation_mesh(MeshType::STK, aMeshData);
}

Mesh* generate_single_element_quad4()
{
	// Set up a simple mtk mesh
	uint tNumDim = 2;

	// Just a single Quad4
	uint aNumElemTypes = 1;
	Matrix< IdMat > tElemConnQuad4 = {{1, 2, 3, 4}};
	Matrix< IdMat > tElemLocaltoGlobalQuad4 = {{1}};

	Matrix< DDRMat > tCoords = {{0.0, 0.0},
								{1.0, 0.0},
								{1.0, 1.0},
								{0.0, 1.0}};

	Matrix< IdMat > tNodeLocaltoGlobal {{1, 2, 3, 4}};

	// Define block sets
	MtkSetsInfo tMtkMeshSets;

	Matrix< IdMat > tCellIdsBS({{1}});
	MtkBlockSetInfo tBlockSet;
	tBlockSet.mCellIdsInSet = &tCellIdsBS;
	tBlockSet.mBlockSetName = "blockset";
	tBlockSet.mBlockSetTopo = CellTopology::QUAD4;

	tMtkMeshSets.add_block_set(&tBlockSet);

	// Define side sets
	Matrix< IdMat > tElemIdsAndSideOrdsSS = {{1, 2},
											 {1, 3}};

	MtkSideSetInfo tSideSetStruct;
	tSideSetStruct.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS;
	tSideSetStruct.mSideSetName = "sideset";

	tMtkMeshSets.add_side_set(&tSideSetStruct);

	// Define node sets
	Matrix< IdMat > tNodeIdsNS = {{1}, {2}};

	MtkNodeSetInfo tNodeSet;
	tNodeSet.mNodeIds = &tNodeIdsNS;
	tNodeSet.mNodeSetName = "nodeset";

	tMtkMeshSets.add_node_set(&tNodeSet);

	// Create scalar field
	Scalar_Field_Info< DDRMat > tNodeField;

	Matrix < DDRMat > tNodeFieldData = {{0.0, 1.0, 0.5, 0.5}};

	tNodeField.set_field_name("nodefield");
	tNodeField.set_field_entity_rank(EntityRank::NODE);
	tNodeField.add_field_data(&tNodeLocaltoGlobal, &tNodeFieldData);

	// Collect field info
	MtkFieldsInfo tFieldsInfo;
	tFieldsInfo.mRealScalarFields.push_back(&tNodeField);

	// Coalesce mesh data
	MtkMeshData aMeshData(aNumElemTypes);
	aMeshData.SpatialDim = &tNumDim;
	aMeshData.ElemConn(0) = &tElemConnQuad4;
	aMeshData.NodeCoords = &tCoords;
	aMeshData.SetsInfo = &tMtkMeshSets;
	aMeshData.LocaltoGlobalElemMap(0) = &tElemLocaltoGlobalQuad4;
	aMeshData.LocaltoGlobalNodeMap = &tNodeLocaltoGlobal;
	aMeshData.FieldsInfo = &tFieldsInfo;

	// Create mtk mesh
	return create_interpolation_mesh(MeshType::STK, aMeshData);
}

Mesh* generate_single_element_tri3()
{
	// Set up a simple mtk mesh
	uint tNumDim = 2;

	// Just a single Quad4
	uint aNumElemTypes = 1;
	Matrix< IdMat > tElemConnQuad4 = {{1, 2, 3}};
	Matrix< IdMat > tElemLocaltoGlobalQuad4 = {{1}};

	Matrix< DDRMat > tCoords = {{0.0, 0.0},
								{1.0, 0.0},
								{0.0, 1.0}};

	Matrix< IdMat > tNodeLocaltoGlobal {{1, 2, 3}};

	// Define block sets
	MtkSetsInfo tMtkMeshSets;

	Matrix< IdMat > tCellIdsBS({{1}});
	MtkBlockSetInfo tBlockSet;
	tBlockSet.mCellIdsInSet = &tCellIdsBS;
	tBlockSet.mBlockSetName = "blockset";
	tBlockSet.mBlockSetTopo = CellTopology::TRI3;

	tMtkMeshSets.add_block_set(&tBlockSet);

	// Define side sets
	Matrix< IdMat > tElemIdsAndSideOrdsSS = {{1, 2},
											 {1, 1}};

	MtkSideSetInfo tSideSetStruct;
	tSideSetStruct.mElemIdsAndSideOrds = &tElemIdsAndSideOrdsSS;
	tSideSetStruct.mSideSetName = "sideset";

	tMtkMeshSets.add_side_set(&tSideSetStruct);

	// Define node sets
	Matrix< IdMat > tNodeIdsNS = {{1}, {2}};

	MtkNodeSetInfo tNodeSet;
	tNodeSet.mNodeIds = &tNodeIdsNS;
	tNodeSet.mNodeSetName = "nodeset";

	tMtkMeshSets.add_node_set(&tNodeSet);

	// Create scalar field
	Scalar_Field_Info< DDRMat > tNodeField;

	Matrix < DDRMat > tNodeFieldData = {{0.0, 1.0, 0.5}};

	tNodeField.set_field_name("nodefield");
	tNodeField.set_field_entity_rank(EntityRank::NODE);
	tNodeField.add_field_data(&tNodeLocaltoGlobal, &tNodeFieldData);

	// Collect field info
	MtkFieldsInfo tFieldsInfo;
	tFieldsInfo.mRealScalarFields.push_back(&tNodeField);

	// Coalesce mesh data
	MtkMeshData aMeshData(aNumElemTypes);
	aMeshData.SpatialDim = &tNumDim;
	aMeshData.ElemConn(0) = &tElemConnQuad4;
	aMeshData.NodeCoords = &tCoords;
	aMeshData.SetsInfo = &tMtkMeshSets;
	aMeshData.LocaltoGlobalElemMap(0) = &tElemLocaltoGlobalQuad4;
	aMeshData.LocaltoGlobalNodeMap = &tNodeLocaltoGlobal;
	aMeshData.FieldsInfo = &tFieldsInfo;

	// Create mtk mesh
	return create_interpolation_mesh(MeshType::STK, aMeshData);
}

void check_normals( Mesh* aMesh, const std::shared_ptr< Cell_Info >& aConnectivity, Matrix< DDRMat >* aExpectedNormals, uint aNumFaces )
{
	// Access first cell
	Cell & tCell = aMesh->get_mtk_cell(0);
	Matrix< IndexMat > tCellNodes = tCell.get_vertex_inds();

	// Loop through all faces and check normals
	for( uint i = 0; i < aNumFaces; i++)
	{
		// Get face normal
		Matrix< DDRMat > tFaceNormal = tCell.compute_outward_side_normal(i);

		// Check face normal against expected vectors
		CHECK( all_true( aExpectedNormals[i] == tFaceNormal ) );
	}
}

void check_normals_2D( Mesh* aMesh, const std::shared_ptr< Cell_Info >& aConnectivity, Matrix< DDRMat >* aExpectedNormals, uint aNumEdges )
{
	// Access first cell
	Cell & tCell = aMesh->get_mtk_cell(0);
	Matrix< IndexMat > tCellNodes = tCell.get_vertex_inds();

	// Loop through all faces and check normals
	for( uint i = 0; i < aNumEdges; i++)
	{
		// Get node ordinals to determine the edge normal with
		Matrix< IndexMat > tNodeOrdinals = aConnectivity->get_node_map_outward_normal(i);

		// Get global coordinates for each node
		Matrix< DDRMat > tNodeCoords0 = aMesh->get_node_coordinate(tCellNodes(tNodeOrdinals(0)));
		Matrix< DDRMat > tNodeCoords1 = aMesh->get_node_coordinate(tCellNodes(tNodeOrdinals(1)));

		// Calculate vector between relevant nodes
		Matrix< DDRMat > tVec = tNodeCoords1 - tNodeCoords0;

		// Convert {x, y} to {-y, x} to determine right-handed normal
		Matrix< DDRMat > tFaceNormal = {{-tVec(1), tVec(0)}};

		// Check edge normal against expected vectors
		CHECK( all_true( aExpectedNormals[i] == tFaceNormal ) );
	}
}

void check_faces( Mesh* aMesh, const std::shared_ptr< Cell_Info >& aConnectivity, uint aNumFaces, uint aNumNodesPerFace )
{
	// Access first cell
	Cell & tCell = aMesh->get_mtk_cell(0);
	Matrix< IndexMat > tCellNodes = tCell.get_vertex_inds();

	// Get local indices of the faces on the element
	Matrix< IndexMat > tElementFaces = aMesh->get_faces_connected_to_element_loc_inds(0);

	// Get the matrix containing the full node to face ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllFaces = aConnectivity->get_node_to_face_map();

	// Loop through all faces
	for( uint i = 0; i < aNumFaces; i++ )
	{
		// Get the local indices of the nodes on the current face from the mesh
		Matrix< IndexMat > tExpectedFaceNodeInds = aMesh->get_entity_connected_to_entity_loc_inds(tElementFaces(i), EntityRank::FACE, EntityRank::NODE);

		// Get the node ordinal map for the current face
		Matrix< IndexMat > tNodeOrdinalsOnSingleFace = aConnectivity->get_node_to_face_map(i);

		// Change node ordinals from both functions to local indices
		Matrix< IndexMat > tFaceNodeInds(1,aNumNodesPerFace);
		Matrix< IndexMat > tFaceNodeIndsFromAll(1,aNumNodesPerFace);
		for( uint j = 0; j < aNumNodesPerFace; j++ )
		{
			tFaceNodeInds(j) = tCellNodes(tNodeOrdinalsOnSingleFace(j));
			tFaceNodeIndsFromAll(j) = tCellNodes(tNodeOrdinalsOnAllFaces(i,j));
		}

		// Sort local index vectors in ascending order
		sort(tExpectedFaceNodeInds, tExpectedFaceNodeInds, "ascend", 1);
		sort(tFaceNodeInds, tFaceNodeInds, "ascend", 1);
		sort(tFaceNodeIndsFromAll, tFaceNodeIndsFromAll, "ascend", 1);

		// Check local indices from ordinal map against the expected local indices from the mesh
		CHECK( all_true( tExpectedFaceNodeInds == tFaceNodeInds ) );
		CHECK( all_true( tExpectedFaceNodeInds == tFaceNodeIndsFromAll ) );
	}
}

void check_edges( Mesh* aMesh, const std::shared_ptr< Cell_Info >& aConnectivity, uint aNumEdges )
{
	// Access first cell
	Cell & tCell = aMesh->get_mtk_cell(0);
	Matrix< IndexMat > tCellNodes = tCell.get_vertex_inds();

	// Get local indices of the edges on the element
	Matrix< IndexMat > tElementEdges = aMesh->get_edges_connected_to_element_loc_inds(0);

	// Get the matrix containing the full node to edge ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllEdges = aConnectivity->get_node_to_edge_map();

	// Loop through all twelve edges
	for( uint i = 0; i < aNumEdges; i++ )
	{
		// Get the local indices of the nodes on the current edge from the mesh
		Matrix< IndexMat > tExpectedEdgeNodeInds = aMesh->get_entity_connected_to_entity_loc_inds(tElementEdges(i), EntityRank::EDGE, EntityRank::NODE);

		// Get the node ordinal map for the current edge
		Matrix< IndexMat > tNodeOrdinalsOnSingleEdge = aConnectivity->get_node_to_edge_map(i);

		// Change node ordinals from both functions to local indices
		Matrix< IndexMat > tEdgeNodeInds(1,2);
		Matrix< IndexMat > tEdgeNodeIndsFromAll(1,2);
		for( uint j = 0; j < 2; j++ )
		{
			tEdgeNodeInds(j) = tCellNodes(tNodeOrdinalsOnSingleEdge(j));
			tEdgeNodeIndsFromAll(j) = tCellNodes(tNodeOrdinalsOnAllEdges(i,j));
		}

		// Sort local index vectors in ascending order
		sort(tExpectedEdgeNodeInds, tExpectedEdgeNodeInds, "ascend", 1);
		sort(tEdgeNodeInds, tEdgeNodeInds, "ascend", 1);
		sort(tEdgeNodeIndsFromAll, tEdgeNodeIndsFromAll, "ascend", 1);

		// Check local indices from ordinal map against the expected local indices from the mesh
		CHECK( all_true( tExpectedEdgeNodeInds == tEdgeNodeInds ) );
		CHECK( all_true( tExpectedEdgeNodeInds == tEdgeNodeIndsFromAll ) );
	}
}

void check_facets( const std::shared_ptr< Cell_Info >& aConnectivity, uint aNumFacets )
{
	// Get the matrix containing the full node to face ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllFaces = aConnectivity->get_node_to_face_map();

	// Get the matrix containing the full node to facet ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllFacets = aConnectivity->get_node_to_facet_map();

	// Check that the facet map is the same as the face map for 3D elements
	CHECK( all_true( tNodeOrdinalsOnAllFaces == tNodeOrdinalsOnAllFacets ) );

	// Loop through all facets (faces)
	for( uint i = 0; i < aNumFacets; i++ )
	{
		// Get node ordinal map on current facet and current face
		Matrix< IndexMat > tNodeOrdinalsOnSingleFacet = aConnectivity->get_node_to_facet_map(i);
		Matrix< IndexMat > tNodeOrdinalsOnSingleFace = aConnectivity->get_node_to_face_map(i);

		// Check that the facet map is the same as the face map for 3D elements
		CHECK( all_true( tNodeOrdinalsOnSingleFace == tNodeOrdinalsOnSingleFacet ) );
	}
}

void check_facets_2D( const std::shared_ptr< Cell_Info >& aConnectivity, uint aNumFacets )
{
	// Get the matrix containing the full node to edge ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllEdges = aConnectivity->get_node_to_edge_map();

	// Get the matrix containing the full node to facet ordinal map
	Matrix< IndexMat > tNodeOrdinalsOnAllFacets = aConnectivity->get_node_to_facet_map();

	// Check that the facet map is the same as the edge map for 2D elements
	CHECK( all_true( tNodeOrdinalsOnAllEdges == tNodeOrdinalsOnAllFacets ) );

	// Loop through all facets (edges)
	for( uint i = 0; i < aNumFacets; i++ )
	{
		// Get node ordinal map on current facet and current edge
		Matrix< IndexMat > tNodeOrdinalsOnSingleFacet = aConnectivity->get_node_to_facet_map(i);
		Matrix< IndexMat > tNodeOrdinalsOnSingleEdge = aConnectivity->get_node_to_edge_map(i);

		// Check that the facet map is the same as the face map for 2D elements
		CHECK( all_true( tNodeOrdinalsOnSingleEdge == tNodeOrdinalsOnSingleFacet ) );
	}
}

TEST_CASE("Hex8 Connectivity Test", "[MTK_CONN_HEX8]")
{
	// Create test mesh
	Mesh* aMesh = generate_single_element_hex8();

	// Get a hex8 connectivity object
	std::shared_ptr<Cell_Info> aConnectivity = std::make_shared<Cell_Info_Hex8>();

	// Test outward normal node map
	// Expected normal vectors
	Matrix< DDRMat > aExpectedNormals[] = { {{ 0.0}, {-1.0}, { 0.0}},
						{{ 1.0}, { 0.0}, { 0.0}},
						{{ 0.0}, { 1.0}, { 0.0}},
						{{-1.0}, { 0.0}, { 0.0}},
						{{ 0.0}, { 0.0}, {-1.0}},
						{{ 0.0}, { 0.0}, { 1.0}}};

	check_normals(aMesh, aConnectivity, aExpectedNormals, 6);

	// Test node to face map
	check_faces(aMesh, aConnectivity, 6, 4);

	// Test node to edge map
	check_edges(aMesh, aConnectivity, 12);

	// Test node to facet map
	check_facets(aConnectivity, 6);

	// Delete used mesh
	delete aMesh;
}

TEST_CASE("Tet4 Connectivity Test", "[MTK_CONN_TET4]")
{
	// Create test mesh
	Mesh* aMesh = generate_single_element_tet4();

	// Get a tet4 connectivity object
	std::shared_ptr<Cell_Info> aConnectivity = std::make_shared<Cell_Info_Tet4>();

	// Test outward normal node map
	// Expected normal vectors
	Matrix< DDRMat > aExpectedNormals[] = {{{ 0.0 },                 { -1.0 },                  { 0.0 }},
					       {{+5.773502691896258e-01},{ +5.773502691896258e-01 },{ +5.773502691896258e-01 }},
					       {{-1.0 },{  0.0 },{ 0.0 }},
					       {{ 0.0 },{  0.0 },{-1.0 }}};

	check_normals(aMesh, aConnectivity, aExpectedNormals, 4);

	// Test node to face map
	check_faces(aMesh, aConnectivity, 4, 3);

	// Test node to edge map
	check_edges(aMesh, aConnectivity, 6);

	// Test node to facet map
	check_facets(aConnectivity, 4);

	// Delete used mesh
	delete aMesh;
}

TEST_CASE("Quad4 Connectivity Test", "[MTK_CONN_QUAD4]")
{
	// Create test mesh
	Mesh* aMesh = generate_single_element_quad4();

	// Get a quad4 connectivity object
	std::shared_ptr<Cell_Info> aConnectivity = std::make_shared<Cell_Info_Quad4>();

	// Test outward normal node map
	// Expected normal vectors
	Matrix< DDRMat > aExpectedNormals[] = {{{ 0.0, -1.0}},
									       {{ 1.0,  0.0}},
									       {{ 0.0,  1.0}},
									       {{-1.0,  0.0}}};

//	check_normals_2D(aMesh, aConnectivity, aExpectedNormals, 4);

	// Test node to edge map
	check_edges(aMesh, aConnectivity, 4);

	// Test node to facet map
	check_facets_2D(aConnectivity, 4);

	// Delete used mesh
	delete aMesh;
}

TEST_CASE("Tri3 Connectivity Test", "[MTK_CONN_TRI3]")
{
	// Create test mesh
	Mesh* aMesh = generate_single_element_tri3();

	// Get a quad4 connectivity object
	std::shared_ptr<Cell_Info> aConnectivity = std::make_shared<Cell_Info_Tri3>();

	// Test outward normal node map
	// Expected normal vectors
	Matrix< DDRMat > aExpectedNormals[] = {{{ 0.0, -1.0}},
									       {{ 1.0,  1.0}},
									       {{-1.0,  0.0}}};

	check_normals_2D(aMesh, aConnectivity, aExpectedNormals, 3);

	// Test node to edge map
	check_edges(aMesh, aConnectivity, 3);

	// Test node to facet map
	check_facets_2D(aConnectivity, 3);

	// Delete used mesh
	delete aMesh;
}

}
