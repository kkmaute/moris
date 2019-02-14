/*
 * cl_XTK_Child_Mesh_RegSub_2D.cpp
 *
 *  Created on: Jan 25, 2019
 *      Author: ryan
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "op_times.hpp"
#include "op_minus.hpp"

#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "xtk/cl_XTK_Child_Mesh.hpp"


using namespace moris;
namespace xtk
{
TEST_CASE("Regular Subdivision QUAD4","[REG_SUB_TEMPLATE_QUAD4]")
{
	// Spatial dimension
	uint tSpatialDim = 2;

	// Number of nodes in template
	uint tNumNodesInTemplate = 5;

	// Node Coordinates
	Matrix<DDRMat> tNodeCoords(tNumNodeInTemplate, tSpatialDim);
	tNodeCoords(0,0) = 0.0; tNodeCoords(0,1) = 0.0;
	tNodeCoords(1,0) = 1.0; tNodeCoords(1,1) = 0.0;
	tNodeCoords(2,0) = 0.0; tNodeCoords(2,1) = 1.0;
	tNodeCoords(3,0) = 1.0; tNodeCoords(3,1) = 1.0;
	tNodeCoords(4,0) = 0.5; tNodeCoords(4,1) = 0.5;

	// All nodes in template
    Matrix< moris::IndexMat > tNodeIndex({{0, 1, 3, 2, 4}});

    // All node ids in template
    Matrix< moris::IndexMat > tNodeIds({{1,2,4,3,5}});

    // Add ancestry information (on board we called this template parent entity parents
    // Do for nodes, edges, cells


//     Setup mesh modification template
//     Initialize Template
//    Mesh_Modification_Template tRegSubTemplate(tElementsAncestry(0,0),
//                                               0,
//                                               tNodeIndex,
//                                               tParentNodeInds,
//                                               tParentNodeRanks,
//                                               tParentEdgeInds,
//                                               tParentEdgeRanks,
//                                               tParentFaceInds,
//                                               tParentFaceRanks,
//                                                Change enum TemplateType::REGULAR_SUBDIVISION_HEX8);



}
}
