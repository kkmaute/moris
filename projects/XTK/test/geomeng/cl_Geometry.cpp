/*
 * cl_Geometry.cpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */



#include "geometry/cl_Geometry.hpp"

#include "geometry/cl_Composite_Fiber.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geometry/cl_Sphere.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp" // For print
#include "catch.hpp"
#include "mesh/cl_Mesh_Builder.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"

// XTK Includes
#include "xtk/cl_XTK_Pending_Node.hpp"

// Topology
#include "topology/cl_XTK_Topology.hpp"
#include "topology/cl_XTK_Edge_Topology.hpp"
#include "topology/cl_XTK_Quad_4_Topology.hpp"
#include "topology/cl_XTK_Hexahedron_8_Topology.hpp"

#include "linalg_typedefs.hpp"

namespace xtk
{
TEST_CASE("Discretized Level Set with 2 spheres","[DISCRETE_LEVELSET][SPHERE]")
{

    // Create Mesh----------------------------------
    std::string tMeshInput = "generated:1x1x1";
    Cell<std::string> tScalarFieldNames = {"LEVELSET_01","LEVELSET_02"};
    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
    // Setup geometry and discretize onto a levelset mesh--------------------
    //TODO: LOAD AN EXODUS FILE INSTEAD OF USING ANALYTIC LEVELSET EXPRESSION
    real tRadius =  0.25;
    real tXCenter = 1;
    real tYCenter = 1;
    real tZCenter = 1;
    Sphere tSphere1(tRadius, tXCenter, tYCenter, tZCenter);
    tRadius =  0.25;
    tXCenter = 0;
    tYCenter = 0;
    tZCenter = 0;
    Sphere tSphere2(tRadius, tXCenter, tYCenter, tZCenter);
    // ------------------------------------------------------------------------
    Cell<Geometry<real, size_t, moris::DDRMat, moris::DDSTMat>*> tFunctionsToDiscretize = {&tSphere1,&tSphere2};
    Discrete_Level_Set<real,size_t, moris::DDRMat, moris::DDSTMat> tLevelSetMesh(tFunctionsToDiscretize,tMeshInput,tScalarFieldNames,tMeshBuilder);
    std::shared_ptr<mesh::Mesh_Data<real,size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tLevelSetMesh.get_level_set_mesh();
    // ------------------------------------------------------------------------
    CHECK(!tLevelSetMesh.is_analytic());
    CHECK(tLevelSetMesh.get_active_level_set_field_name()=="LEVELSET_01");
    CHECK(tLevelSetMesh.get_level_set_field_name(0) == "LEVELSET_01");
    CHECK(tLevelSetMesh.get_level_set_field_name(1) == "LEVELSET_02");
    CHECK(tLevelSetMesh.get_num_levelset() == 2);

//    CHECK_THROWS(tLevelSetMesh.get_level_set_field_name(2));

    // These may fail if the mesh implementation uses a different node numbering
    // Check to see that the negative level set value corresponds to the node with
    // coordinates (1,1,1).
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(0,EntityRank::NODE) == Approx(2.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(1,EntityRank::NODE) == Approx(1.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(2,EntityRank::NODE) == Approx(1.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(3,EntityRank::NODE) == Approx(0.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(4,EntityRank::NODE) == Approx(1.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(5,EntityRank::NODE) == Approx(0.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(6,EntityRank::NODE) == Approx(0.9375));
    CHECK(tLevelSetMesh.access_field_value_with_entity_index(7,EntityRank::NODE) == Approx(-0.0625));



}

}
