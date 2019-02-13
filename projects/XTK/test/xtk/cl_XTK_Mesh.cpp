/*
 * cl_XTK_Background_Mesh.hpp
 *
 *  Created on: Jul 21, 2017
 *      Author: ktdoble
 */

#include <utility>
#include <iostream>

#include "cl_XTK_Background_Mesh.hpp"
#include "catch.hpp"

#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_Sphere.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Mesh.hpp"

#include "mesh/cl_Mesh_Data.hpp"
#include "cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"


#include "cl_Cell.hpp"


namespace xtk
{

TEST_CASE("XTK Mesh with Element Downward Inheritance",
          "[Background_Mesh]""[Downward_Inheritance]")
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    /*
     * This section tests the following:
     * 1.) Hardcoding of the element to simple mesh pairs in serial
     */
    SECTION("Hardcoded downward inheritance pairs")
    {
        if(tProcSize==1)
        {
        /*
         * Setup Mesh
         */
        std::string tMeshFileName = "generated:5x5x5";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);


        /*
         * Construct XTK Mesh
         */
        Background_Mesh<real, size_t, moris::DDRMat, moris::DDSTMat> tXTKMesh(tMeshData);

        /*
         * Setup Cut Mesh
         */
        Cut_Mesh tCutMesh(3);

        /*
         * Say Element with index 18 is paired up with child mesh index 0;
         */
        moris::Cell<std::pair<size_t,size_t>> tElementToSimpleMeshPairs(2);
        tElementToSimpleMeshPairs(0) = std::pair<size_t,size_t>(18,0);
        tElementToSimpleMeshPairs(1) = std::pair<size_t,size_t>(14,24);

        /*
         * Register the Pair
         */
        tXTKMesh.register_new_downward_inheritance(tElementToSimpleMeshPairs);

        /*
         * Ask if elements  18, 14 and 23 have children
         * 18 and 14 should have children, 23 should not
         */
        CHECK(tXTKMesh.entity_has_children(18,EntityRank::ELEMENT));
        CHECK(tXTKMesh.child_mesh_index(18,EntityRank::ELEMENT) == 0);
        CHECK(tXTKMesh.entity_has_children(14,EntityRank::ELEMENT));
        CHECK(tXTKMesh.child_mesh_index(14,EntityRank::ELEMENT) == 24);
        CHECK(!tXTKMesh.entity_has_children(23,EntityRank::ELEMENT));
        }
    }

    /*
     * This section tests the following:
     * 1.) XTK decomposition determines the inheritance pairs of the element to simple mesh pairs
     */
    SECTION("XTK downward inheritance pairs")
    {
        /*
         * Use a level set sphere to intersect mesh,
         * The sphere is intentionally set to intersect 2 of the 3 elements
         */
        real tRadius = 0.25;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 1;
        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

        /*
         * Setup Mesh
         */
        std::string tMeshFileName = "generated:1x1x3";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);


        /*
         * Decompose
         * Note: this establishes the link that is hardcoded in the above pair
         */
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
        tXTKModel.decompose(tDecompositionMethods);

        // Access the XTK Mesh-------------------------------------------------------------
        Background_Mesh<real,size_t, moris::DDRMat, moris::DDSTMat> & tXTKMesh = tXTKModel.get_background_mesh();

        CHECK(tXTKMesh.entity_has_children(0,EntityRank::ELEMENT));
        CHECK(tXTKMesh.child_mesh_index(0,EntityRank::ELEMENT)==0);
        CHECK(tXTKMesh.entity_has_children(1,EntityRank::ELEMENT));
        CHECK(tXTKMesh.child_mesh_index(1,EntityRank::ELEMENT)==1);
        CHECK(!tXTKMesh.entity_has_children(2,EntityRank::ELEMENT));

    }

}

}
