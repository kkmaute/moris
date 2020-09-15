


/*
 * cl_XTK_Model_Recursion.hpp
 *
 *  Created on: Oct 20, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_Plane.hpp"
using namespace moris;
namespace xtk
{

TEST_CASE("2 Intersecting Geometries","[2_Phase][OVER]")
{

    if(par_size() == 1 || par_size() ==2)
    {
        Matrix< IndexMat > tPhaseTableData (
                {{0,0},
            {0,1},
            {1,0},
            {1,1}});

        Cell<std::string> tPhaseNames = {"m0","m1","m2","m3"};

        // Geometry Engine Setup -----------------------
        real tXc1 = 0.4;
        real tYc1 = 0.4;
        real tZc1 = 0.4;
        real tXn1 = 0.0;
        real tYn1 = 1.0;
        real tZn1 = 0.0;

        Plane<3> tPlane1(tXc1,tYc1,tZc1,tXn1,tYn1,tZn1);


        real tXc2 = 0.55;
        real tYc2 = 0.55;
        real tZc2 = 0.55;
        real tXn2 = 0.0;
        real tYn2 = 0.0;
        real tZn2 = 1.0;

        Plane<3> tPlane2(tXc2,tYc2,tZc2,tXn2,tYn2,tZn2);

        real tXc3 = 0.7;
        real tYc3 = 0.7;
        real tZc3 = 0.7;
        real tXn3 = 0.0;
        real tYn3 = 0.0;
        real tZn3 = 1.0;

        Plane<3> tPlane3(tXc3,tYc3,tZc3,tXn3,tYn3,tZn3);


        Phase_Table tPhaseTable (3);
        Geometry_Engine tGeometryEngine({&tPlane1, &tPlane2,&tPlane3},tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:2x2x2";
        Cell<std::string> tScalarFields(0);
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose  =  false;

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);


        // Verify that the tets created have correct topology
        Cut_Mesh   & tCutMesh   = tXTKModel.get_cut_mesh();
        Child_Mesh & tChildMesh = tCutMesh.get_child_mesh(0);

        bool tValidTopo = verify_tet4_topology(tChildMesh.get_element_to_node(),
                                               tChildMesh.get_element_to_edge(),
                                               tChildMesh.get_element_to_face(),
                                               tChildMesh.get_edge_to_node(),
                                               tChildMesh.get_face_to_node());

        CHECK(tValidTopo);


        Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = false;
        tOutputOptions.mInternalUseFlag = true;
        mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

        std::string tPrefix = std::getenv("MORISOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/unit_recursive_intersect.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);


        delete tMeshData;
        delete tCutMeshData;


    }
}
}
