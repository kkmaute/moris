


/*
 * cl_XTK_Model_Recursion.hpp
 *
 *  Created on: Oct 20, 2017
 *      Author: ktdoble
 */

#include <memory>
#include <mpi.h>
#include "catch.hpp"
#include "cl_Logger.hpp"


#include "xtk/cl_XTK_Child_Mesh.hpp"


#include "cl_Sphere.hpp"
#include "geometry/cl_Plane.hpp"
// XTKL: Mesh Includes
#include "cl_Mesh_Enums.hpp"
#include "fn_verify_tet_topology.hpp"

// XTKL: Geometry  Include
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

#include "cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Phase_Table.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Face_Registry.hpp"

#include "fn_local_child_mesh_flood_fill.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_generate_element_to_element.hpp"

using namespace moris;
namespace xtk
{
TEST_CASE("Phase Table","[Phase_Table]")
        {

    Matrix< IndexMat > tPhaseTableData (
            {{0,0},
        {0,1},
        {1,0},
        {1,1}});

    Cell<std::string> tPhaseNames = {"m0","m1","m2"};

    // TODO: figure out how to check a throw in this constructor
    //    CHECK_THROWS(Phase_Table(tPhaseTableData,tPhaseNames));

    tPhaseNames = {"m0","m1","m2","m3"};

    Phase_Table tPhaseTable (tPhaseTableData,tPhaseNames);
    Matrix< IndexMat > tRow(0,0);

    size_t tIndex = 0;
    for(size_t iR = 0; iR<tPhaseTableData.n_rows(); iR++ )
    {
        tRow = tPhaseTableData.get_row(iR);

        tIndex = tPhaseTable.get_phase_index(tRow);

        CHECK(tIndex == iR);
    }
//
//#ifdef DEBUG
//    (tRow)(0,0) = 2;
//    CHECK_THROWS(tPhaseTable.get_phase_index(tRow));
//    (tRow)(0,0) = 0;
//    (tRow)(0,1) = 2;
//
//    CHECK_THROWS(tPhaseTable.get_phase_index(tRow));
//#endif


    // Check a 3 phase problem

    tPhaseTableData = Matrix< IndexMat >(
            {{0,0,0},
        {0,0,1},
        {0,1,0},
        {0,1,1},
        {1,0,0},
        {1,0,1},
        {1,1,0},
        {1,1,1}});


    // TODO: figure out how to check a throw in this constructor
    //    CHECK_THROWS(Phase_Table(tPhaseTableData,tPhaseNames));

    tPhaseNames = {"m0","m1","m2","m3","m4","m5","m6","m7"};

    Phase_Table tPhaseTable2 (tPhaseTableData,tPhaseNames);

    //Check indices are correct
    tRow = Matrix< IndexMat >(1,tPhaseTableData.n_cols());
    for(size_t iR = 0; iR<tPhaseTableData.n_rows(); iR++ )
    {
        tRow = tPhaseTableData.get_row(iR);

        tIndex = tPhaseTable2.get_phase_index(tRow);

        CHECK(tIndex == iR);
    }

#ifdef DEBUG
    (tRow)(0,0) = 2;

    CHECK_THROWS(tPhaseTable2.get_phase_index(tRow));
    (tRow)(0,0) = 0;
    (tRow)(0,1) = 2;

    CHECK_THROWS(tPhaseTable2.get_phase_index(tRow));
#endif


        }


TEST_CASE("Autogenerate Exponential Base 2 Table","[AUTO_PHASE_TABLE]")
{
    Phase_Table tPhaseTable (3, Phase_Table_Structure::EXP_BASE_2);


    Matrix< IndexMat > tExpectedPhaseTableData (
            {{0,0,0},
        {0,0,1},
        {0,1,0},
        {0,1,1},
        {1,0,0},
        {1,0,1},
        {1,1,0},
        {1,1,1}});

    CHECK(all_true(tExpectedPhaseTableData == tPhaseTable.get_phase_table_data()));
}
TEST_CASE("2 Nonintersecting geometries","[2_Phase],[NO_OVER]")
{

    // FIXME: THIS EXAMPLE IS NOT HANDLING AN MPI CALL CORRECTLY
    //    Matrix< IndexMat > tPhaseTableData (
    //            {{0,0},
    //             {0,1},
    //             {1,0},
    //             {1,1}});
    //
    //    Cell<std::string> tPhaseNames = {"m0","m1","m2","m3"};
    //
    //    // Geometry Engine Setup -----------------------
    //    // Using a 2 Levelset Spheres as the Geometry
    //    real tRadius1  = 2.70;
    //    real tXCenter1 = 4.51;
    //    real tYCenter1 = 4.51;
    //    real tZCenter1 = 4.51;
    //    Sphere tLevelsetSphere1(tRadius1, tXCenter1, tYCenter1, tZCenter1);
    //
    //    real tRadius2  = 2.80;
    //    real tXCenter2 = 4.51;
    //    real tYCenter2 = 4.51;
    //    real tZCenter2 = 5.51;
    //    Sphere tLevelsetSphere2(tRadius2, tXCenter2, tYCenter2, tZCenter2);
    //
    //
    //    Phase_Table tPhaseTable (2,  Phase_Table_Structure::EXP_BASE_2);
    //    Geometry_Engine tGeometryEngine(
    //        {&tLevelsetSphere1,
    //         &tLevelsetSphere2},
    //        tPhaseTable);
    //
    //    // Create Mesh ---------------------------------
    //    std::string tMeshFileName = "generated:10x10x10";
    //    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName );
    //
    //    // Setup XTK Model -----------------------------
    //    size_t tModelDimension = 3;
    //    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    //
    //    tXTKModel.mVerbose=true;
    //
    //    //Specify your decomposition methods and start cutting
    //    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    //    tXTKModel.decompose(tDecompositionMethods);
    //
    //    Output_Options tOutputOptions;
    //
    //    // Write output mesh
    //    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
    //
    //    std::string tPrefix = std::getenv("MORISOUTPUT");
    //    std::string tMeshOutputFile = tPrefix + "/unit_recursive_no_intersect.e";
    //    tCutMeshData->create_output_mesh(tMeshOutputFile);
    //
    //
    //    delete tMeshData;
    //    delete tCutMeshData;


}


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

        Plane tPlane1(tXc1,tYc1,tZc1,tXn1,tYn1,tZn1);


        real tXc2 = 0.55;
        real tYc2 = 0.55;
        real tZc2 = 0.55;
        real tXn2 = 0.0;
        real tYn2 = 0.0;
        real tZn2 = 1.0;

        Plane tPlane2(tXc2,tYc2,tZc2,tXn2,tYn2,tZn2);

        real tXc3 = 0.7;
        real tYc3 = 0.7;
        real tZc3 = 0.7;
        real tXn3 = 0.0;
        real tYn3 = 0.0;
        real tZn3 = 1.0;

        Plane tPlane3(tXc3,tYc3,tZc3,tXn3,tYn3,tZn3);


        Phase_Table tPhaseTable (3,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine({&tPlane1, &tPlane2,&tPlane3},tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:2x2x2";
        Cell<std::string> tScalarFields(0);
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

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
