/*
 * UT_XTK_Model_2D.cpp
 *
 *  Created on: Aug 7, 2019
 *      Author: ryan
 */


#include "catch.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_XTK_Model.hpp"

#include "cl_Circle.hpp"
#include "cl_MGE_Geometry_Engine.hpp"

using namespace moris;

namespace xtk
{

TEST_CASE("2D Regular Subdivision Method","[REGULAR_SUBDIVISION_MODEL_2D]")
        {
    int tProcRank = par_rank();
    int tProcSize = par_size();

    if(tProcSize==1)
    {
        // Geometry Engine Setup -----------------------
        // Using a Levelset Circle as the Geometry

        real tRadius = 0.7;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        Circle tLevelsetCircle(tRadius, tXCenter, tYCenter);

        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelsetCircle,tPhaseTable,2);

        // Create Mesh ---------------------------------
        // Generate data for test
        uint aNumDim = 2;
        Matrix< DDRMat >  aCoords(6,2);
        aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
        aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
        aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
        aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
        aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
        aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
        Matrix< IdMat >     aElemConn( 2, 4 );

        // 0D to 3D connectivity (node to element)
        aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
        aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

        Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

        // No need of an element map since elements in connectivity table are assumed to be contiguous

        // Create MORIS mesh using MTK database
        moris::mtk::MtkMeshData aMeshData;
        aMeshData.CreateAllEdgesAndFaces = true;
        aMeshData.SpatialDim = &aNumDim;
        aMeshData.ElemConn(0)= &aElemConn;
        aMeshData.NodeCoords = &aCoords;
        aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, aMeshData );

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 2;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
        tXTKModel.decompose(tDecompositionMethods);

        // Access the decomposed XTK Mesh
        Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

        // Do some testing
        size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
        size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

        // 2 element was subdivided
        std::cout<<"tNumNodesAfterDecompositionXTK = "<<tNumNodesAfterDecompositionXTK<<std::endl;

        CHECK(tNumNodesAfterDecompositionXTK == 10); /* two duplicates from the shared nodes*/
        CHECK(tNumElementsAfterDecompositionXTK == 8);

        moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

        moris::Matrix< moris::DDRMat > tExpectedNodeCoordinates = {{0.0, 0.0},
                                                                   {1.0, 0.0},
                                                                   {1.0, 1.0},
                                                                   {0.0, 1.0},
                                                                   {2.0, 0.0},
                                                                   {2.0, 1.0},
                                                                   {0.5, 0.5},
                                                                   {1.5, 0.5}};
        CHECK(equal_to(tNodeCoordinates,tExpectedNodeCoordinates));



        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
        std::string tMeshOutputFile    = "./xtk_exo/xtk_test_output_regular_subdivision_2d.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);
        delete tMeshData;
        delete tCutMeshData;

    }
        }

TEST_CASE("2D Conformal Subdivision","[CONFORMAL_MODEL_2D]")
{
    if(par_size()==1)
    {
        // Geometry Engine Setup -----------------------
        // Using a Levelset Circle as the Geometry
        real tRadius = 0.7;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        Circle tLevelsetCircle(tRadius, tXCenter, tYCenter);

        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tLevelsetCircle,tPhaseTable,2);

        // Create Mesh ---------------------------------
        // Generate data for test
        uint aNumDim = 2;
        Matrix< DDRMat >  aCoords(6,2);
        aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
        aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
        aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
        aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
        aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
        aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
        Matrix< IdMat >     aElemConn( 2, 4 );

        // 0D to 3D connectivity (node to element)
        aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
        aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

        Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

        // No need of an element map since elements in connectivity table are assumed to be contiguous

        // Create MORIS mesh using MTK database
        moris::mtk::MtkMeshData aMeshData;
        aMeshData.CreateAllEdgesAndFaces = true;
        aMeshData.SpatialDim = &aNumDim;
        aMeshData.ElemConn(0)= &aElemConn;
        aMeshData.NodeCoords = &aCoords;
        aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, aMeshData );

        // Setup XTK Model ----------------------------------------------------------------
        size_t tModelDimension = 2;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify decomposition Method and Cut Mesh ---------------------------------------
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
        tXTKModel.decompose(tDecompositionMethods);

        // output to exodus file ----------------------------------------------------------
        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

        std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_conformal_2d.e";
        tCutMeshData->create_output_mesh(tMeshOutputFile);

        delete tCutMeshData;
        delete tMeshData;
    }
}


}
