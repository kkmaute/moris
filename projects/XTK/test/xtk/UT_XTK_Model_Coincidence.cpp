/*
 * cl_XTK_Model_Coincidence.cpp
 *
 *  Created on: Oct 26, 2018
 *      Author: doble
 */

#include "catch.hpp"


// XTKL: Mesh Includes
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Enums.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes
#include "cl_Matrix.hpp"

#include "geometry/cl_Plane.hpp"
#include "geometry/cl_Multi_Cylinder.hpp"

#include "cl_XTK_Model.hpp"

using namespace moris;
namespace xtk
{


TEST_CASE("Cylinders coincident with each other ","[COINCIDENT]")
{


    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize==1)
    {
            // Geometry Engine Setup -----------------------
            // Using a Levelset Sphere as the Geometry
            Cell<Cell<real>> tCenter = {{2,2,5}};
            Cell<real> tR =  {2.0};
            Cell<real> tL =  {9.0};
            Cell<Cell<real>> tAxis = {{0.0,0.0,1.0}};
            Multi_Cylinder tCylinder1(tCenter,{tR},{tL},tAxis);

            tCenter = {{2,2,5}};
            tR =  {2.0};
            tL =  {15.0};
            tAxis = {{0.0,0.0,1.0}};
            Multi_Cylinder tCylinder2(tCenter,{tR},{tL},tAxis);

            moris::Cell<Geometry*> tGeometryVector = {&tCylinder1, &tCylinder2};

            Phase_Table tPhaseTable (2,  Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine tGeometryEngine(tGeometryVector,tPhaseTable);

            // Create Mesh ---------------------------------
            std::string tMeshFileName = "generated:5x5x10";
            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );

            // Setup XTK Model -----------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
            tXTKModel.mVerbose = true;

            //Specify your decomposition methods and start cutting
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                   Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

            std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_coincident.e";
            std::cout<<"Output mesh: "<<tMeshOutputFile<<std::endl;
            tCutMeshData->create_output_mesh(tMeshOutputFile);
            delete tCutMeshData;
            delete tMeshData;
    }
}

TEST_CASE("Plane coincident with background mesh ","[COINCIDENT_PLANE]")
{


    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize==1)
    {
        real tXc = 0.0;
        real tYc = 0.0;
        real tZc = 0.0;
        real tXn = 0.0;
        real tYn = 1.0;
        real tZn = 0.0;

        Plane tPlane(tXc,tYc,tZc,tXn,tYn,tZn);

        Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
        Geometry_Engine tGeometryEngine(tPlane,tPhaseTable);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:1x1x10";
        moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
        tXTKModel.mVerbose = true;

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                               Subdivision_Method::C_HIERARCHY_TET4};
        tXTKModel.decompose(tDecompositionMethods);

        moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

        std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_coincident_plane.e";
        std::cout<<"Output mesh: "<<tMeshOutputFile<<std::endl;
        tCutMeshData->create_output_mesh(tMeshOutputFile);
        delete tCutMeshData;
        delete tMeshData;
    }
}

}
