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


TEST_CASE("Plane coincident to regular subdivision plane","[COINCIDENT]")
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
            moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

            // Setup XTK Model -----------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

            //Specify your decomposition methods and start cutting
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                   Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

            std::string tPrefix = std::getenv("XTKOUTPUT");
            std::string tMeshOutputFile = tPrefix + "/xtk_test_output_coincident.e";
            tCutMeshData->create_output_mesh(tMeshOutputFile);
            delete tCutMeshData;
            delete tMeshData;
    }
}

}
