/*
 * UT_XTK_Model_Mesh_Clusters.cpp
 *
 *  Created on: May 24, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_Sphere.hpp"

namespace xtk
{

TEST_CASE("Mesh Cluster Output","[XTK] [XTK_CLUSTER]")
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize<=4)
    {
            // Geometry Engine Setup ---------------------------------------------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius  = 0.25;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0.0;
            Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
            Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:1x1x4";
            moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
            tXTKModel.mVerbose = true;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            // output to exodus file ----------------------------------------------------------
            Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = true;
            tOutputOptions.mAddSideSets = false;
            tOutputOptions.mAddClusters = true;
            tOutputOptions.mAddParallelFields = true;

            moris::mtk::Integration_Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

            std::string tMeshOutputFile = "./xtk_exo/xtk_cluster_output.e";
            tCutMeshData->create_output_mesh(tMeshOutputFile);

            delete tCutMeshData;
            delete tMeshData;
        }
    }
}
