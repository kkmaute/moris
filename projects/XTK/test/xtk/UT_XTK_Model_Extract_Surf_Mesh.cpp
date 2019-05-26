/*
 * UT_XTK_Model_Extract_Surf_Mesh.cpp
 *
 *  Created on: May 17, 2019
 *      Author: doble
 */

#include "catch.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_Plane.hpp"
#include "cl_Gyroid.hpp"

namespace xtk
{
TEST_CASE("Extract a surface mesh from XTK","[Extract_Surf]")
        {
//            real tXCenter = 4.2;
//            real tYCenter = 4.2;
//            real tZCenter = 4.2;
//
//            Plane tLevelsetSphere(tXCenter,tYCenter,tZCenter,1.0,1.0,1.0);
            Gyroid tGyroid;
            Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
            Geometry_Engine tGeometryEngine(tGyroid,tPhaseTable);

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:20x20x20|sideset:xXyYzZ";
            moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName );

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
            tXTKModel.mVerbose = true;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            moris::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                          Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            // select the bounding box sets
            moris::Cell<std::string> tBoundingSets = {"surface_1","surface_2","surface_3","surface_4","surface_5","surface_6"};

            // extract the surface to an obj file
            std::string tOutputObj = "./extract_surf.obj";
            tXTKModel.extract_surface_mesh_to_obj(tOutputObj,0,tBoundingSets);

            // dump full mesh
            moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
            std::string tPrefix = std::getenv("MORISOUTPUT");
            std::string tMeshOutputFile = tPrefix + "/surface_test_full.e";
            tCutMeshData->create_output_mesh(tMeshOutputFile);

            delete tCutMeshData;
        }
}



