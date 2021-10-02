/*
 * cl_XTK_Model.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "paths.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Diagnostics.hpp"


#include "cl_GEN_Sphere.hpp"
#include "fn_GEN_Triangle_Geometry.hpp" // For surface normals

// Linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"
#include "fn_norm.hpp"
#include "op_minus.hpp"

#include "cl_Profiler.hpp"
#include "Child_Mesh_Verification_Utilities.hpp"


#include "fn_compute_interface_surface_area.hpp"

namespace xtk
{

    TEST_CASE("Regular Subdivision Method","[XTK] [REGULAR_SUBDIVISION_MODEL]")
    {
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
        MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

        if(tProcSize==1)
        {
            // Geometry Engine Setup -----------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius = 0.7;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0;
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
            tGeometry(0) = std::make_shared<moris::ge::Sphere>(tXCenter, tYCenter, tZCenter, tRadius);

            // Create Mesh ---------------------------------
            std::string tMeshFileName = "generated:1x1x1";
            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

            // Setup XTK Model -----------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);

            //Specify your decomposition methods and start cutting
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
            tXTKModel.decompose(tDecompositionMethods);

            Cut_Integration_Mesh* tCutMesh  = tXTKModel.get_cut_integration_mesh();

            // Do some testing
            size_t tNumNodesAfterDecompositionXTK = tCutMesh->get_num_entities(EntityRank::NODE,0);
            size_t tNumElementsAfterDecompositionXTK = tCutMesh->get_num_entities(EntityRank::ELEMENT,0);

            // 1 element was subdivided
            CHECK(tNumNodesAfterDecompositionXTK == 15);
            CHECK(tNumElementsAfterDecompositionXTK == 25); // 24 from template plus the original background cell 

            moris::Cell<std::shared_ptr<Matrix<DDRMat>>>* tCoords = tCutMesh->get_all_vertex_coordinates_loc_inds();

            moris::Matrix< moris::DDRMat > tExpectedNodeCoordinates(
                {{0, 0, 0},
                {1, 0, 0},
                {0, 1, 0},
                {1, 1, 0},
                {0, 0, 1},
                {1, 0, 1},
                {0, 1, 1},
                {1, 1, 1},
                {0.5, 0, 0.5},
                {1, 0.5, 0.5},
                {0.5, 1, 0.5},
                {0, 0.5, 0.5},
                {0.5, 0.5, 0},
                {0.5, 0.5, 1},
                {0.5, 0.5, 0.5}});

            CHECK(tCoords->size() == 15);
            for( moris::uint iNode = 0; iNode < tCoords->size(); iNode++ ) 
            {
                Matrix<DDRMat> tNodeCoords = *((*tCoords)(iNode));
                Matrix<DDRMat> tGoldNodeCoordinates = tExpectedNodeCoordinates.get_row(iNode);
                CHECK(equal_to(tNodeCoords,tGoldNodeCoordinates));
            }

            CHECK(interpolated_coordinate_check(tCutMesh));
         
         
         }
    }

    TEST_CASE("Regular Subdivision and Nodal Hierarchy Subdivision","[XTK] [CONFORMAL_MODEL]")
    {
        int tProcSize = moris::par_size();

        if(tProcSize<=2)
        {
            // Geometry Engine Setup ---------------------------------------------------------
            // Using a Levelset Sphere as the Geometry

            real tRadius  = 0.25;
            real tXCenter = 1.0;
            real tYCenter = 1.0;
            real tZCenter = 0.0;
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
            tGeometry(0) = std::make_shared<moris::ge::Sphere>(tXCenter, tYCenter, tZCenter, tRadius);

            // Create Mesh --------------------------------------------------------------------
            std::string tMeshFileName = "generated:1x1x4";
            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, NULL );
            std::string tBMOutputFile ="./xtk_exo/xtk_test_output_conformal_bm.e";
            tMeshData->create_output_mesh(tBMOutputFile);

            moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

            // Setup XTK Model ----------------------------------------------------------------
            size_t tModelDimension = 3;
            Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
            tXTKModel.mVerbose  =  false;

            //Specify decomposition Method and Cut Mesh ---------------------------------------
            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            // Access the Cut Mesh-------------------------------------------------------------
            Cut_Integration_Mesh* tCutMesh  = tXTKModel.get_cut_integration_mesh();

            // Do some testing
            size_t tNumNodesAfterDecompositionXTK = tCutMesh->get_num_entities(EntityRank::NODE,0);
            size_t tNumElementsAfterDecompositionXTK = tCutMesh->get_num_entities(EntityRank::ELEMENT,0);

            if(par_size() == 1)
            {   
                CHECK(tNumNodesAfterDecompositionXTK == 34 );
                CHECK(tNumElementsAfterDecompositionXTK == 46);
            }

            CHECK(interpolated_coordinate_check(tCutMesh));

            moris::Matrix<moris::DDRMat> tIsoContourThreshold = {{0.0}};
            moris::Matrix<moris::DDRMat> tIsoContourTolerance = {{1e-12}};
            CHECK(verify_interface_vertices(&tXTKModel,tIsoContourThreshold,tIsoContourTolerance));
        }
    }
}
