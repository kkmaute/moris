
#include "catch.hpp"

// Xtk includes
#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"

// Geometry
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"

// Linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

using namespace moris;
using namespace xtk;

TEST_CASE("XTK Model Unzip Interface","[unzip_xtk]")
{
     real tRadius  = 0.75;
     real tXCenter = 1.0;
     real tYCenter = 1.0;
     real tZCenter = 1.0;
     Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
     Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
     Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

     // Create Mesh --------------------------------------------------------------------
     std::string tMeshFileName = "generated:1x1x2";
     moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName );

     // Setup XTK Model ----------------------------------------------------------------
     size_t tModelDimension = 3;
     Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
     tXTKModel.mVerbose = true;

     //Specify decomposition Method and Cut Mesh ---------------------------------------
     xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
                                                                 Subdivision_Method::C_HIERARCHY_TET4};
     tXTKModel.decompose(tDecompositionMethods);

     // Test interface information prior to unzipping
     Background_Mesh const & tBM = tXTKModel.get_background_mesh();

     // get the interface nodes
     moris::Matrix< moris::IndexMat > tGoldInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35}});
     moris::Matrix< moris::IndexMat > tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);
     CHECK(all_true(tGoldInterfaceNodes == tInterfaceNodes));

     // Unzip interface
     tXTKModel.unzip_interface();

     // Now there should be twice as many nodes
     moris::Matrix< moris::IndexMat > tGoldUnzippedInterfaceNodes({{25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
                                                                    36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46}});
     tInterfaceNodes = tBM.get_interface_nodes_loc_inds(0);
     CHECK(all_true(tGoldUnzippedInterfaceNodes == tInterfaceNodes));


     // Output to file
     moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
     std::string tPrefix = std::getenv("MORISOUTPUT");
     std::string tMeshOutputFile = tPrefix + "/xtk_test_unzip.e";
     tCutMeshData->create_output_mesh(tMeshOutputFile);
     delete tCutMeshData;
     delete tMeshData;
}
