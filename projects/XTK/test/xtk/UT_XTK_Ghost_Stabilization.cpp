/*
 * UT_XTK_Ghost_Stabilization.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"

#include "../projects/GEN/src/geometry/cl_GEN_Geometry.hpp"
#include "../projects/GEN/src/geometry/cl_GEN_Plane.hpp"
//#include "cl_MGE_Geometry_Engine.hpp"
//
//#include "cl_Plane.hpp"
#include "cl_Mesh_Factory.hpp"

namespace xtk
{

TEST_CASE("Face oriented ghost stabilization","[GHOST]")
{
    moris::Matrix<moris::DDRMat> tCenters = {{ 2.0,2.0,2.0 }};
    moris::Matrix<moris::DDRMat> tNormals = {{ 1.0,1.0,1.0 }};
    moris::ge::Plane<3> tPlane(tCenters,tNormals);

    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tPlane,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:1x1x8|sideset:z";
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );

    std::string tPrefix = std::getenv("MORISOUTPUT");
    std::string tBackgroundFile = tPrefix + "/xtk_test_ghost_background.e";
    tMeshData->create_output_mesh(tBackgroundFile);

    // create model
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mVerbose  =  false;

    // decompose
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.construct_face_oriented_ghost_penalization_cells();

    // output to exodus file ----------------------------------------------------------
    Output_Options tOutputOptions;
    tOutputOptions.mAddNodeSets = true;
    tOutputOptions.mAddSideSets = true;

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh(tOutputOptions);

    std::string tMeshOutputFile = tPrefix + "/xtk_test_ghost.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);

    delete tMeshData;
    delete tCutMeshData;


}
}
