/*
 * UT_XTK_Ghost_Stabilization.cpp
 *
 *  Created on: Mar 26, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"

#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_MTK_Writer_Exodus.hpp"


#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Plane.hpp"
//#include "cl_MGE_Geometry_Engine.hpp"
//
//#include "cl_Plane.hpp"
#include "cl_Mesh_Factory.hpp"

namespace xtk
{

TEST_CASE("Face oriented ghost stabilization","[GHOST]")
{
    if(par_size() == 1)
    {
    moris::Matrix<moris::DDRMat> tCenters = {{ 2.0,2.0,2.1 }};
    moris::Matrix<moris::DDRMat> tNormals = {{ 1.0,1.0,1.0 }};
    moris::ge::Plane<3> tPlane(tCenters,tNormals);

    moris::ge::GEN_Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    moris::ge::GEN_Geometry_Engine tGeometryEngine(tPlane,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:1x1x8|sideset:z";
    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName );


    // create model
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mVerbose  =  false;

    // decompose
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    tXTKModel.perform_basis_enrichment(EntityRank::NODE,0);

    tXTKModel.construct_face_oriented_ghost_penalization_cells();

    tXTKModel.get_ghost_stabilization().visualize_ghost_on_mesh(0);

    // Write mesh
    Writer_Exodus writer(&tXTKModel.get_enriched_integ_mesh());
    writer.write_mesh("", "./xtk_exo/xtk_test_ghost.exo");

    // Write the fields
    writer.set_time(0.0);
    writer.close_file();


    delete tMeshData;
    }


}
}
