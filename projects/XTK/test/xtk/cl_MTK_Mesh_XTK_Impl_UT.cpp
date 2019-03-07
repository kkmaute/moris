/*
 * cl_MTK_Mesh_XTK_Impl_UT.cpp
 *
 *  Created on: Feb 21, 2019
 *      Author: doble
 */

#include "catch.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_Sphere.hpp"

namespace xtk
{

TEST_CASE("Unit test of XTK's Interface with MTK","[XTK_MTK]")
{
if(par_size() < 2)
{

    real tRadius  = 0.55;
    real tXCenter = 2.0;
    real tYCenter = 2.0;
    real tZCenter = 2.0;
    Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh --------------------------------------------------------------------
    std::string tMeshFileName = "generated:2x2x2";
    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );

    // Setup XTK Model ----------------------------------------------------------------
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mVerbose = true;

    //Specify decomposition Method and Cut Mesh ---------------------------------------
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);
    tXTKModel.unzip_interface();
    tXTKModel.perform_basis_enrichment();

    moris::mtk::Mesh* tXTKToMTK = tXTKModel.get_xtk_as_mtk();

    for(moris::moris_index i = 0; i<(moris::moris_index)tXTKToMTK->get_num_entities(EntityRank::NODE); i++)
    {
//        std::cout<<"i = "<<i<<std::endl;
//        moris::print(tXTKToMTK->get_bspline_inds_of_node_loc_ind(i,EntityRank::INVALID),"basis functions");
//        moris::print(tXTKToMTK->get_t_matrix_of_node_loc_ind(i,EntityRank::INVALID),"basis functions");
    }

    moris::print(tXTKToMTK->get_entity_connected_to_entity_loc_inds(0,EntityRank::ELEMENT,EntityRank::NODE),"elem 0 nodes");

    delete tXTKToMTK;
    delete tMeshData;
}
                }
}

