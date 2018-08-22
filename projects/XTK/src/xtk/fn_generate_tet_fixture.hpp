/*
 * fn_generate_tet_fixture.hpp
 *
 *  Created on: Feb 1, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_FN_GENERATE_TET_FIXTURE_HPP_
#define SRC_XTK_FN_GENERATE_TET_FIXTURE_HPP_

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"


// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg/cl_XTK_MatrixBase_Utilities.hpp"

#include "xtk/cl_XTK_Model.hpp"

#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"


namespace xtk
{

//std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>>
//generate_tet_fixture(std::string aGeneratedString)
//{
//
//    Matrix_Factory<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMatrixFactory;
//    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine;
//
//    // Create Mesh ---------------------------------
//
//    Cell<std::string> tScalarFields(0);
//    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
//    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string(tMatrixFactory,
//                                                                                                                                                aGeneratedString,
//
//    // Setup XTK Model -----------------------------
//    size_t tModelDimension = 3;
//    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,*tMeshData,tGeometryEngine,tMatrixFactory);
//
//    //Specify your decomposition methods and start cutting
//    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::T_GENERATE_TET_MESH_FROM_HEX};
//    tXTKModel.decompose(tDecompositionMethods,false);
//
//    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> & tCutMesh = tXTKModel.get_cut_mesh();
//
//    size_t tNumNodePerElement = 4;
//    size_t tNumElements = tCutMesh.get_num_entities(EntityRank::ELEMENT);
//    std::shared_ptr<xtk::Matrix_Base<xtk::size_t, Default_Matrix_Integer>> tElementConnectivity = tMatrixFactory.create_integer_type_matrix( tNumNodePerElement,tNumElements );
//
//    size_t tNumChildMeshes = tCutMesh.get_num_simple_meshes();
//    size_t tCount = 0;
//    for(size_t i = 0; i<tNumChildMeshes; i++)
//    {
//        std::shared_ptr<Matrix_Base<size_t, Default_Matrix_Integer>> tChildElemNode =  tCutMesh.get_full_connectivity(i,EntityRank::ELEMENT,EntityRank::NODE,2);
//
//        for(size_t j = 0; j<tChildElemNode->get_num_rows(); j++)
//        {
//            replace_row(j,*tChildElemNode,tCount,*tElementConnectivity);
//            tCount++;
//        }
//    }
//
//    xtk::print(tElementConnectivity,"Output");
//
//    // Nodes connected to Elements
//
//    // Element ownership ( both 0 because there is only one processor)
//    std::shared_ptr<xtk::Matrix_Base<xtk::size_t, Default_Matrix_Integer>> tElemProcs = tMatrixFactory.create_integer_type_matrix( 1,tNumElements,0 );
//
//    // Get node index to id map
//    xtk::Matrix_Base<xtk::size_t, Default_Matrix_Integer> const & tNodeLocaltoGlobal = tMeshData->get_local_to_global_map(EntityRank::NODE);
//
//    // Get node coordinates
//    std::shared_ptr<xtk::Matrix_Base<xtk::real, Default_Matrix_Real>> tCoords = tMeshData->get_all_node_coordinates_loc_inds();
//
//    // Construct mesh from data
//    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tTetMeshData = tMeshBuilder.build_mesh_from_data(tMatrixFactory, tModelDimension, *tElementConnectivity, *tCoords, tNodeLocaltoGlobal,*tElemProcs, {"block_1"});
//
//    return tTetMeshData;
//}

}



#endif /* SRC_XTK_FN_GENERATE_TET_FIXTURE_HPP_ */
