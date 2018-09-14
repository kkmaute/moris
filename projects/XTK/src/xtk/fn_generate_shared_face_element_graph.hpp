/*
 * fn_generate_shared_face_element_graph.hpp
 *
 *  Created on: May 16, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_
#define SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_

// XTKL: Linalg Includes
#include "linalg/cl_XTK_Matrix.hpp"
#include "fn_trans.hpp"


// XTK includes
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Face_Registry.hpp"


namespace xtk
{

/*
 * Given a shared face from the background mesh and two children meshes with the shared face, construct cross parent element
 * boundary.
 * @param[in] aFaceIndex       - Face index shared across parent element boundary
 * @param[in] aChildMeshIndex0 - Child mesh index 0 connected to boundary
 * @param[in] aChildMeshIndex1 - Child mesh index 1 connected to boundary
 * @param[in] aCutMesh         - Mesh containing elements around the interface
 * @param[in] aBackgroundMesh  - Background mesh (Lagrangian Mesh)
 * @param[in] aMatrixFactory   - Means of creating matrix objects
 * @param[out] Cross shared face child element pairs (row 0 element indices from child mesh index 0,
 *                                                    row 1 element indices from child mesh index 1)
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
moris::Matrix< Integer_Matrix >
generate_shared_face_element_pairs(Integer const & aFaceIndex,
                                   Integer const & aChildMeshIndex0,
                                   Integer const & aChildMeshIndex1,
                                   Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix>    const & aCutMesh)
                                   {
    // Get references to the child meshes and needed connectivities
    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> const & tChildMesh0 = aCutMesh.get_child_mesh(aChildMeshIndex0);
    moris::Matrix< Integer_Matrix > const & tFaceToNode0    = tChildMesh0.get_face_to_node();
    moris::Matrix< Integer_Matrix > const & tElementToFace0 = tChildMesh0.get_element_to_face();

    Child_Mesh_Test<Real, Integer, Real_Matrix, Integer_Matrix> const & tChildMesh1 = aCutMesh.get_child_mesh(aChildMeshIndex1);
    moris::Matrix< Integer_Matrix > const & tFaceToNode1    = tChildMesh0.get_face_to_node();
    moris::Matrix< Integer_Matrix > const & tElementToFace1 = tChildMesh0.get_element_to_face();

    // Allocate Matrixes
    moris::Matrix< Integer_Matrix > tFaceOrdinals(1,1);
    moris::Matrix< Integer_Matrix > tChildrenElementCMInds(1,1);
    moris::Matrix< Integer_Matrix > tChildrenElementIds(1,1);
    moris::Matrix< Integer_Matrix > tChildrenElementPhaseIndex(1,1);

    // Get children elements attached to aFaceIndex on the side of child mesh index 0
    tChildMesh0.get_child_elements_connected_to_parent_face(aFaceIndex,
                                                            tChildrenElementIds,
                                                            tChildrenElementCMInds,
                                                            tFaceOrdinals);

    // Allocate output where top row is a child element indices in first child mesh and second row is for second child mesh index
    moris::Matrix< Integer_Matrix > tChildElementPairs(2,tChildrenElementCMInds.n_cols());

    // Allocate downward map
    std::unordered_map<Integer,Integer> tChildElement0Map;
    for(Integer i = 0; i<tChildrenElementCMInds.n_cols(); i++)
    {
        tChildElement0Map[tChildrenElementCMInds(0,i)] = i;
    }

    // Get the face to nodes of face on the interface
    moris::Matrix< Integer_Matrix > tFaceNodes(tChildrenElementCMInds.n_cols(),tFaceToNode0.n_cols());

    for(Integer i = 0; i<tChildrenElementCMInds.n_cols(); i++)
    {
        Integer tFaceIndex = tElementToFace0(tChildrenElementCMInds(0,i),tFaceOrdinals(0,i));
        replace_row(tFaceIndex,tFaceToNode0,i,tFaceNodes);
    }

    //Allocate face registry and set first parent element's children on shared face information
    Integer tNumFaces = tChildrenElementIds.n_cols();
    Integer_Matrix tChildrenElementCMIndsTrans = trans(tChildrenElementCMInds);
    tChildrenElementCMInds.matrix_data() = tChildrenElementCMIndsTrans;
    Face_Registry<Real, Integer, Real_Matrix, Integer_Matrix> tFaceRegistry(tNumFaces,
                                                                           tFaceNodes,
                                                                           tChildrenElementCMInds);

    moris::Matrix< Integer_Matrix > & tFaceToElement = tFaceRegistry.get_face_to_element();


    for(Integer j = 0; j<tChildrenElementCMInds.n_rows(); j++)
    {
        (tChildElementPairs)(0,j) = tChildrenElementCMInds(j,0);
    }

    // Get children elements attached to aFaceIndex on the side of child mesh index 0
    tChildMesh1.get_child_elements_connected_to_parent_face(aFaceIndex,
                                                            tChildrenElementIds,
                                                            tChildrenElementCMInds,
                                                            tFaceOrdinals);

    for(Integer i = 0; i<tChildrenElementCMInds.n_cols(); i++)
    {
        Integer tFaceIndex = tElementToFace1(tChildrenElementCMInds(0,i),tFaceOrdinals(0,i));
        replace_row(tFaceIndex,tFaceToNode1,i,tFaceNodes);
    }

    moris::Matrix< Integer_Matrix > tFaceInds = tFaceRegistry.get_face_indices(tFaceNodes,true);


    // Establish pairs
    for(Integer i = 0; i<tChildElementPairs.n_cols(); i++)
    {
        Integer tElementFrom0OnFace = (tFaceToElement)(i,0);
        Integer tIndexInPairs = tChildElement0Map[tElementFrom0OnFace];

        (tChildElementPairs)(1,tIndexInPairs) = tChildrenElementCMInds(0,i);
    }


    return tChildElementPairs;

                                   }


}
#endif /* SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_ */
