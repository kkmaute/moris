/*
 * fn_generate_shared_face_element_graph.hpp
 *
 *  Created on: May 16, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_
#define SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_

// XTKL: Linalg Includes
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"


// XTK includes
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Face_Registry.hpp"


namespace xtk
{


//FIXME: THIS FUNCTION NEEDS TO GO ONCE THE NEIGHBORHOOD IS GENERATED IN XTK
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
inline
moris::Matrix< moris::IndexMat >
generate_shared_face_element_pairs(
        moris::moris_index const & aFaceIndex,
        moris::moris_index const & aChildMeshIndex0,
        moris::moris_index const & aChildMeshIndex1,
        Cut_Mesh    const & aCutMesh,
        Matrix<IndexMat>  & aElementPairOrdinals)
{
    // Get references to the child meshes and needed connectivities
    Child_Mesh const & tChildMesh0 = aCutMesh.get_child_mesh(aChildMeshIndex0);
    moris::Matrix< moris::IndexMat > const & tFaceToNode0    = tChildMesh0.get_facet_to_node();
    moris::Matrix< moris::IndexMat > const & tElementToFace0 = tChildMesh0.get_element_to_facet();

    Child_Mesh const & tChildMesh1 = aCutMesh.get_child_mesh(aChildMeshIndex1);
    moris::Matrix< moris::IndexMat > const & tFaceToNode1    = tChildMesh1.get_facet_to_node();
    moris::Matrix< moris::IndexMat > const & tElementToFace1 = tChildMesh1.get_element_to_facet();

    // estimate maximum number of elements on face
     const uint tMaxElemOnFace = 100;

    // Allocate Matrixes
    moris::Matrix< moris::IndexMat > tFaceOrdinals(1,tMaxElemOnFace);
    moris::Matrix< moris::IndexMat > tChildrenElementCMInds(1,tMaxElemOnFace);
    moris::Matrix< moris::IdMat >    tChildrenElementIds(1,tMaxElemOnFace);

    // define variable for actual number of child elements on face
    uint tNumberOfChildElemsOnFace;

    // Get children elements attached to aFaceIndex on the side of child mesh index 0
    tChildMesh0.get_child_elements_connected_to_parent_facet(
            aFaceIndex,
            tNumberOfChildElemsOnFace,
            tChildrenElementIds,
            tChildrenElementCMInds,
            tFaceOrdinals);

    // Allocate output where top row is a child element indices in first child mesh and second row is for second child mesh index
    moris::Matrix< moris::IndexMat > tChildElementPairs(2,tNumberOfChildElemsOnFace);
    aElementPairOrdinals.resize(2,tNumberOfChildElemsOnFace);

    // Allocate downward map
    std::unordered_map<moris::moris_index,moris::moris_index> tChildElement0Map;
    for(moris::size_t i = 0; i<tNumberOfChildElemsOnFace; i++)
    {
        tChildElement0Map[tChildrenElementCMInds(0,i)] = i;
    }

    // Get the face to nodes of face on the interface
    moris::Matrix< moris::IndexMat > tFaceNodes(tNumberOfChildElemsOnFace,tFaceToNode0.n_cols());

    for(moris::size_t i = 0; i<tNumberOfChildElemsOnFace; i++)
    {
        moris::moris_index tFaceIndex = tElementToFace0(tChildrenElementCMInds(0,i),tFaceOrdinals(0,i));
        replace_row(tFaceIndex,tFaceToNode0,i,tFaceNodes);
    }

    //Allocate face registry and set first parent element's children on shared face information
    moris::Matrix< moris::IndexMat > tChildrenElementCMIndsTrans = moris::trans(tChildrenElementCMInds);

    Face_Registry tFaceRegistry(tNumberOfChildElemsOnFace,
                                tFaceNodes,
                                tChildrenElementCMIndsTrans);

    moris::Matrix< moris::IndexMat > & tFaceToElement = tFaceRegistry.get_face_to_element();


    for(moris::size_t j = 0; j<tNumberOfChildElemsOnFace; j++)
    {
        tChildElementPairs(0,j)   = tChildrenElementCMInds(j);
        aElementPairOrdinals(0,j) = tFaceOrdinals(j);
    }

    // Get children elements attached to aFaceIndex on the side of child mesh index 0
    tChildMesh1.get_child_elements_connected_to_parent_facet(
            aFaceIndex,
            tNumberOfChildElemsOnFace,
            tChildrenElementIds,
            tChildrenElementCMInds,
            tFaceOrdinals);

    for(moris::size_t i = 0; i<tNumberOfChildElemsOnFace; i++)
    {
        moris::size_t tFaceIndex = tElementToFace1(tChildrenElementCMInds(0,i),tFaceOrdinals(0,i));
        replace_row(tFaceIndex,tFaceToNode1,i,tFaceNodes);
    }

    moris::Matrix< moris::IndexMat > tFaceInds = tFaceRegistry.get_face_indices(tFaceNodes,true);

    // Establish pairs
    for(moris::size_t i = 0; i<tChildElementPairs.n_cols(); i++)
    {
        moris::moris_index tElementFrom0OnFace = (tFaceToElement)(i,0);
        moris::size_t tIndexInPairs = tChildElement0Map[tElementFrom0OnFace];

        (tChildElementPairs)(1,tIndexInPairs) = tChildrenElementCMInds(i);
        aElementPairOrdinals(1,tIndexInPairs) = tFaceOrdinals(i);
    }

    return tChildElementPairs;

                                   }


}
#endif /* SRC_XTK_FN_GENERATE_SHARED_FACE_ELEMENT_GRAPH_HPP_ */
