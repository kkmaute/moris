/*
 * fn_generate_cross_boundary_pairs.hpp
 *
 *  Created on: May 25, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_FN_GENERATE_CROSS_BOUNDARY_PAIRS_HPP_
#define SRC_XTK_FN_GENERATE_CROSS_BOUNDARY_PAIRS_HPP_



// XTKL: Linalg Includes
#include "linalg/cl_XTK_Matrix_Base.hpp"


// XTK includes
#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Mesh.hpp"

namespace xtk
{

/*
 * Given a shared face between two parent elements, constructs a boundary pair matrix which contains element indices. Also, constructs  If the parent element has children elements,
 * then the indices are in the child mesh local frame and if the parent element has no elements the index is processor local. The child elements are in their local
 * frame because this allows for efficient use of the internal flood-fill data which is also in this indexing system. Only creates a boundary pair if the bulk
 * phase is the same.
 * @param[in] aSharedFaceIndex - Shared Face Index
 * @param[in] aParentElement0  - Parent Element Index 0
 * @param[in] aParentElement1  - Parent Element Index 1
 * @param[in] aCutMesh         - Cut Mesh
 * @param[in] aXTKMesh         - XTK Mesh
 * @param[out] aBoundaryPairPhase - Bulk phase of the boundary pair
 * @param[out] Element pairs across shared face index (indexed according to the description above)
 *
 */
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
std::shared_ptr<Matrix_Base<Integer, Integer_Matrix>> generate_cross_boundary_pairs(Integer const & aSharedFaceIndex,
                                                                                    Integer const & aParentElement0,
                                                                                    Integer const & aParentElement1,
                                                                                    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix>    & aCutMesh,
                                                                                    XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix>    & aXTKMesh,
                                                                                    Matrix_Factory<Real,Integer,Real_Matrix,Integer_Matrix> & aMatrixFactory,
                                                                                    Matrix_Base<Integer,Integer_Matrix> & aBoundaryPairPhase)
{

    // Determine if both parent elements have children, or one has children or none
    bool tElement0HasChildren = aXTKMesh.entity_has_children(aParentElement0,EntityRank::ELEMENT);
    bool tElement1HasChildren = aXTKMesh.entity_has_children(aParentElement1,EntityRank::ELEMENT);

    // Parent element 1 variable initialization
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tElementInds0  = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceInds0     = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tElemPhaseInd0 = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceNodes0    = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceOrdinals0 = aMatrixFactory.create_integer_type_matrix_base(1,1);

    // Parent element 1 variable initialization
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tElementInds1  = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceInds1     = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tElemPhaseInd1 = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceNodes1    = aMatrixFactory.create_integer_type_matrix_base(1,1);
    std::shared_ptr<Matrix_Base<Integer,Integer_Matrix>> tFaceOrdinals1 = aMatrixFactory.create_integer_type_matrix_base(1,1);



    if(tElement0HasChildren && tElement1HasChildren)
    {
        Integer tChildMeshIndex0 = aXTKMesh.child_mesh_index(aParentElement0,EntityRank::ELEMENT);
        Integer tChildMeshIndex1 = aXTKMesh.child_mesh_index(aParentElement1,EntityRank::ELEMENT);

        std::cout<<"Parent Element 0 = "<< aParentElement0<<std::endl;
        std::cout<<"Parent Element 1 = "<< aParentElement1<<std::endl;
        std::cout<<"Shared Face = " <<aSharedFaceIndex<<std::endl;



        // Get children elements attached to aFaceIndex on the side of child mesh index 0
        aCutMesh.get_child_elements_connected_to_parent_face( tChildMeshIndex0,
                                                              aSharedFaceIndex,
                                                             *tElementInds0,
                                                             *tElemPhaseInd0,
                                                             *tFaceOrdinals0,
                                                             *tFaceInds0,
                                                             *tFaceNodes0);

        // Get children elements attached to aFaceIndex on the side of child mesh index 1
        aCutMesh.get_child_elements_connected_to_parent_face( tChildMeshIndex1,
                                                              aSharedFaceIndex,
                                                             *tElementInds1,
                                                             *tElemPhaseInd1,
                                                             *tFaceOrdinals1,
                                                             *tFaceInds1,
                                                             *tFaceNodes1);

        // Initialize Face Registry (used to generate face to element connectivity and assign unique face indices for the tied connectivity
        Integer_Matrix tElementIndsBase0 = transpose(tElementInds0->matrix_data());
        tElementInds0 = aMatrixFactory.create_integer_type_matrix_base(tElementIndsBase0);
        Face_Registry<Real,Integer, Real_Matrix, Integer_Matrix> tFaceRegistry(tFaceInds0->n_cols(),
                                                                               *tFaceNodes0,
                                                                               *tElementInds0,
                                                                               aMatrixFactory);

        // Get the face indices from child element 1 side of things
        std::shared_ptr<Matrix_Base<Integer, Integer_Matrix>> tFaceInds = tFaceRegistry.get_face_indices(*tFaceNodes1,true);
        print(tFaceInds,"tFaceInds");



    }

    else
    {
        std::cout<<"Mixed children parent element case not implemented"<<std::endl;
    }
    return aMatrixFactory.create_integer_type_matrix_base(1,1);
}

}

#endif /* SRC_XTK_FN_GENERATE_CROSS_BOUNDARY_PAIRS_HPP_ */
