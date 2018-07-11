/*
 * cl_Hierarchical_Mesh_TMatrix.hpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_TMATRIX_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_TMATRIX_HPP_
#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Cell.hpp" // CON/src
#include "cl_Hierarchical_Mesh_Element.hpp"
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Bspline.hpp" // MOD/src
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Hierarchical_Mesh_TMatrix
    {
    protected:

    public:
        //Create Object of Element
        Hierarchical_Mesh_Element mHMRElement;
        //Create Object of base Element
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_TMatrix()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_TMatrix() = default;

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... in a reordered fashion (classical FEM, counter clockwise numbering!)
         *
         * @param[in] aElementId                       Element Id.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Tmatrix_and_IdField_Reorder(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the truncated T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... in a reordered fashion (classical FEM, counter clockwise numbering!)
         *
         * @param[in] aElementId                       Element Id.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Truncated_Tmatrix_and_IdField_Reorder(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the Id field, to know the active basis functions, which have support in the element
         *
         * @param[in] aElementId                       Element Id.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         *
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_IdField(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<uint> & aIdField);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,....
         *
         * @param[in] aElementId                       Element number.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Tmatrix_and_IdField(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the truncated T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,....
         *
         * @param[in] aElementId                       Element Id.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Truncated_Tmatrix_and_IdField(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... and projects this T-Matrix of the design variables to the FEM mesh
         *
         * @param[in] aElementId                       Element number.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Tmatrix_and_IdField_DesignToFEMProjection(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the truncated T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... and projects this T-Matrix of the design variables to the FEM mesh
         *
         * @param[in] aElementId                       Element number.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void give_Truncated_Tmatrix_and_IdField_DesignToFEMProjection(uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... and projects this T-Matrix of the design variables to the FEM mesh
         *
         * @param[in] aElementId                       Element number.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         * @param[in] aNaturalCoordinate               A list of natural coordinates (Rows defines the number of natural coordinates, 0 <= xi <= 1, and the columns the direction x,y,z)
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Tmatrix_and_IdField_SpecificPoint(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aBasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Mat<real> const & aNaturalCoordinate,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates the T-matrix for a 1,2 or 3 dimensional problem for p = 1,2,.... and projects this T-Matrix of the design variables to the FEM mesh
         *
         * @param[in] aElementId                       Element number.
         * @param[in] aModelDim                        Dimension of model
         * @param[in] aPolynomial                      Polynomial degree
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aBasisActive                     Bitset with active basis functions
         * @param[in] aTMatrixParentChild              T-Matrices of a parent/child relations
         * @param[in] aNaturalCoordinate               A list of natural coordinates (Rows defines the number of natural coordinates and the columns the direction x,y,z)
         *
         * @param[out] aTMatrix        T-Matrix of the element
         * @param[out] aIdField        List of basis functions, which have support in this element
         *
         */
        void
        give_Truncated_Tmatrix_and_IdField_SpecificPoint(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & BasisActive,
                Mat<real> & aTMatrix,
                Mat<uint> & aIdField,
                Mat<real> const & aNaturalCoordinate,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Creates a list of elements (children) and their respective T-Matrix to calculate the node property on a coarser level (L2-Projection needs this)
         *
         * @param[in] aElementId      Element number.
         * @param[in] aModelDim              Dimension of model
         * @param[in] aPolynomial       Polynomial degree
         * @param[in] aLevelFem         Level of FEM analysis from last optimization step
         * @param[in] aLevelDesign      Level of Design variables from last optimization step
         * @param[in] aNumberOfElementsPerDirection      Number of elements in each direction
         * @param[in] aElementActiveLastStep   Bitset with active elements of last optimization step
         *
         * @param[out] aTMatrixOfChildren        T-Matrix of the children elements
         * @param[out] aListOfChildren        List of elements, which are an old active children of the current element
         *
         */
        void
        give_Tmatrix_and_IdField_Design_L2Projection_Coarsening(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                uint const & aLevelFem,
                uint const & aLevelDesign,
                Mat<uint> const & aNumberOfElementsPerDirection,
                BoostBitset const & aElementActiveLastStep,
                Cell<Mat<real>> & aTMatrixOfChildren,
                Mat<uint> & aListOfChildren,
                Cell<Mat<real>> const & aTMatrixParentChild);

        /**
         * Give the projection matrix for the design variables to the FEM mesh
         *
         * @param[in] aModelDim              Number of Dimensions
         * @param[in] aPolynomial       Polynomial degree
         *
         * @param[out] tProject     Projection matrix
         *
         */
        static Mat<real>
        give_projection_matrix(
                uint const & aModelDim,
                uint const & aPolynomial);

        /**
         * Give the projection matrix for the design variables to the FEM mesh
         *
         * @param[in] aModelDim              Number of Dimensions
         * @param[in] aPolynomialLagrange    Polynomial degree of Lagrange basis
         * @param[in] aPolynomial            Polynomial degree of Bspline basis
         *
         * @param[out] tProject     Projection matrix
         *
         */
        Mat<real>
        give_projection_matrix_new(
                uint const & aModelDim,
                uint const & aPolynomialLagrange,
                uint const & aPolynomial) const;

        /**
         * Gives the T-matrices of all children for 1D, 2D or 3D, for the polynomial degrees p = 1, 2, 3
         *
         * @param[in] aPolynomial   Polynomial degree of the b-splines for each direction p_1 = p_2 = p_3.
         * @param[in] aModelDim     Number of dimensions.
         *
         * @param[out] T_child            Gives the Tmatrices of all childs in a Cell. The order is:\n
         *                                2D:...------------------... 3D:...-----plain 1------.....-----plain 2------\n
         *                                ......|Child 3..Child 4|..........|Child 3..Child 4|.....|Child 7..Child 8|\n
         *                                ......|................|..........|................|.....|................|\n
         *                                ......|Child 1..Child 2|..........|Child 1..Child 2|.....|Child 5..Child 6|\n
         *                                ......|----------------|..........|----------------|.....|----------------|
         *
         */
        static Cell<Mat<real>>
        give_Tmatrices_of_childs(
                uint const & aPolynomial,
                uint const & aModelDim);

        /**
         * Provides a vector for reordering (needed for paraview)
         *
         */
        Mat<uint>
        give_vector_for_reorder(
                uint const & aModelDim,
                uint const & aPolynomial) const;
    };

}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_TMATRIX_HPP_ */
