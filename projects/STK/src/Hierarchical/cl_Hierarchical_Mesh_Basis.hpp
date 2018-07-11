/*
 * cl_Hierarchical_Mesh_Basis.hpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_BASIS_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_BASIS_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_BoostBitset.hpp" // CON/src

namespace moris
{

    class Hierarchical_Mesh_Basis
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Basis()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Basis() = default;

        /**
         * Provides the number of basis functions within all levels from level 0 until aLevel
         *
         * @param[in] aLevel              Level of the basis functions.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] basis_number        Number of basis functions within all Levels until aLevel.
         *
         */
        static uint
        give_number_of_basis(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection);

        /**
         * Provides the neighbours of an basis
         *
         * @param[in] aBasisId              Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial         Polynomial degree
         * @param[in] aBuffer             Provides the number of layers around the basis
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Basis neighbour        Basis Id plus neighbor Ids (9 for aModelDim = 2, 27 for aModelDim = 3).\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElement, 3D: Element(13) = aElement;
         *
         */
        Mat<uint>
        give_neighbor_of_basis(uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                uint const & aBuffer,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the level of the basis function
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] level        Level of the element from the hierarchical mesh.
         *
         */
        uint
        give_basis_level(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the position of a basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Basis_position        Position I,J,K.
         *
         */
        Mat<uint>
        give_position_of_basis(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the basis function Id with the position i,j,k
         *
         * @param[in] aLevel             Level of the basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] tBasis       Basis function Id.
         *
         */
        uint
        give_basis_of_position(
                uint const & aLevel,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aIJKPosition) const;

        /**
         * Provides the basis function Id of the parent (Works only for a linear polynomial degree) (Returns UINT_MAX if there is no parent)
         *
         * @param[in] aBasisId            Basis function number.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] tBasis            Basis function Id of the parent.
         *
         */
        uint
        give_basis_of_parent(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the elements, which have support with the basis function from a tensorial grid
         *
         * @param[in] aBasisId            Basis function Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] tElements        Elements with support of basis function Id.
         *
         */
        Mat<uint>
        give_element_of_basis(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the coordinates for a specific basis function Id
         *
         * @param[in] aBasisId                Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection          NumElements=[Number of elements in x,y,z direction].
         * @param[in] aModelDimensions           Dimensions of the whole domain
         * @param[in] aModelDimensions_Offset    Offset of the problem, which is embedded in the outer layer of elements (it is needed to determine the point where the first basis functions sits, Point of origin (default = 0,0,0))
         *
         * @param[out] tCoordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_coordinate_from_basis(uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<real> const & aModelDimensions,
                Mat<real> const & aModelDimensions_Offset) const;

        /**
         * Defines owner for all basis function Ids from the vetor "aNodalLocaltoGlobal"
         *
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aNodalLocaltoGlobal   A list of node Ids, which are sitting on a proc (Needed for MTK/STK)
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         * @param[in] aBasisActive          A bitset with active basis functions
         *
         * @param[out] tNodeProcs           A vector with ownership of the node Ids
         *
         */
        Mat<uint>
        give_basis_proc_owner(
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aNodalLocaltoGlobal,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp,
                uint const & aNumberBasis,
                BoostBitset const & aBasisActive) const;

        /**
         * Defines owner of a basis function Id
         *
         * @param[in] aBasisId                  Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         *
         * @param[out] tProcOwner           Owner of the basis function id
         *
         */
        uint
        give_basis_proc_owner(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp) const;

        /**
         *  Defines a list of processors, which share a basis function Id
         *
         * @param[in] aBasisId                  Basis function Id
         * @param[in] aModelDim                  Number of dimensions.
         * @param[in] aPolynomial           Polynomial degree of the basis functions.
         * @param[in] aProcNeighbours       A list of proc neighbors
         * @param[in] aDecomp               Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aNumberBasis          Number of basis functions for the current level
         *
         * @param[out] tProcShare           Who shares this basis function
         *
         */
        Mat<uint>
        give_basis_share(
                uint const & aBasisId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aProcNeighbour,
                Mat<uint> const & aDecomp ) const;

        // FIXME: add descriptions
        uint
        give_symmetry_master_of_basis(
                const uint      & aBasisId,
                const uint      & aModelDim,
                const uint      & aPolynomial,
                const uint      & aSymmetryPlane,
                const uint      & aSymmetryIndex,
                const Mat<uint> & aNumberOfElementsPerDirection,
                const Mat<uint> & aNumberOfBasisPerDirection ) const;

        // FIXME: add descriptions
        uint
        give_symmetry_slave_of_basis(
                const uint      & aBasisId,
                const uint      & aModelDim,
                const uint      & aPolynomial,
                const uint      & aSymmetryPlane,
                const uint      & aSymmetryIndex,
                const Mat<uint> & aNumberOfElementsPerDirection,
                const Mat<uint> & aNumberOfBasisPerDirection ) const;

        // FIXME: add descriptions
        bool
        basis_is_symmetry_master(
                const uint      & aBasisId,
                const uint      & aModelDim,
                const uint      & aPolynomial,
                const uint      & aSymmetryPlane,
                const uint      & aSymmetryIndex,
                const Mat<uint> & aNumberOfElementsPerDirection ) const;
    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_BASIS_HPP_ */
