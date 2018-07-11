/*
 * cl_Base_Mesh_Element.hpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_BASE_MESH_ELEMENT_HPP_
#define SRC_MESH_CL_BASE_MESH_ELEMENT_HPP_

#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_unique.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_BoostBitset.hpp" // CON/src

namespace moris
{

    class Base_Mesh_Element
    {
    protected:

    public:
        /**
         * Hierarchical_Mesh constructor
         */
        Base_Mesh_Element()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Base_Mesh_Element() = default;

        /**
         * Provides the element with the position i,j,k
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomialDegree   Polynomial degree of the b-spline.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         * @param[in] aIJKPosition        Position I,J,K.
         *
         * @param[out] Bspline            Element id.
         *
         */
        uint
        give_element_of_position(
                uint const & aLevel,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> const & aIJKPosition) const;

        /**
         * Provides the position of an element from a tensorial grid
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Element_position        Position I,J,K.
         *
         */
        Mat<uint>
        give_position_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the level of the element from the hierarchical mesh.
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] level        Level of the elemen Id t from the hierarchical mesh.
         *
         */
        uint
        give_element_level(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the parent of an element
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Parent        Parent of an element Id.
         *
         */
        uint
        give_parent_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the parent of an element for a specific level
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aWhichLevel         User defined level on which the parent needs to be determined
         *
         * @param[out] Parent        Parent of an element.
         *
         */
        uint
        give_parent_of_level_x(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                uint const & aWhichLevel);

        /**
         * Provides the parent of an element and the relation to the child
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Parent          A vector with two entries. First entry is the parent of the element Id and second entry is the child of the parent relation (Which child is element Id of parent)
         *
         */
        Mat<uint>
        give_parent_child_realation_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection);

        /**
         * Provides the number of elements within all levels from level 0 until "aLevel"
         *
         * @param[in] aLevel              Level of the element.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] element_number        Number of elements within all Levels until aLevel.
         *
         */
        static uint
        give_number_of_elements(
                uint const & aLevel,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection);

        /**
         * Provides the children of an element
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Children        Childrens of an element (4 children for Dim=2 and 8 children for Dim=3).
         *
         */
        Mat<uint>
        give_children_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the neighbors of an element (including diagonal neighbors)
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aBuffer             Provides the number of layers around the element
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] Element_neighbour        Element plus neighbors (9 for aModelDim = 2, 27 for aModelDim = 3) if a buffer of one is used.\n
         *                                2D:...----------------------------------...3D:....-----plain 1----------------------.....-----plain 2-------------------------\n.....-----plain 3-------------------------\n
         *                                ......|Element(6) Element(7) Element(8)|..........|Element(6) Element(7) Element(8)|.....|Element(15) Element(16) Element(17)|\n.....|Element(24) Element(25) Element(26)|\n
         *                                ......|Element(3) Element(4) Element(5)|..........|Element(3) Element(4) Element(5)|.....|Element(12) Element(13) Element(14)|\n.....|Element(21) Element(22) Element(23)|\n
         *                                ......|Element(0) Element(1) Element(2)|..........|Element(0) Element(1) Element(2)|.....|Element(9) Element(10) Element(11) |\n.....|Element(18) Element(19) Element(20)|\n
         *                                ......|--------------------------------|..........|--------------------------------|.....|-----------------------------------|  .....|-----------------------------------|\n
         *                                2D: Element(4) = aElementId, 3D: Element(13) = aElementId;
         *
         */
        Mat<uint>
        give_neighbor_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aBuffer,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the ownership of an element
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aProcRange          Vector with a element range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aProcNeighbours     A list of proc neighbors
         *
         * @param[out] Rank        Proc rank
         *
         */
        uint
        give_element_owner(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> & aProcRange,
                Mat<uint> const & aProcNeighbours);

        /**
         * Provides a list, which share an element (owner plus aura procs) - An proc can identify owner and share only if the element lives on the proc or in the aura of that proc
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
         * @param[in] aProcNeighbours     A list of proc neighbors
         *
         * @param[out] Share        Vector with procs, which share the element
         *
         */
        Mat<uint>
        give_element_share(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<uint> & aProcRange,
                Mat<uint> const & aProcNeighbors);

        /**
         * Provides the neighbours of an element, which are connected through the faces.
         * (Therefore, no diagonal neighbors)
         * This function uses the function "give_neighbour_of_element" and extracts only the necessary elements
         *
         * @param[in] aElementId            Element number.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Element_neighbour        Element plus neighbours (4 for aModelDim = 2, 6 for aModelDim = 3).
         *
         */
        Mat<uint>
        give_face_neighbor(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

    };
}

#endif /* SRC_MESH_CL_BASE_MESH_ELEMENT_HPP_ */
