/*
 * cl_Hierarchical_Mesh_Element.hpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_ELEMENT_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_ELEMENT_HPP_

#include "algorithms.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_unique.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_BoostBitset.hpp" // CON/src
#include "cl_Hierarchical_Mesh_Basis.hpp"

#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Hierarchical_Mesh_Element
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Element()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Element() = default;

        /**
         * Provides the basis function IDs of an element from a tensorial grid
         *
         * @param[in] aElementId            Element Id.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aPolynomial          Polynomial degree of the basis functions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         *
         * @param[out] basis        Basis function IDs of an element.
         *
         */
        Mat<uint>
        give_basis_of_element(
                uint const & aElementId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;

        /**
         * Provides the coordinate in the middle of an element, it is needed to have a rough idea where the element sits, needed for "Delete Domain in Hierarchical_Mesh_Main::create_mesh()"
         *
         * @param[in] aElementId        A specific element.
         * @param[in] aModelDim            Dimension of the mesh
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aDomainRange        Domain range of the problem (length, height, thickness in x,y,z direction)
         * @param[in] aDomainOffset     Offset of the problem, which is embedded in the outer layer of elements (it is needed to determine the point where the first basis functions sits, Point of origin (default = 0,0,0))
         *
         * @param[out] coordinates       Coordinates in x,y,z direction.
         *
         */
        Mat<real>
        give_middlecoordinate_from_element(
                uint const & aElementId,
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<real> const & aDomainRange,
                Mat<real> const & aDomainOffset) const;

        /**
         * Provides the element on level zero, which has support for the specific coordinate (independent if it is active or deactive, needed for filter)
         *
         * @param[in] aModelDim               Number of Dimensions
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
         * @param[in] aDomainRange        Domain range of the problem (length, height, thickness in x,y,z direction)
         * @param[in] aDomainOffset     Offset of the problem, which is embedded in the outer layer of elements (it is needed to determine the point where the first basis functions sits, Point of origin (default = 0,0,0))
         * @param[in] aPointOfOrigin     Point of origin of the active elements (not of the layer)
         * @param[in] aCoordinate       A user specific coordinate
         *
         * @param[out] aElementId         Element Id of the active element
         *
         */
        uint
        give_element_for_coordinate_on_level_zero(
                uint const & aModelDim,
                Mat<uint> const & aNumberOfElementsPerDirection,
                Mat<real> const & aDomainRange,
                Mat<real> const & aDomainOffset,
                Mat<real> const & aPointOfOrigin,
                Mat<real> const & aCoordinate) const;

        /**
           * Provides all active neighbors (includes all levels)
           *
           * @param[in] aElementId            Element number.
           * @param[in] aModelDim                Number of dimensions.
           * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
           * @param[in] aLevel              Highest level for estimating the possible number of neighbors
           * @param[in] aElementIdActive     Bitset with active elements
           *
           * @param[out] Element_neighbor   Active neighbors
           *
           */
          Mat<uint>
          give_active_face_neighbor_of_element(
                  uint const & aElementId,
                  uint const & aModelDim,
                  Mat<uint> const & aNumberOfElementsPerDirection,
                  uint const & aLevel,
                  BoostBitset const & aElementActive) const;

          /**
             * Provides all active neighbors (includes all levels)
             *
             * @param[in] aElementId            Element number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
             * @param[in] aLevel              Highest level for estimating the possible number of neighbors
             * @param[in] aElementIdActive     Bitset with active elements
             *
             * @param[out] Element_neighbor   Active neighbors
             *
             */
            Mat<uint>
            give_active_neighbor_of_element(
                    uint const & aElementId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection,
                    uint const & aLevel,
                    BoostBitset const & aElementActive) const;

          /**
            * Provides the active element, which has support for the specific coordinate
            *
            * @param[in] aModelDim               Number of Dimensions
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
            * @param[in] aDomainRange        Domain range of the problem (length, height, thickness in x,y,z direction)
            * @param[in] aDomainOffset     Offset of the problem, which is embedded in the outer layer of elements (it is needed to determine the point where the first basis functions sits, Point of origin (default = 0,0,0))
            * @param[in] aPointOfOrigin     Point of origin of the active elements (not of the layer)
            * @param[in] aCoordinate       A user specific coordinate
            * @param[in] aElementIdActive    Bitset with active elements
            *
            * @param[out] aElementId         Element Id of the active element
            *
            */
           uint
           give_active_element_for_coordinate(
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   Mat<real> const & aDomainRange,
                   Mat<real> const & aDomainOffset,
                   Mat<real> const & aPointOfOrigin,
                   Mat<real> const & aCoordinate,
                   BoostBitset const & aElementActive) const;

           /**
            * Provides a stencil of neighbors of an element (Idea comes from SIMP procedure and feature size resolution)
            *
            * @param[in] aFeatureResolution    Size of a feature resolution
            * @param[in] aElementId            Element Id.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aPolynomial          Polynomial degree of the basis functions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
            * @param[in] ElementActive       Bitset with active elements
            *
            * @param[out] Element_stencil     Gives a matrix with neighbour elements. In each row is one specific direction ( 2D: 2 axes + 2 diagonals, 3D: 3 axes + 9 diagonals)
            *
            */
           Mat<uint>
           give_neighbor_stencil_of_element(
                   uint const & aFeatureResolution,
                   uint const & aElementId,
                   uint const & aModelDim,
                   uint const & aPolynomial,
                   Mat<uint> const & aNumberOfElementsPerDirection,
                   BoostBitset & aElementActive);
    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_ELEMENT_HPP_ */
