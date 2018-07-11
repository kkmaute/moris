/*
 * cl_Base_Mesh_Edge.hpp
 *
 *  Created on: Feb 23, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_BASE_MESH_EDGE_HPP_
#define SRC_MESH_CL_BASE_MESH_EDGE_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src

#include "cl_Communication_Tools.hpp" // COM/src
namespace moris
{

    class Base_Mesh_Edge
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Base_Mesh_Edge()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Base_Mesh_Edge() = default;

        /**
           * Provides the number of edges in x-direction
           *
           * @param[in] aLevel              Level of the element.
           * @param[in] aModelDim                Number of dimensions.
           * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
           *
           * @param[out] tEdgeNumber        Number of edges in x-direction.
           *
           */
          uint
          give_number_of_edges_x(
                  uint const & aLevel,
                  uint const & aModelDim,
                  Mat<uint> const & aNumberOfElementsPerDirection) const;

          /**
            * Provides the number of edges in y-direction
            *
            * @param[in] aLevel              Level of the element.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
            *
            * @param[out] tEdgeNumber        Number of edges in y-direction.
            *
            */
           uint
           give_number_of_edges_y(
                   uint const & aLevel,
                   uint const & aModelDim,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

           /**
             *  Provides the number of edges in z-direction
             *
             * @param[in] aLevel              Level of the element.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] tEdgeNumber        Number of edges in z-direction.
             *
             */
            uint
            give_number_of_edges_z(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the number of edges within all levels from level 0 until aLevel
             *
             * @param[in] aLevel              Level of the basis functions.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] edges_number        Number of edges within all Levels until aLevel.
             *
             */
            uint
            give_number_of_edges(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
              * Provides the edge number in x-direction of the position i,j,k
              *
              * @param[in] aLevel              Level of the element.
              * @param[in] aModelDim                Number of dimensions.
              * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
              * @param[in] aIJKPosition        Position I,J,K.
              *
              * @param[out] tEdgeNumber        Number of edge in x-direction.
              *
              */
            uint
            give_edge_x_of_position(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection,
                    Mat<uint> const & aIJKPosition) const;

            /**
              * Provides the edge number in y-direction of the position i,j,k
              *
              * @param[in] aLevel              Level of the element.
              * @param[in] aModelDim                Number of dimensions.
              * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
              * @param[in] aIJKPosition        Position I,J,K.
              *
              * @param[out] tEdgeNumber        Number of edge in y-direction.
              *
              */
            uint
            give_edge_y_of_position(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection,
                    Mat<uint> const & aIJKPosition) const;

            /**
              * Provides the edge number in z-direction of the position i,j,k
              *
              * @param[in] aLevel              Level of the element.
              * @param[in] aModelDim                Number of dimensions.
              * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
              * @param[in] aIJKPosition        Position I,J,K.
              *
              * @param[out] tEdgeNumber        Number of edge in z-direction.
              *
              */
            uint
            give_edge_z_of_position(
                    uint const & aLevel,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection,
                    Mat<uint> const & aIJKPosition) const;

            /**
             * Provides the level of the edge from the hierarchical mesh.
             *
             * @param[in] aEdgeId            Edge number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] level        Level of the edge from the hierarchical mesh.
             *
             */
            uint
            give_edge_level(
                    uint const & aEdgeId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the edges, which are connected to an element.
             *
             * @param[in] aElement            Element number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] Edges        Edges, which are connected to an element. In 2D: Edges(4,1) = [left edge_X, right edge_X, bottom edge_Y, top_edge_Y], In 2D: Edges(12,1) = [bottom left edge_X, bottom right edge_X, top left edge_X, top right edge_X, front bottom edge_Y, front top_edge_Y, back bottom edge_Y, back top_edge_Y, frontleft edge_z, front right edge_z, back left edge_z, back right edge_z]
             *
             */
            Mat<uint>
            give_element_edges(
                    uint const & aElement,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the elements, which are connected to an edge.
             *
             * @param[in] aEdgeId            Edge Id.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] Element        Elements, which are connected to an edge.
             *
             */
            Mat<uint>
            give_elements_of_edge(
                    uint const & aEdgeId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the position of an edge
             *
             * @param[in] aEdgeId            Edge number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             *
             * @param[out] Position         Position I,J,K.
             *
             */
            Mat<uint>
            give_edge_position(
                    uint const & aEdgeId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection) const;

            /**
             * Provides the ownership of an edge
             *
             * @param[in] aEdgeId            Edge number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
             * @param[in] aProcNeighbors    Vector with a list of proc neighbors
             *
             * @param[out] Rank        Proc rank
             *
             */
            uint
            give_edge_owner(
                    uint const & aEdgeId,
                    uint const & aModelDim,
                    Mat<uint> const & aNumberOfElementsPerDirection,
                    Mat<uint> const & aProcRange,
                    Mat<uint> const & aProcNeighbors) const;

            /**
             * Provides the a list, which share an edge
             *
             * @param[in] aEdgeId            Edge number.
             * @param[in] aModelDim                Number of dimensions.
             * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
             * @param[in] aProcRange        Vector with a range of the proc. 2D: (x_start, x_end, y_start, y_end), 3D: (x_start, x_end, y_start, y_end, z_start, z_end)
             * @param[in] aProcNeighbors    Vector with a list of proc neighbors
             *
             * @param[out] Share        Vector with procs, which share the edge
             *
             */
            Mat<uint>
            give_edge_share(
                    uint const& aEdgeId,
                    uint const& aModelDim,
                    Mat<uint> const& aNumberOfElementsPerDirection,
                    Mat<uint> const & aProcRange,
                    Mat<uint> const & aProcNeighbors) const;
    };
}

#endif /* SRC_MESH_CL_BASE_MESH_EDGE_HPP_ */
