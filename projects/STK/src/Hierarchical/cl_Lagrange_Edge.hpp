/*
 * cl_Lagrange_Edge.hpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_LAGRANGE_EDGE_HPP_
#define SRC_MESH_CL_LAGRANGE_EDGE_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "cl_Base_Mesh_Edge.hpp"
#include "cl_Lagrange_Basis.hpp"
#include "cl_Lagrange_Element.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src

#include "cl_Communication_Tools.hpp" // COM/src
namespace moris
{

    class Lagrange_Edge
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        //Create Object of Basis
        Lagrange_Basis mLagrangeBasis;
        //Create Object of lagrange element
        Lagrange_Element mLagrangeElement;
        //Create Object of base edge
        Base_Mesh_Edge mBaseEdge;
        /**
         * Hierarchical_Mesh constructor
         */
        Lagrange_Edge()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Lagrange_Edge() = default;

           /**
            * Provides the nodes of an edge
            *
            * @param[in] aEdgeId            Edge number.
            * @param[in] aModelDim                Number of dimensions.
            * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction].
            *
            * @param[out] Nodes         Node id's, which are connected to an edge (Only for linear possible)
            *
            */
           Mat<uint>
           give_edge_nodes(
                   uint const & aEdgeId,
                   uint const & aModelDim,
                   uint const & aPolynomial,
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

    };
}
#endif /* SRC_MESH_CL_LAGRANGE_EDGE_HPP_ */
