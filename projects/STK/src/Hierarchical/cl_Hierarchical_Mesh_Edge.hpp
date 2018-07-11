/*
 * cl_Hierarchical_Mesh_Edge.hpp
 *
 *  Created on: Dec 12, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_EDGE_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_EDGE_HPP_

#include "algorithms.hpp"
#include "cl_Base_Mesh_Element.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Base_Mesh_Edge.hpp"
#include "cl_Communication_Tools.hpp" // COM/src
namespace moris
{

    class Hierarchical_Mesh_Edge
    {
    protected:

    public:
        //Create Object of BaseElement
        Base_Mesh_Element mBaseElement;
        //Create Object of Basis
        Hierarchical_Mesh_Basis mBasis;
        //Create Object of Basis Edge
        Base_Mesh_Edge mBaseEdge;
        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Edge()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Edge() = default;

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
                   Mat<uint> const & aNumberOfElementsPerDirection) const;

    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_EDGE_HPP_ */
