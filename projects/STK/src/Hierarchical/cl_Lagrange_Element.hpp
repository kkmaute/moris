/*
 * cl_Lagrange_Element.hpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_LAGRANGE_ELEMENT_HPP_
#define SRC_MESH_CL_LAGRANGE_ELEMENT_HPP_

#include "cl_Lagrange_Basis.hpp"
#include "cl_Base_Mesh_Element.hpp"

namespace moris
{

    class Lagrange_Element
    {
    protected:

    public:
        //Create Object of base Element
        Base_Mesh_Element mBaseElement;
        /**
         * Hierarchical_Mesh constructor
         */
        Lagrange_Element()
    {
    }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Lagrange_Element() = default;

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

    };
}

#endif /* SRC_MESH_CL_LAGRANGE_ELEMENT_HPP_ */
