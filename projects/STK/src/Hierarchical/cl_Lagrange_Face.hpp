/*
 * cl_Lagrange_Face.hpp
 *
 *  Created on: Mar 1, 2018
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_LAGRANGE_FACE_HPP_
#define SRC_MESH_CL_LAGRANGE_FACE_HPP_

#include "algorithms.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Base_Mesh_Face.hpp"
#include "cl_Lagrange_Element.hpp"

namespace moris
{

    class Lagrange_Face
    {
    protected:

    public:
        //Create Object of Base Face
        Base_Mesh_Face mBaseFace;
        //Create Object of lagrange element
        Lagrange_Element mLagrangeElement;

        /**
         * Hierarchical_Mesh constructor
         */
        Lagrange_Face()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Lagrange_Face() = default;

        /**
         * Provides the basis, which influence a face
         *
         * @param[in] aFaceId            Face number.
         * @param[in] aModelDim                Number of dimensions.
         * @param[in] aNumberOfElementsPerDirection        NumElements=[Number of elements in x,y,z direction] 1x3.
         *
         * @param[out] Nodes         Node id's, which are connected to a face
         *
         */
        Mat<uint>
        give_face_basis(
                uint const & aFaceId,
                uint const & aModelDim,
                uint const & aPolynomial,
                Mat<uint> const & aNumberOfElementsPerDirection) const;
    };
}

#endif /* SRC_MESH_CL_LAGRANGE_FACE_HPP_ */
