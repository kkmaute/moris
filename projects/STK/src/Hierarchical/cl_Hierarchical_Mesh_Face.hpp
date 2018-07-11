/*
 * cl_Hierarchical_Mesh_Face.hpp
 *
 *  Created on: Dec 12, 2017
 *      Author: gleim
 */

#ifndef SRC_MESH_CL_HIERARCHICAL_MESH_FACE_HPP_
#define SRC_MESH_CL_HIERARCHICAL_MESH_FACE_HPP_

#include "algorithms.hpp"
#include "cl_Hierarchical_Mesh_Basis.hpp"
#include "linalg.hpp"
#include "cl_Base_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Base_Mesh_Face.hpp"

namespace moris
{

    class Hierarchical_Mesh_Face
    {
    protected:

    public:
        //Create Object of Base Face
           Base_Mesh_Face mBaseFace;
           //Create Object of B-spline basis
              Hierarchical_Mesh_Basis mBasis;

        /**
         * Hierarchical_Mesh constructor
         */
        Hierarchical_Mesh_Face()
        {
        }

        /**
         * Hierarchical_Mesh destructor.
         */
        ~Hierarchical_Mesh_Face() = default;

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
                   Mat<uint> const & aNumberOfElementsPerDirection) const;
    };
}

#endif /* SRC_MESH_CL_HIERARCHICAL_MESH_FACE_HPP_ */
