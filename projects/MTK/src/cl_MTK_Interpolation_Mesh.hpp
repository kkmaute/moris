/*
 * cl_MTK_Interpolation_Mesh.hpp
 *
 *  Created on: Apr 15, 2019
 *      Author: doble
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_

#include "cl_MTK_Mesh_Core.hpp"

#include "assert.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace mtk
{
class Interpolation_Mesh: public virtual Mesh
{
    // Functions only valid for interpolation mIntegrationMeshes
public:
    Interpolation_Mesh(){};

    /*
     * Get elements interpolated into by a basis function. For a Lagrange mesh,
     * the elements in support of basis is equivalent to the elements connected
     * to a node. Therefore, a call to get_elements
     */
    virtual
    Matrix< IndexMat >
    get_elements_in_support_of_basis(moris_index aBasisIndex,
                                     moris_index aInterpIndex = 0)
                                     {
        MORIS_ERROR(0,"get_elements_in_support_of_basis not implemented");
        return Matrix<IndexMat>(0,0);
                                     }

};
}
}



#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERPOLATION_MESH_HPP_ */
